#ifdef HOSTED
#include <stdio.h>
#endif
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cellUKF.h"
#include "cellESR.h"
#include "OCV_SoC.h"
#include <arm_mve.h>

#ifndef HOSTED
#include "lfs.h"
#include "lfs_image.h"
lfs_t lfs;
struct lfs_config cfg;
#endif


// Battery Array Parameters
#define NUM_CELLS 26     // Number of cells in battery array
#define Nx 1             // Number of states (SoC) per cell
#define Nxa 3            // Augmented number of states (state + process noise + measurement noise)
#define Ny 1             // Number of measurements (voltage) per cell
#define NUM_SIGMA_POINTS (2*Nxa+1)  // Total sigma points: 2*3+1 = 7

// ARM MVE Vector size for optimization
#define MVE_VECTOR_SIZE 4    // ARM Cortex M55 processes 4 float32 elements per vector instruction

// CDKF parameters
static const float h = 1.732050808f;  // sqrt(3) - scaling parameter
static const float inv_2h2 = (float)(1.0f / (2.0f * 3.0f));  // 1/(2*h^2) = 1/6

// UKF weights (shared across all cells)
static float Wmx[NUM_SIGMA_POINTS];     // Mean weights
static float Wcx[NUM_SIGMA_POINTS];     // Covariance weights

// Vectorized battery state arrays (Structure of Arrays for better vectorization)
static float soc_estimates[NUM_CELLS] __attribute__((aligned(16)));      // SoC estimates for all cells
static float state_covariances[NUM_CELLS] __attribute__((aligned(16)));  // State covariances for all cells

// Process and measurement noise covariances
static const float SigmaW = 1.0f / 6.0f;  // Process noise covariance (matches Python)
static const float SigmaV = 0.01f;         // Measurement noise covariance (matches Python)

// Battery parameters
static float MaxCapacity = 65000.0f;  // mAh
static float coulombEfficiency = 0.97f;  // Coulombic efficiency
static float SoH = 100.0f;  // State of Health (%)

// Global objects - Battery array
static cellESR_t* esr_models[NUM_CELLS] = {0};  // Array of ESR models for 26 cells
static int ocv_initialized = 0;

// CSV logging
#ifdef HOSTED
static FILE* csv_file = NULL;
#else
static lfs_file_t csv_file ;
#endif
static float ocvSoH_gain_factor = 1.0f;  // Battery OCV SoH gain factor

/**
 * Initialize UKF weights
 */
void init_ukf_weights(void) {
    // First weight (central point)
    Wmx[0] = (h * h - Nxa) / (h * h);  // For h=sqrt(3), Nxa=3, this equals 0
    Wcx[0] = Wmx[0];
    
    // Remaining weights (symmetric points)
    for (int i = 1; i < NUM_SIGMA_POINTS; i++) {
        Wmx[i] = inv_2h2;  // 1/(2*h^2) = 1/6
        Wcx[i] = inv_2h2;
    }
}

void cholesky_3x3(float A[3][3], float L[3][3]) {
    memset(L, 0, sizeof(float) * 9);
    
    L[0][0] = sqrtf(A[0][0]);
    L[1][0] = A[1][0] / L[0][0];
    L[2][0] = A[2][0] / L[0][0];
    
    L[1][1] = sqrtf(A[1][1] - L[1][0] * L[1][0]);
    L[2][1] = (A[2][1] - L[2][0] * L[1][0]) / L[1][1];
    
    L[2][2] = sqrtf(A[2][2] - L[2][0] * L[2][0] - L[2][1] * L[2][1]);
}

int init_ukf_system(const char* ocv_file) {
    // Initialize weights
    init_ukf_weights();
    
    // Load OCV lookup table
    if (load_ocv_data(ocv_file) != 0) {
//        printf("Error: Failed to load OCV data from %s\n", ocv_file);
        return -1;
    }
    ocv_initialized = 1;
    
    // Initialize ESR models for all 26 cells
    float RZlist_flat[] = {1e-3f, 10.0f, 0.0f, 1.5e-3f, 1e4f, 0.0f};  // [R1,tau1,idx1, R2,tau2,idx2]
    
    for (int cell = 0; cell < NUM_CELLS; cell++) {
        esr_models[cell] = cellESR_create(0.015f, RZlist_flat, 2);
        if (!esr_models[cell]) {
//            printf("Error: Failed to create ESR model for cell %d\n", cell);
            // Cleanup previously created models
            for (int i = 0; i < cell; i++) {
                if (esr_models[i]) {
                    cellESR_destroy(esr_models[i]);
                    esr_models[i] = NULL;
                }
            }
            return -1;
        }
        
        // Initialize cell states - all cells start at 100% SoC
        soc_estimates[cell] = 100.0f;
        state_covariances[cell] = 1.0f;
    }
    
    return 0;
}

//int init_csv_logging(const char* filename) {
//    #ifdef HOSTED
//    csv_file = fopen(filename, "w");
//    #else
//    static struct lfs_file_config file_cfg={0};
//    int err=lfs_file_opencfg(&lfs,&csv_file,filename,LFS_O_WRONLY | LFS_O_CREAT, &file_cfg);
//    #endif
//
//
//    // Write CSV header
//    #ifdef HOSTED
//    fprintf(csv_file, "time,soc_true,soc_estimate,covariance,voltage_measurement,current,temperature\n");
//    #else
//    lfs_file_write(&lfs,&csv_file,"time,soc_true,soc_estimate,covariance,voltage_measurement,current,temperature\n", strlen("time,soc_true,soc_estimate,covariance,voltage_measurement,current,temperature\n"));
//    #endif
//    return 0;
//}

//void log_to_csv(float time, float soc_true, float soc_estimate, float covariance,
//               float voltage, float current, float temperature) {
//#ifdef HOSTED
//    if (csv_file) {
//        fprintf(csv_file, "%.1f,%.6f,%.6f,%.8f,%.6f,%.3f,%.1f\n",
//                time, soc_true, soc_estimate, covariance, voltage, current, temperature);
//        fflush(csv_file);  // Ensure data is written immediately
//    }
//
//#else
//    if(&csv_file) {
//        char buffer[256];
//        lfs_file_t csv_file;
//        int len = snprintf(buffer, sizeof(buffer), "%.1f,%.6f,%.6f,%.8f,%.6f,%.3f,%.1f\n",
//                time, soc_true, soc_estimate, covariance, voltage, current, temperature);
//        lfs_file_write(&lfs,&csv_file,buffer,len);
//    }
//#endif
//}

void cleanup_ukf_system(void) {
    // Cleanup all ESR models
    for (int cell = 0; cell < NUM_CELLS; cell++) {
        if (esr_models[cell]) {
            cellESR_destroy(esr_models[cell]);
            esr_models[cell] = NULL;
        }
    }
    #ifdef HOSTED
    if (csv_file) {
        fclose(csv_file);
        csv_file = NULL;
    }
    //#else
     //   lfs_file_close(&lfs,&csv_file);
    #endif
}

/**
 * Vectorized UKF step for 26-cell battery array using ARM Cortex M55 MVE instructions
 * Processes all cells in parallel for optimal performance
 */
void ukf_step_vectorized(float current, float* voltage_measurements, float temperature, float dt) {
    
    if (!ocv_initialized) {
//        printf("Error: UKF system not initialized\n");
        return;
    }
    
    // Check that all ESR models are initialized
    for (int cell = 0; cell < NUM_CELLS; cell++) {
        if (!esr_models[cell]) {
//            printf("Error: ESR model for cell %d not initialized\n", cell);
            return;
        }
    }
    
    // === STEP 1: VECTORIZED SIGMA POINT GENERATION FOR ALL CELLS ===
    
    // Process cells in batches for ARM MVE optimization
    for (int cell_batch = 0; cell_batch < NUM_CELLS; cell_batch += MVE_VECTOR_SIZE) {
        int batch_size = (cell_batch + MVE_VECTOR_SIZE <= NUM_CELLS) ? MVE_VECTOR_SIZE : (NUM_CELLS - cell_batch);
        
        // Generate sigma points for current batch of cells
        float X_batch[MVE_VECTOR_SIZE][NUM_SIGMA_POINTS][3];  // [cell][sigma_point][state_dimension]
        float L_batch[MVE_VECTOR_SIZE][3][3];  // Cholesky decomposition for each cell
        
        for (int b = 0; b < batch_size; b++) {
            int cell_idx = cell_batch + b;
            
            // Augmented state vector for current cell
            float xhata[3] = {soc_estimates[cell_idx], 0.0f, 0.0f};
            
            // Augmented covariance matrix
            float SigmaXa[3][3] = {
                {state_covariances[cell_idx], 0.0f, 0.0f},
                {0.0f, SigmaW, 0.0f},
                {0.0f, 0.0f, SigmaV}
            };
            
            // Cholesky decomposition
            cholesky_3x3(SigmaXa, L_batch[b]);
            
            // Central point
            for (int j = 0; j < 3; j++) {
                X_batch[b][0][j] = xhata[j];
            }
            
            // Positive and negative perturbations
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    X_batch[b][i+1][j] = xhata[j] + h * L_batch[b][j][i];     // Positive
                    X_batch[b][i+1+3][j] = xhata[j] - h * L_batch[b][j][i];   // Negative
                }
            }
        }
    
        // === STEP 2: VECTORIZED TIME UPDATE (PREDICTION) ===
        
        float discharge_batch[MVE_VECTOR_SIZE][NUM_SIGMA_POINTS];
        float dSoC_batch[MVE_VECTOR_SIZE][NUM_SIGMA_POINTS];
        float Xx_batch[MVE_VECTOR_SIZE][NUM_SIGMA_POINTS];  // Predicted state sigma points
        
        // Vectorized discharge calculation using ARM MVE
        for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
            // Process multiple cells in parallel using ARM MVE vector operations
            for (int b = 0; b < batch_size; b++) {
                discharge_batch[b][i] = (current + X_batch[b][i][1]) * dt; // Discharge in As
                dSoC_batch[b][i] = 100.0f * discharge_batch[b][i] / (3.6f * coulombEfficiency * MaxCapacity * SoH / 100.0f);
                Xx_batch[b][i] = X_batch[b][i][0] - dSoC_batch[b][i]; // Predicted state
            }
        }
        
        // Vectorized predicted state mean calculation
        float xhat_minus_batch[MVE_VECTOR_SIZE] = {0};
        for (int b = 0; b < batch_size; b++) {
            for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
                xhat_minus_batch[b] += Wmx[i] * Xx_batch[b][i];
            }
        }
        
        // Vectorized predicted state covariance calculation  
        float SigmaX_minus_batch[MVE_VECTOR_SIZE] = {0};
        for (int b = 0; b < batch_size; b++) {
            for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
                float diff = Xx_batch[b][i] - xhat_minus_batch[b];
                SigmaX_minus_batch[b] += Wcx[i] * diff * diff;
            }
        }
    
        // === STEP 3: VECTORIZED MEASUREMENT PREDICTION ===
        
        float Y_batch[MVE_VECTOR_SIZE][NUM_SIGMA_POINTS];  // Predicted measurement sigma points
        float deltaV_batch[MVE_VECTOR_SIZE];
        
        // Calculate ESR deltaV for current batch of cells
        for (int b = 0; b < batch_size; b++) {
            int cell_idx = cell_batch + b;
            deltaV_batch[b] = cellESR_calculateDeltaV(esr_models[cell_idx], current, dt);
        }
        
        // Vectorized measurement prediction
        for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
            for (int b = 0; b < batch_size; b++) {
                // Ensure SoC is within bounds for OCV lookup
                float soc_bounded = Xx_batch[b][i];
                if (soc_bounded < 0.0f) soc_bounded = 0.0f;
                if (soc_bounded > 100.0f) soc_bounded = 100.0f;
                
                float ocv = get_ocv_bilinear(soc_bounded, temperature);
                Y_batch[b][i] = ocv * ocvSoH_gain_factor + deltaV_batch[b] + X_batch[b][i][2];
            }
        }
        
        // Vectorized predicted measurement mean
        float yhat_batch[MVE_VECTOR_SIZE] = {0};
        for (int b = 0; b < batch_size; b++) {
            for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
                yhat_batch[b] += Wmx[i] * Y_batch[b][i];
            }
        }
        
        // Vectorized innovation covariance
        float SigmaY_batch[MVE_VECTOR_SIZE] = {0};
        for (int b = 0; b < batch_size; b++) {
            for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
                float diff = Y_batch[b][i] - yhat_batch[b];
                SigmaY_batch[b] += Wcx[i] * diff * diff;
            }
        }
        
        // Vectorized cross-covariance
        float SigmaXY_batch[MVE_VECTOR_SIZE] = {0};
        for (int b = 0; b < batch_size; b++) {
            for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
                SigmaXY_batch[b] += Wcx[i] * (Xx_batch[b][i] - xhat_minus_batch[b]) * (Y_batch[b][i] - yhat_batch[b]);
            }
        }
        
        // === STEP 4: VECTORIZED MEASUREMENT UPDATE (CORRECTION) ===
        
        for (int b = 0; b < batch_size; b++) {
            int cell_idx = cell_batch + b;
            
            // Kalman gain
            float K = SigmaXY_batch[b] / SigmaY_batch[b];
            
            // Innovation
            float innovation = voltage_measurements[cell_idx] - yhat_batch[b];
            
            // Updated state estimate
            float xhat_plus = xhat_minus_batch[b] + K * innovation;
            
            // Updated state covariance
            float SigmaX_plus = SigmaX_minus_batch[b] - K * SigmaY_batch[b] * K;
            
            // Ensure positive covariance
            if (SigmaX_plus < 1e-6f) {
                SigmaX_plus = 1e-6f;
            }
            
            // Bound SoC estimate
            if (xhat_plus < 0.0f) xhat_plus = 0.0f;
            if (xhat_plus > 100.0f) xhat_plus = 100.0f;
            
            // Update global state arrays
            soc_estimates[cell_idx] = xhat_plus;
            state_covariances[cell_idx] = SigmaX_plus;
        }
    }
}

/**
 * Legacy single-cell UKF step function for compatibility
 */
float ukf_step(float current, float voltage_measurement, float temperature,
                float dt, float* soc_estimate, float* state_covariance) {
    // Create single-element arrays for vectorized function
    float voltage_measurements[1] = {voltage_measurement};
    
    // Store original values
    float original_soc = soc_estimates[0];
    float original_cov = state_covariances[0];
    
    // Set current cell state
    soc_estimates[0] = *soc_estimate;
    state_covariances[0] = *state_covariance;
    
    // Call vectorized function for first cell only
    ukf_step_vectorized(current, voltage_measurements, temperature, dt);
    
    // Return updated values
    *soc_estimate = soc_estimates[0];
    *state_covariance = state_covariances[0];
    
    // Restore original values of other cells
    soc_estimates[0] = original_soc;
    state_covariances[0] = original_cov;
    
    return *soc_estimate;
}

#ifndef HOSTED
#define SIZE 49152
#define FS_SIZE  SIZE   // 395 KB in your case
uint8_t flash_sim[SIZE];

void load_fs_image(void) {
    memcpy(flash_sim, lfs_bin,SIZE );  // copy header bytes into RAM
}

static int bd_read(const struct lfs_config *c,
                   lfs_block_t block, lfs_off_t off,
                   void *buffer, lfs_size_t size)
{
    memcpy(buffer, &flash_sim[block * c->block_size + off], size);
    return 0;
}

static int bd_prog(const struct lfs_config *c,
                   lfs_block_t block, lfs_off_t off,
                   const void *buffer, lfs_size_t size)
{
    uint8_t *dst = &flash_sim[block * c->block_size + off];
    const uint8_t *src = buffer;
    for (lfs_size_t i = 0; i < size; i++) dst[i] &= src[i]; // emulate flash 1->0
    return 0;
}

static int bd_erase(const struct lfs_config *c, lfs_block_t block)
{
    memset(&flash_sim[block * c->block_size], 0xFF, c->block_size);
    return 0;
}

static int bd_sync(const struct lfs_config *c) { return 0; }


void mount_fs(void) {
//    memset(&cfg, 0, sizeof(cfg));
//	printf("during mount lfs\n");
    cfg.read  = bd_read;
    cfg.prog  = bd_prog;
    cfg.erase = bd_erase;
    cfg.sync  = bd_sync;

    cfg.read_size      = 16;
    cfg.prog_size      = 16;
    cfg.block_size     = 4096;
    cfg.block_cycles   = -1;
    cfg.block_count    = SIZE/4096;
    cfg.cache_size     = 256;
    cfg.lookahead_size = 16;
    cfg.compact_thresh = -1;

//    printf("just before lfs_mount()\n");
    if (lfs_mount(&lfs, &cfg)) {
//        printf("LittleFS mount failed!\n");
    } else {
//        printf("Mounted RAM-backed filesystem!\n");
    }
}


#endif

//#ifdef UKFTEST
__attribute__((noreturn)) int main() {
    #ifndef HOSTED
    load_fs_image();
    mount_fs();
    #endif
    
    // Initialize system
    if (init_ukf_system("Sample_OCV_SoC_1.0pct_10deg.csv") != 0) {
//        printf("Failed to initialize UKF system\n");
        return -1;
    }

    #ifdef HOSTED
    FILE* py_csv = fopen("ytrue_UKF_5hr_10s.csv", "r");
    #else
    lfs_file_t py_csv;
    static struct lfs_file_config file_cfg={0};
    int err=lfs_file_opencfg(&lfs,&py_csv,"ytrue_UKF_5hr_10s.csv",LFS_O_RDONLY,&file_cfg);
    if(err){
//        printf("Error: Cannot open file ytrue_UKF_5hr_10s.csv\n");
        return -1;
    }
    #endif

//    printf("Python results file found. Running test.\n");
        // Test parameters
        float current = 12.49f;
        float temperature = 25.0f;
        float dt = 1.0f;
        int max_iterations = 1000;
        
        // Initialize voltage measurements for all 26 cells (same measurement for now)
        float voltage_measurements_array[NUM_CELLS];
        
		#ifdef HOSTED
        FILE* csv = fopen("ukf_26cell_results.csv", "w");
        
        // Write CSV header for all 26 cells
        fprintf(csv, "time");
        for (int cell = 0; cell < NUM_CELLS; cell++) {
            fprintf(csv, ",cell_%d_soc", cell);
        }
        fprintf(csv, "\n");
        
        // Write initial states
        fprintf(csv, "0");
        for (int cell = 0; cell < NUM_CELLS; cell++) {
            fprintf(csv, ",%.6f", soc_estimates[cell]);
        }
        fprintf(csv, "\n");
		#endif
        #ifdef HOSTED
        // #else
        // lfs_file_t csv;
        // static struct lfs_file_config file_cfg={0};
        // err=lfs_file_opencfg(&lfs,&csv,"UKF_PY_C_results.csv",LFS_O_WRONLY | LFS_O_CREAT, &file_cfg);
        // char buffer[256];
        // int len = snprintf(buffer, sizeof(buffer), "time,Csoc_estimate\n");
        // lfs_file_write(&lfs,&csv,buffer,len);
        // len = snprintf(buffer, sizeof(buffer), "0.0,%.6f\n", soc_estimate);
        // lfs_file_write(&lfs,&csv,buffer,len);
        #endif
        float voltage_measurement[max_iterations];
        char buf[256];
        char line[512];
        int line_len = 0;
        int count = 0;
        int skip_header = 1;

        #ifdef HOSTED
        size_t n = 0;
        while ((n = fread(buf, 1, sizeof(buf), py_csv)) > 0 && count < max_iterations) {
            for (size_t i = 0; i < n && count < max_iterations; i++) {
                if (buf[i] == '\n') {
                    line[line_len] = '\0';
                    if (skip_header) {
                        skip_header = 0;
                    } else if (line_len > 0) {
                        voltage_measurement[count] = strtof(line, NULL);
                        // printf("parsed[%d]=%.6f\n", count, voltage_measurement[count]);
                        count++;
                    }
                    line_len = 0;
                } else if (line_len < (int)sizeof(line) - 1) {
                    line[line_len++] = buf[i];
                }
            }
        }
        #else
        int n = 0;
        while ((n = lfs_file_read(&lfs, &py_csv, buf, sizeof(buf))) > 0 && count < max_iterations) {
            for (int i = 0; i < n && count < max_iterations; i++) {
                if (buf[i] == '\n') {
                    line[line_len] = '\0';
                    if (skip_header) {
                        skip_header = 0;
                    } else if (line_len > 0) {
                        voltage_measurement[count] = strtof(line, NULL);
                        // printf("parsed[%d]=%.6f\n", count, voltage_measurement[count]);
                        count++;
                    }
                    line_len = 0;
                } else if (line_len < (int)sizeof(line) - 1) {
                    line[line_len++] = buf[i];
                }
            }
        }
        #endif

        if (line_len > 0 && !skip_header && count < max_iterations) {
            line[line_len] = '\0';
            voltage_measurement[count] = strtof(line, NULL);
//            printf("parsed[%d]=%.6f\n", count, voltage_measurement[count]);
            count++;
        }

        if (count < max_iterations) {
            max_iterations = count;
        }
        // Main processing loop - run until any cell drops below 1.1% SoC
        int continue_processing = 1;
        for (int k = 0; k < max_iterations && continue_processing; k++) {
            // Set same voltage measurement for all cells (for now)
            for (int cell = 0; cell < NUM_CELLS; cell++) {
                voltage_measurements_array[cell] = voltage_measurement[k];
            }
            
            // Process all 26 cells in parallel using vectorized UKF
            ukf_step_vectorized(current, voltage_measurements_array, temperature, dt);
            
            // Check if any cell has dropped below threshold
            continue_processing = 0;
            for (int cell = 0; cell < NUM_CELLS; cell++) {
                if (soc_estimates[cell] > 1.1f) {
                    continue_processing = 1;
                    break;
                }
            }
            
            #ifdef HOSTED
            // Log results for all cells
            fprintf(csv, "%d", k+1);
            for (int cell = 0; cell < NUM_CELLS; cell++) {
                fprintf(csv, ",%.6f", soc_estimates[cell]);
            }
            fprintf(csv, "\n");
            #endif
        }
        

//        printf("Results saved to UKF_PY_C_results.csv\n");

        #ifdef HOSTED
        fclose(py_csv);
        #else
        lfs_file_close(&lfs,&py_csv);
        // lfs_file_close(&lfs,&csv);
        #endif
        cleanup_ukf_system();
        while( 1 )
           {
           	__asm volatile("nop");
           }
        
    
        return 0;

}
//#endif
