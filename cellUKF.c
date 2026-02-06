#include "y_true.h"
#include "ocv_table.h"
#include "cellESR.h"

#ifdef DEBUGG
#include <stdio.h>
#endif

float sqrtf(float x) {
    if (x < 0.0f) return 0.0f / 0.0f; // NaN
    if (x == 0.0f) return 0.0f;

    float guess = x;
    for (int i = 0; i < 10; i++) {
        guess = 0.5f * (guess + x / guess);
    }
    return guess;
}

float fabs_f(float x) {
    return x < 0.0f ? -x : x;
}

// max number of RZ elements for embedded systems
#define MAX_RZ_ELEMENTS 20

// Structure for cellESR class
struct cellESR {
    // basic params
    float R0;
    float R0multiplier;
    float esrDC;
    float deltaV;
    int time;
    
    // arrays for RZ elements
    float R_values[MAX_RZ_ELEMENTS];
    float tau_values[MAX_RZ_ELEMENTS];
    int indices[MAX_RZ_ELEMENTS];
    float ix_values[MAX_RZ_ELEMENTS];
    int num_elements;
    
};

// for shared library
#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

// Create and initialize cellESR structure
cellESR_t cell[sizeof(cellESR_t)] = {0};
EXPORT cellESR_t* cellESR_create(float R0, float* RZlist_flat, int num_RZ) {
    // cellESR_t* cell= (cellESR_t*)(malloc ? malloc(sizeof(cellESR_t)) : cell);
    // if (!cell) return NULL;
    
    cell->R0 = R0;
    cell->R0multiplier = 1.0f;
    cell->esrDC = R0;
    cell->deltaV = 0.0f;
    cell->time = 0;
    cell->num_elements = (num_RZ < MAX_RZ_ELEMENTS) ? num_RZ : MAX_RZ_ELEMENTS;
    
    // Process RZlist (flattened array: [R0, Z0, type0, R1, Z1, type1, ...])
    for (int i = 0; i < cell->num_elements; i++) {
        float R = RZlist_flat[i * 3];
        float Z = RZlist_flat[i * 3 + 1];
        int type = (int)RZlist_flat[i * 3 + 2];  // 0 or 1
        
        cell->R_values[i] = R;
        cell->ix_values[i] = 0.0f;
        cell->indices[i] = type;
        
        if (type == 0) {
            // tau = R * Z, and add R to esrDC
            cell->tau_values[i] = R * Z;
            cell->esrDC += R;
        } else {
            // tau = Z / R
            cell->tau_values[i] = Z / R;
        }
    }
    
    return cell;
}

// Destroy cellESR structure
EXPORT void cellESR_destroy(cellESR_t* cell) {
    if (cell) {
        // free(cell);
    }
}

// Calculate delta V using ESR model
EXPORT float cellESR_calculateDeltaV(cellESR_t* cell, float i_new, float dt) {
    if (!cell) return 0.0f;
    
    // Update time
    cell->time += dt;
    
    // Calculate ESR contribution
    float ESRout = cell->R0 * cell->R0multiplier;
    float dV = i_new * ESRout;
    
    // Process each RZ element
    for (int i = 0; i < cell->num_elements; i++) {
        // Update ix_values using: ix_next = (ix_previous + i_new * dt / tau) / (1 + dt / tau)
        cell->ix_values[i] = (cell->ix_values[i] + i_new * dt / cell->tau_values[i]) / 
                             (1.0f + dt / cell->tau_values[i]);
        
        // Clamp values to prevent overshoot
        if (fabs_f(cell->ix_values[i]) > fabs_f(i_new)) {
            cell->ix_values[i] = i_new;
        }
        
        // Calculate deltaV contribution based on type
        if (cell->indices[i] == 1) {
            // For indices == 1: dV += R * (i_new - ix_values)
            dV += cell->R_values[i] * (i_new - cell->ix_values[i]);
        } else {
            // For indices == 0: dV += R * ix_values  
            dV += cell->R_values[i] * cell->ix_values[i];
        }
    }
    
    cell->deltaV = dV;
    return dV;
}

// Getter functions for accessing private data
EXPORT float cellESR_get_R0(cellESR_t* cell) {
    return cell ? cell->R0 : 0.0f;
}

EXPORT float cellESR_get_esrDC(cellESR_t* cell) {
    return cell ? cell->esrDC : 0.0f;
}

EXPORT int cellESR_get_time(cellESR_t* cell) {
    return cell ? cell->time : 0;
}

EXPORT float cellESR_get_deltaV(cellESR_t* cell) {
    return cell ? cell->deltaV : 0.0f;
}

// Additional getter functions for arrays
EXPORT int cellESR_get_num_elements(cellESR_t* cell) {
    return cell ? cell->num_elements : 0;
}

EXPORT void cellESR_get_R_values(cellESR_t* cell, float* output) {
    if (cell && output) {
        for (int i = 0; i < cell->num_elements; i++) {
            output[i] = cell->R_values[i];
        }
    }
}

EXPORT void cellESR_get_tau_values(cellESR_t* cell, float* output) {
    if (cell && output) {
        for (int i = 0; i < cell->num_elements; i++) {
            output[i] = cell->tau_values[i];
        }
    }
}

EXPORT void cellESR_get_indices(cellESR_t* cell, int* output) {
    if (cell && output) {
        for (int i = 0; i < cell->num_elements; i++) {
            output[i] = cell->indices[i];
        }
    }
}

EXPORT void cellESR_get_ix_values(cellESR_t* cell, float* output) {
    if (cell && output) {
        for (int i = 0; i < cell->num_elements; i++) {
            output[i] = cell->ix_values[i];
        }
    }
}



typedef struct{
	float soc;
	float temp;
	float val;
}data;

// Cache for last lookup cell
static struct {
    data luc;
    data ldc;
    data ruc;
    data rdc;
} last_ocv = {{-1,0.0f},{-1, 0.0f},{-1, 0.0f},{-1, 0.0f}};  // Initialize to impossible values

// Temperature array: -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
static const float temperatures[TEMP_POINTS] = {
    -40.0f, -30.0f, -20.0f, -10.0f, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f
};
// float get_ocv_bilinear(float soc, float temp) {
//     if (soc < 0.0f) soc = 0.0f;
//     if (soc > 100.0f) soc = 100.0f;
//     if (temp < -40.0f) temp = -40.0f;
//     if (temp > 90.0f) temp = 90.0f;

//     // Use int soc as index
//     int soc_lower = (int)soc;
//     int temp_index = (int)(temp*0.1f + 4.0f);
//     float soc_weight = 0.0f;
//     float temp_weight = 0.0f;
//     printf("%.4f %.4f",soc,temp);
//     if(soc>=last_ocv.luc.soc && soc<=last_ocv.ldc.soc && temp<=last_ocv.luc.temp && temp>=last_ocv.ruc.temp){
//         soc_weight = soc - last_ocv.luc.soc;
//         temp_weight = (temp - last_ocv.luc.temp)*0.1; // 1/(temp_at_upper - temp_at_lower)
//         // return (last_ocv.luc.val*(1.0f-temp_weight)+last_ocv.ruc.val*temp_weight)*(1.0f-soc_weight)+\
//         			(last_ocv.ldc.val*(1.0f-temp_weight)+last_ocv.rdc.val*temp_weight)*soc_weight;
//         float q11 = last_ocv.luc.val;
//         float q12 = last_ocv.ruc.val;
//         float q21 = last_ocv.ldc.val;
//         float q22 = last_ocv.rdc.val;
//         float r1  = q11*(1.0f-temp_weight)+q12*temp_weight;
//         float r2  = q21*(1.0f-temp_weight)+q22*temp_weight;
//         printf("--hi--");
//         return r1*(1.0f-soc_weight)+r2*soc_weight; 

//     }

//     int soc_upper = soc_lower;

//     // Handle SoC corner cases
//     if (soc_lower < 100) {
//         soc_upper = soc_lower + 1;
//         soc_weight = soc - (float)soc_lower;
//     }
//     int temp_lower = temp_index;
//     int temp_upper = temp_index;

//     // Handle temperature corner cases
//     if (temp_index < 0) {
//         // temp < -40, already clamped above
//         temp_lower = temp_upper = 0;
//     } else if (temp_index >= TEMP_POINTS) {
//         // temp > 90, already clamped above
//         temp_lower = temp_upper = TEMP_POINTS - 1;
//     } else if (temp_index < TEMP_POINTS - 1) {
//         // Normal case: interpolate between two temperature points
//         temp_upper = temp_index + 1;
//         float temp_at_lower = temperatures[temp_lower];
//         float temp_at_upper = temperatures[temp_upper];
//         temp_weight = (temp - temp_at_lower) *0.1f; // 1/(temp_at_upper - temp_at_lower)
//     }
//     float q11 = ocv_table[soc_lower][temp_lower];
//     float q12 = ocv_table[soc_lower][temp_upper];
//     float q21 = ocv_table[soc_upper][temp_lower];
//     float q22 = ocv_table[soc_upper][temp_upper];
//     last_ocv.luc.soc=soc_lower;
//     last_ocv.luc.val=q11;
//     last_ocv.ldc.soc=soc_upper;
//     last_ocv.ldc.val=q12;
//     last_ocv.luc.temp=temperatures[temp_lower];
//     last_ocv.ruc.val=q21;
//     last_ocv.ruc.temp=temperatures[temp_upper];
//     last_ocv.rdc.val=q22;

//     // Perform bilinear interpolation
//     float temp_weight_inv = 1.0f - temp_weight;
//     float soc_weight_inv = 1.0f - soc_weight;

//     float r1 = q11 * temp_weight_inv + q12 * temp_weight;
//     float r2 = q21 * temp_weight_inv + q22 * temp_weight;

//     return r1 * soc_weight_inv + r2 * soc_weight;
// }

float get_ocv_bilinear(float soc, float temp) {
	// int x=get_cycle_count();
    // if (!table_initialized) {
    	// bilinearCounts+=0;
        // return -1.0f;
    // }

    // Clamp inputs for corner cases
    if (soc < 0.0f) soc = 0.0f;
    if (soc > 100.0f) soc = 100.0f;
    if (temp < -40.0f) temp = -40.0f;
    if (temp > 90.0f) temp = 90.0f;

    // Use int soc as index
    int soc_lower = (int)soc;
    int soc_upper = soc_lower;
    float soc_weight = 0.0f;

    // Handle SoC corner cases
    if (soc_lower < 100) {
        soc_upper = soc_lower + 1;
        soc_weight = soc - (float)soc_lower;
    }
    // If soc_lower == 100, soc_upper stays 100 and soc_weight stays 0

    // Use int (temp+40)/10 as temp index
    int temp_index = (int)(temp + 40.0f) / 10;
    int temp_lower = temp_index;
    int temp_upper = temp_index;
    float temp_weight = 0.0f;

    // Handle temperature corner cases
    if (temp_index < 0) {
        // temp < -40, already clamped above
        temp_lower = temp_upper = 0;
    } else if (temp_index >= TEMP_POINTS) {
        // temp > 90, already clamped above
        temp_lower = temp_upper = TEMP_POINTS - 1;
    } else if (temp_index < TEMP_POINTS - 1) {
        // Normal case: interpolate between two temperature points
        temp_upper = temp_index + 1;
        float temp_at_lower = temperatures[temp_lower];
        float temp_at_upper = temperatures[temp_upper];
        temp_weight = (temp - temp_at_lower) *0.1f; // 1/(temp_at_upper - temp_at_lower)
    }
    // If temp_index == TEMP_POINTS-1, temp_upper stays same as temp_lower

    // Get the four corner points for bilinear interpolation
    float q11 = ocv_table[100-soc_lower][temp_lower];
    float q12 = ocv_table[100-soc_lower][temp_upper];
    float q21 = ocv_table[100-soc_upper][temp_lower];
    float q22 = ocv_table[100-soc_upper][temp_upper];
    // printf("socl--%d socu--%d templ--%d tempu--%d\n",soc_lower,soc_upper,temp_lower,temp_upper);
    // printf("q11-%.4f q12-%.4f q21-%.4f q22-%.4f ",q11,q12,q21,q22);

    // Perform bilinear interpolation
    float temp_weight_inv = 1.0f - temp_weight;
    float soc_weight_inv = 1.0f - soc_weight;

    float r1 = q11 * temp_weight_inv + q12 * temp_weight;
    float r2 = q21 * temp_weight_inv + q22 * temp_weight;
    // int y=get_cycle_count();
    // bilinearCounts+=y-x;
    return r1 * soc_weight_inv + r2 * soc_weight;
}

#define Nx 1       
#define Nxa 3      
#define Ny 1
#define NUM_SIGMA_POINTS (2*Nxa+1) 

static const float inv_2h2 = (float)(1.0f / (2.0f * 3.0f));  
static const float h = 1.732050808f; 
static const float ocvSoH_gain_factor = 1.0f;

// UKF weights
static float Wmx[NUM_SIGMA_POINTS];     // Mean weights
static float Wcx[NUM_SIGMA_POINTS];     // Covariance weights

// Process and measurement noise covariances
static const float SigmaW = 1.0f / 6.0f;  // Process noise covariance (matches Python)
static const float SigmaV = 0.01f;         // Measurement noise covariance (matches Python)

// Battery parameters
static float MaxCapacity = 65000.0f;  // mAh
static float coulombEfficiency = 0.97f;  // Coulombic efficiency
static float SoH = 100.0f;  // State of Health (%)

// Global objects
static cellESR_t* esr_model ;
static int ocv_initialized = 0;


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
    // memset(L, 0, sizeof(float) * 9);
    L[0][0] = sqrtf(A[0][0]);
    L[1][0] = A[1][0] / L[0][0];
    L[2][0] = A[2][0] / L[0][0];

    L[1][1] = sqrtf(A[1][1] - L[1][0] * L[1][0]);
    L[2][1] = (A[2][1] - L[2][0] * L[1][0]) / L[1][1];

    L[2][2] = sqrtf(A[2][2] - L[2][0] * L[2][0] - L[2][1] * L[2][1]);
}

int init_ukf_system(void) {
    // Initialize weights
    init_ukf_weights();
    ocv_initialized = 1;

    // Initialize ESR model
    float RZlist_flat[] = {1e-3f, 10.0f, 0.0f, 1.5e-3f, 1e4f, 0.0f};  // [R1,tau1,idx1, R2,tau2,idx2]
    esr_model = cellESR_create(0.015f, RZlist_flat, 2);
    if (!esr_model) {
//        printf("Error: Failed to create ESR model\n");
        return -1;
    }

    return 0;
}
void cleanup_ukf_system(void) {
    if (esr_model) {
        cellESR_destroy(esr_model);
        // esr_model = NULL;
    }
}

float ukf_step(float current, float voltage_measurement, float temperature,
                float dt, float* soc_estimate, float* state_covariance) {
	// int x=get_cycle_count();

    if (!ocv_initialized || !esr_model) {
        return *soc_estimate;
    }

    // Current estimates
    float xhat = *soc_estimate;
    float SigmaX = *state_covariance;

    // === STEP 1: GENERATE SIGMA POINTS ===

    float xhata[3] = {xhat, 0.0, 0.0};

    // Augmented covariance matrix
    float SigmaXa[3][3] = {
        {SigmaX, 0.0f, 0.0f},
        {0.0f, SigmaW, 0.0f},
        {0.0f, 0.0f, SigmaV}
    };

    // Cholesky decomposition
    float L[3][3];
    cholesky_3x3(SigmaXa, L);

    // Generate sigma points
    float X[NUM_SIGMA_POINTS][3];  // [sigma_point][state_dimension]

    // Central point
    for (int j = 0; j < 3; j++) {
        X[0][j] = xhata[j];
    }

    // Positive and negative perturbations
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            X[i+1][j] = xhata[j] + h * L[j][i];     // Positive perturbation
            X[i+1+3][j] = xhata[j] - h * L[j][i];   // Negative perturbation
        }
    }

    // === STEP 2: TIME UPDATE (PREDICTION) ===

    float discharge[NUM_SIGMA_POINTS];
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        discharge[i] = (current + X[i][1]) * dt; // Discharge in As
    }

    float dSoC[NUM_SIGMA_POINTS];
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        dSoC[i] = 100.0f * discharge[i] / (3.6f * coulombEfficiency * MaxCapacity * SoH / 100.0f); // Change in SoC considering SoH
    }

    float Xx[NUM_SIGMA_POINTS];  // Predicted state sigma points
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        Xx[i] = X[i][0] - dSoC[i] ;  // state - discharge + process_noise
        // printf("Xx[%d]=%.6f ",i,Xx[i]);
    }
    // printf("\n");


    // Predicted state mean
    float xhat_minus = 0.0f;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        xhat_minus += Wmx[i] * Xx[i];
    }

    // Predicted state covariance
    float SigmaX_minus = 0.0f;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        float diff = Xx[i] - xhat_minus;
        SigmaX_minus += Wcx[i] * diff * diff;
    }

    // === STEP 3: MEASUREMENT PREDICTION ===

    float Y[NUM_SIGMA_POINTS];  // Predicted measurement sigma points

    float deltaV = cellESR_calculateDeltaV(esr_model, current, dt);

    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        // Ensure SoC is within bounds for OCV lookup
        float soc_bounded = Xx[i];
        if (soc_bounded < 0.0f) soc_bounded = 0.0f;
        if (soc_bounded > 100.0f) soc_bounded = 100.0f;

        float ocv = get_ocv_bilinear(soc_bounded, temperature);
        // printf("soc-%.4f temp-%.4f\n",soc_bounded,temperature);

        Y[i] = ocv * ocvSoH_gain_factor + deltaV + X[i][2];
    }

    // Predicted measurement mean
    float yhat = 0.0;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        yhat += Wmx[i] * Y[i];
    }
    // printf("--yhat--%.4f-- ",yhat);

    // Innovation covariance
    float SigmaY = 0.0f;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        float diff = Y[i] - yhat;
        SigmaY += Wcx[i] * diff * diff;
    }

    // Cross-covariance
    float SigmaXY = 0.0f;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        SigmaXY += Wcx[i] * (Xx[i] - xhat_minus) * (Y[i] - yhat);
    }

    // === STEP 4: MEASUREMENT UPDATE (CORRECTION) ===

    // Kalman gain
    float K = SigmaXY / SigmaY;

    // Innovation
    float innovation = voltage_measurement - yhat;

    // Updated state estimate
    float xhat_plus = xhat_minus + K * innovation;

    // Updated state covariance
    float SigmaX_plus = SigmaX_minus - K * SigmaY * K;

    // Ensure positive covariance
    if (SigmaX_plus < 1e-6f) {
        SigmaX_plus = 1e-6f;
    }

    // Bound SoC estimate
    if (xhat_plus < 0.0f) xhat_plus = 0.0f;
    if (xhat_plus > 100.0f) xhat_plus = 100.0f;

    // Update outputs
    *soc_estimate = xhat_plus;
    *state_covariance = SigmaX_plus;
    // int y=get_cycle_count();
    // UKFcounts+=y-x;
    return xhat_plus;
}
// #if 0
__attribute__((noreturn)) int main(void) {
        init_ukf_system();
        float soc_estimate = 100.0f;
        float state_covariance = 1.0f;
        float current = 12.49f;
        float temperature = 25.0f;
        float dt = 1.0f;
        int max_iterations = 1000;
        for (int k = 0; k < max_iterations && soc_estimate > 1.1f; k++) {
            ukf_step(current, voltage_measurement[k], temperature, dt, &soc_estimate, &state_covariance);
            #ifdef DEBUGG
            printf("%.4f %.4f %.4f\n",voltage_measurement[k],soc_estimate,current);
            #endif
        }
        #ifdef ARMCM55
        while( 1 )
           {
           	__asm volatile("nop");
           }
        #endif

}

// #endif

// int main(void){
//     printf("%.4f %.4f %.4f %.4f\n",ocv_table[1][6],ocv_table[1][7],ocv_table[0][6],ocv_table[0][7]);
// }