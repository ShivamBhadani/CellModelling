#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cellUKF.h"
#include "cellESR.h"
#include "OCV_SoC.h"

// UKF Parameters
#define Nx 1         // Number of states (SoC)
#define Nxa 3        // Augmented number of states (state + process noise + measurement noise)
#define Ny 1         // Number of measurements (voltage)
#define NUM_SIGMA_POINTS (2*Nxa+1)  // Total sigma points: 2*3+1 = 7

// CDKF parameters
static const double h = 1.732050808;  // sqrt(3) - scaling parameter
static const double inv_2h2 = 1.0 / (2.0 * 3.0);  // 1/(2*h^2) = 1/6

// UKF weights
static double Wmx[NUM_SIGMA_POINTS];     // Mean weights
static double Wcx[NUM_SIGMA_POINTS];     // Covariance weights

// Process and measurement noise covariances
static const double SigmaW = 1.0 / 6.0;  // Process noise covariance (matches Python)
static const double SigmaV = 0.01;         // Measurement noise covariance (matches Python)

// Battery parameters
static double MaxCapacity = 65000;  // mAh
static double coulombEfficiency = 0.97;  // Coulombic efficiency
static double SoH = 100.0;  // State of Health (%)

// Global objects
static cellESR_t* esr_model = NULL;
static int ocv_initialized = 0;

// CSV logging
static FILE* csv_file = NULL;
static double ocvSoH_gain_factor = 1.0;  // Battery OCV SoH gain factor

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

void cholesky_3x3(double A[3][3], double L[3][3]) {
    memset(L, 0, sizeof(double) * 9);
    
    L[0][0] = sqrt(A[0][0]);
    L[1][0] = A[1][0] / L[0][0];
    L[2][0] = A[2][0] / L[0][0];
    
    L[1][1] = sqrt(A[1][1] - L[1][0] * L[1][0]);
    L[2][1] = (A[2][1] - L[2][0] * L[1][0]) / L[1][1];
    
    L[2][2] = sqrt(A[2][2] - L[2][0] * L[2][0] - L[2][1] * L[2][1]);
}

int init_ukf_system(const char* ocv_file) {
    // Initialize weights
    init_ukf_weights();
    
    // Load OCV lookup table
    if (load_ocv_data(ocv_file) != 0) {
        printf("Error: Failed to load OCV data from %s\n", ocv_file);
        return -1;
    }
    ocv_initialized = 1;
    
    // Initialize ESR model
    double RZlist_flat[] = {1e-3, 10.0, 0.0, 1.5e-3, 1e4, 0.0};  // [R1,tau1,idx1, R2,tau2,idx2]
    esr_model = cellESR_create(0.015, RZlist_flat, 2);
    if (!esr_model) {
        printf("Error: Failed to create ESR model\n");
        return -1;
    }
    
    return 0;
}

int init_csv_logging(const char* filename) {
    csv_file = fopen(filename, "w");
    if (!csv_file) {
        printf("Error: Cannot create CSV file %s\n", filename);
        return -1;
    }
    
    // Write CSV header
    fprintf(csv_file, "time,soc_true,soc_estimate,covariance,voltage_measurement,current,temperature\n");
    return 0;
}

void log_to_csv(double time, double soc_true, double soc_estimate, double covariance, 
               double voltage, double current, double temperature) {
    if (csv_file) {
        fprintf(csv_file, "%.1f,%.6f,%.6f,%.8f,%.6f,%.3f,%.1f\n", 
                time, soc_true, soc_estimate, covariance, voltage, current, temperature);
        fflush(csv_file);  // Ensure data is written immediately
    }
}


void cleanup_ukf_system(void) {
    if (esr_model) {
        cellESR_destroy(esr_model);
        esr_model = NULL;
    }
    if (csv_file) {
        fclose(csv_file);
        csv_file = NULL;
    }
}

double ukf_step(double current, double voltage_measurement, double temperature, 
                double dt, double* soc_estimate, double* state_covariance) {
    
    if (!ocv_initialized || !esr_model) {
        printf("Error: UKF system not initialized\n");
        return *soc_estimate;
    }
    
    // Current estimates
    double xhat = *soc_estimate;
    double SigmaX = *state_covariance;
    
    // === STEP 1: GENERATE SIGMA POINTS ===

    double xhata[3] = {xhat, 0.0, 0.0};
    
    // Augmented covariance matrix
    double SigmaXa[3][3] = {
        {SigmaX, 0.0, 0.0},
        {0.0, SigmaW, 0.0},
        {0.0, 0.0, SigmaV}
    };
    
    // Cholesky decomposition
    double L[3][3];
    cholesky_3x3(SigmaXa, L);
    
    // Generate sigma points
    double X[NUM_SIGMA_POINTS][3];  // [sigma_point][state_dimension]
    
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
    
    double discharge[NUM_SIGMA_POINTS];
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        discharge[i] = (current + X[i][1]) * dt; // Discharge in As
    }
    
    double dSoC[NUM_SIGMA_POINTS];
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        dSoC[i] = 100.0 * discharge[i] / (3.6 * coulombEfficiency * MaxCapacity * SoH / 100.0); // Change in SoC considering SoH
    }

    double Xx[NUM_SIGMA_POINTS];  // Predicted state sigma points
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        Xx[i] = X[i][0] - dSoC[i] ;  // state - discharge + process_noise
        // printf("Xx[%d]=%.6f ",i,Xx[i]);
    }
    // printf("\n");

    
    // Predicted state mean
    double xhat_minus = 0.0;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        xhat_minus += Wmx[i] * Xx[i];
    }
    
    // Predicted state covariance
    double SigmaX_minus = 0.0;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        double diff = Xx[i] - xhat_minus;
        SigmaX_minus += Wcx[i] * diff * diff;
    }
    
    // === STEP 3: MEASUREMENT PREDICTION ===
    
    double Y[NUM_SIGMA_POINTS];  // Predicted measurement sigma points
    
    double deltaV = cellESR_calculateDeltaV(esr_model, current, dt);
    
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        // Ensure SoC is within bounds for OCV lookup
        double soc_bounded = Xx[i];
        if (soc_bounded < 0.0) soc_bounded = 0.0;
        if (soc_bounded > 100.0) soc_bounded = 100.0;
        
        double ocv = get_ocv_bilinear(soc_bounded, temperature);
        
        Y[i] = ocv * ocvSoH_gain_factor + deltaV + X[i][2];
    }
    
    // Predicted measurement mean
    double yhat = 0.0;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        yhat += Wmx[i] * Y[i];
    }
    
    // Innovation covariance
    double SigmaY = 0.0;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        double diff = Y[i] - yhat;
        SigmaY += Wcx[i] * diff * diff;
    }
    
    // Cross-covariance
    double SigmaXY = 0.0;
    for (int i = 0; i < NUM_SIGMA_POINTS; i++) {
        SigmaXY += Wcx[i] * (Xx[i] - xhat_minus) * (Y[i] - yhat);
    }
    
    // === STEP 4: MEASUREMENT UPDATE (CORRECTION) ===
    
    // Kalman gain
    double K = SigmaXY / SigmaY;
    
    // Innovation
    double innovation = voltage_measurement - yhat;
    
    // Updated state estimate
    double xhat_plus = xhat_minus + K * innovation;
    
    // Updated state covariance
    double SigmaX_plus = SigmaX_minus - K * SigmaY * K;
    
    // Ensure positive covariance
    if (SigmaX_plus < 1e-6) {
        SigmaX_plus = 1e-6;
    }
    
    // Bound SoC estimate
    if (xhat_plus < 0.0) xhat_plus = 0.0;
    if (xhat_plus > 100.0) xhat_plus = 100.0;
    
    // Update outputs
    *soc_estimate = xhat_plus;
    *state_covariance = SigmaX_plus;
    
    return xhat_plus;
}

#ifdef UKFTEST
int main() {
    
    // Initialize system
    if (init_ukf_system("Sample_OCV_SoC_1.0pct_10deg.csv") != 0) {
        printf("Failed to initialize UKF system\n");
        return -1;
    }
    
    FILE* py_csv = fopen("ytrue_UKF_5hr_10s.csv", "r");
    if (py_csv) {
        printf("Python results file found. Running test.\n");
        double soc_estimate = 100.0;
        double state_covariance = 1.0;
        double current = 12.49;
        double temperature = 25.0;
        double dt = 1.0;
        int max_iterations = 1000;
        
        FILE* csv = fopen("ukf_c_results.csv", "w");
        fprintf(csv, "time,Csoc_estimate\n");
        fprintf(csv, "0.0,%.6f\n", soc_estimate);
        double voltage_measurement[max_iterations];
        char line[256];
            fgets(line, sizeof(line), py_csv); // Skip header
        for(int iter=0; iter<max_iterations; iter++) {
            fgets(line, sizeof(line), py_csv); // Placeholder for actual measurement generation
            voltage_measurement[iter] = strtod(line, NULL);
        }
        for (int k = 0; k < max_iterations && soc_estimate > 1.1; k++) {
            ukf_step(current, voltage_measurement[k], temperature, dt, &soc_estimate, &state_covariance);
            fprintf(csv, "%d,%.6f\n", (int)(k+1), soc_estimate);
        }
        fclose(csv);
        printf("Results saved to ukf_c_results.csv\n");
        cleanup_ukf_system();
        fclose(py_csv);
        
    }
    else {
        printf("file not found\n");
    }
    return 0;

}
#endif 