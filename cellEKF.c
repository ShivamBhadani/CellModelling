/*
 * cellEKF.c
 *
 *  Created on: 13-Feb-2026
 *      Author: shivam
 */


#include <math.h>
#include "y_true.h"
#include "ocv_table.h"
//#include <stdlib.h>
#include "cellESR.h"

#ifdef DEBUGG
#include <stdio.h>
#endif

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

// Create and initialize cellESR structure
cellESR_t cell_t = {0};
cellESR_t* cell=&cell_t;
cellESR_t* cellESR_create(float R0, float* RZlist_flat, int num_RZ) {
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
void cellESR_destroy(cellESR_t* cell) {
    if (cell) {
//        free(cell);
    }
}

// Calculate delta V using ESR model
float cellESR_calculateDeltaV(cellESR_t* cell, float i_new, float dt) {
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
        if (fabs(cell->ix_values[i]) > fabs(i_new)) {
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
        //printf(" dV: %.4f ix: %.4f ",dV, cell->ix_values[i]);
    }

    cell->deltaV = dV;
    return dV;
}


// Temperature array: -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
static const float temperatures[TEMP_POINTS] = {
    -40.0f, -30.0f, -20.0f, -10.0f, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f
};

// Global objects
static cellESR_t* esr_model ;
static cellESR_t* esr_model_list[26];

//static int ocv_initialized = 0;

int init_ukf_system(void) {

    // Initialize ESR model
    float RZlist_flat[] = {1e-3f, 10.0f, 0.0f, 1.5e-3f, 1e4f, 0.0f};  // [R1,tau1,idx1, R2,tau2,idx2]

    esr_model = cellESR_create(0.015f, RZlist_flat, 2);
    if (!esr_model) {
//        printf("Error: Failed to create ESR model\n");
        return -1;
    }

    return 0;
}

// Global variable to store the last computed OCV derivative
static float last_ocv_derivative = 0.0f;


float get_ocv_bilinear(float soc, float temp) {
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

    // Perform bilinear interpolation
    float temp_weight_inv = 1.0f - temp_weight;
    float soc_weight_inv = 1.0f - soc_weight;

    float r1 = q11 * temp_weight_inv + q12 * temp_weight;
    float r2 = q21 * temp_weight_inv + q22 * temp_weight;
    //printf("sl- %d su- %d tl- %d tu- %d q11-%.4f q12-%.4f q21-%.4f q22-%.4f r1-%.4f r2-%.4f ",soc_lower,soc_upper,temp_lower,temp_upper,q11,q12,q21,q22,r1,r2);

    // Compute derivative analytically: d(OCV)/d(soc) = r2 - r1
    // This is the derivative per 1% SOC change
    last_ocv_derivative = r2 - r1;

    // int y=get_cycle_count();
    // bilinearCounts+=y-x;
    return r1 * soc_weight_inv + r2 * soc_weight;
}

// Returns the OCV derivative computed during the last call to get_ocv_bilinear
static float ocv_derivative_soc(float soc, float temp) {
    // Simply return the derivative that was computed in get_ocv_bilinear
    // Note: get_ocv_bilinear must be called first with the same soc and temp
    return last_ocv_derivative;
}

static const float ocvSoH_gain_factor = 1.0f;

// Process and measurement noise covariances
static const float SigmaW = 1.0f / 6.0f;  // Process noise covariance (matches Python)
static const float SigmaV = 0.01f;         // Measurement noise covariance (matches Python)

// Battery parameters
static float MaxCapacity = 65000.0f;  // mAh
static float coulombEfficiency = 0.97f;  // Coulombic efficiency
static float SoH = 100.0f;  // State of Health (%)

// Global objects
static cellESR_t* esr_model ;
static int ocv_initialized = 1;

int init_filter_system(void) {
    return init_ukf_system();
}

float ekf_step(float current, float voltage_measurement, float temperature,
               float dt, float* soc_estimate, float* state_covariance) {
    if (!ocv_initialized || !esr_model) {
        return *soc_estimate;
    }

    // Current estimates
    float xhat = *soc_estimate;
    float SigmaX = *state_covariance;
    if (!isfinite(xhat)) {
        xhat = 0.0f;
    }
    if (!isfinite(SigmaX) || SigmaX < 1e-6f) {
        SigmaX = 1e-6f;
    }

    // === TIME UPDATE ===
    float discharge = current * dt; // As
    float dSoC = 100.0f * discharge / (3.6f * coulombEfficiency * MaxCapacity * SoH / 100.0f);
    float xhat_minus = xhat - dSoC;
    float SigmaX_minus = SigmaX + SigmaW;

    // === MEASUREMENT UPDATE ===
    float soc_bounded = xhat_minus;
    if (soc_bounded < 0.0f) soc_bounded = 0.0f;
    if (soc_bounded > 100.0f) soc_bounded = 100.0f;

    float deltaV = cellESR_calculateDeltaV(esr_model, current, dt);
    float ocv = get_ocv_bilinear(soc_bounded, temperature);
    float yhat = ocv * ocvSoH_gain_factor + deltaV;
    //printf("soc_bounded: %.4f dSoC: %.4f xhat_minus: %.4f SigmaX_minus: %.4f deltaV: %.4f ocv: %.4f yhat: %.4f ",soc_bounded,dSoC,xhat_minus,SigmaX_minus,deltaV,ocv,yhat);

    float H = ocv_derivative_soc(soc_bounded, temperature) * ocvSoH_gain_factor;
    //printf("H: %f ",H);
    float S = H * SigmaX_minus * H + SigmaV;
    float K = (S > 1e-12f) ? (SigmaX_minus * H / S) : 0.0f;

    float innovation = voltage_measurement - yhat;
    float xhat_plus = xhat_minus + K * innovation;
    float SigmaX_plus = (1.0f - K * H) * SigmaX_minus;

    if (SigmaX_plus < 1e-6f) {
        SigmaX_plus = 1e-6f;
    }

    if (xhat_plus < 0.0f) xhat_plus = 0.0f;
    if (xhat_plus > 100.0f) xhat_plus = 100.0f;

    *soc_estimate = xhat_plus;
    *state_covariance = SigmaX_plus;
    return xhat_plus;
}

float filter_step(float current, float voltage_measurement, float temperature,
                  float dt, float* soc_estimate, float* state_covariance) {

return ekf_step(current, voltage_measurement, temperature, dt, soc_estimate, state_covariance);

}
// #if 0
int main(void) {
    init_filter_system();
    float soc_estimate = 100.0f;
    float state_covariance = 1.0f;
    float current = 12.49f;
    float temperature = 25.0f;
    float dt = 1.0f;
    int max_iterations = 1000;
    for (int k = 0; k < max_iterations && soc_estimate > 1.1f; k++) {
        filter_step(current, voltage_measurement[k], temperature, dt, &soc_estimate, &state_covariance);
        #ifdef DEBUGG
        printf("%.4f\n",soc_estimate);
        #endif
    }
    #ifdef ARMCM55
    while( 1 )
       {
       	__asm volatile("nop");
       }
    #endif

    return 0;
}
