/*
 * cellEKF.c
 *
 *  Created on: 13-Feb-2026
 *      Author: shivam
 */


#include <math.h>
#include "y_true.h"
#include "ocv_table.h"

#include <arm_mve.h>

//#define DEBUGG

#ifdef DEBUGG
#include <stdio.h>
#endif

// max number of RZ elements for embedded systems
#define NUM_CELLS 26
#define MAX_RZ_ELEMENTS 20
#define NUM_RZ_ELEMENTS 2


// Structure for cellESR class
struct cellESR {
    // basic params
    float R0[NUM_CELLS];
    float R0multiplier[NUM_CELLS];
    float esrDC[NUM_CELLS];
    float deltaV[NUM_CELLS];
    int time;

    // arrays for RZ elements
    float R_values[MAX_RZ_ELEMENTS][NUM_CELLS];
    float tau_values[MAX_RZ_ELEMENTS][NUM_CELLS];
    int indices[MAX_RZ_ELEMENTS][NUM_CELLS]; // denotes type RC or LR
    float ix_values[MAX_RZ_ELEMENTS][NUM_CELLS];
    int num_elements;

};

typedef struct cellESR cellESR_t;

// Create and initialize cellESR structure
cellESR_t cell_t = {0};
cellESR_t* cell = &cell_t;
void cellESR_create(float R0[], float *RZlist_flat_p, int num_RZ) {
    // cellESR_t* cell= (cellESR_t*)(malloc ? malloc(sizeof(cellESR_t)) : cell);
    // if (!cell) return NULL;
	for(int i=0; i<NUM_CELLS; i+=4){

		mve_pred16_t p = vctp32q(NUM_CELLS-i);
//		cell->R0 = R0;
		float32x4_t r0_vec = vld1q_z_f32(&R0[i], p);   // load 4 elements
		vst1q_p(&cell->R0[i], r0_vec, p);          // store 4 elements safely

//		cell->R0multiplier = 1.0f;
		float32x4_t r0mul_vec = vdupq_n_f32(1.0f);             // all lanes = 1.0
		vst1q_p(&cell->R0multiplier[i], r0mul_vec, p);

//		cell->esrDC = R0;
		float32x4_t esrDC_vec = r0_vec;

//		cell->deltaV = 0.0f;
		float32x4_t deltaV_vec = vdupq_n_f32(0.0f);
		vst1q_p(&cell->deltaV[i], deltaV_vec, p);

		// Process RZlist (flattened array: [R0, Z0, type0, R1, Z1, type1, ...])
		cell->num_elements = (num_RZ < MAX_RZ_ELEMENTS) ? num_RZ : MAX_RZ_ELEMENTS;

		for (int j = 0; j < cell->num_elements; j++) {

			//float R = RZlist_flat[j * 3];
			//cell->R_values[j] = R;
			//Load R from flat SoA array for element j
			float32x4_t R_vec = vld1q_z_f32(&RZlist_flat_p[NUM_CELLS*(j*3+0)+i], p);
			// Store directly into cell->R_values[j] for these 4 cells
			float *R_values_p = (float*)cell->R_values;
			vst1q_p(&R_values_p[NUM_CELLS*j+i], R_vec, p);

//			float Z = RZlist_flat[j * 3 + 1];
			float32x4_t Z_vec = vld1q_z_f32(&RZlist_flat_p[NUM_CELLS*(j*3+1)+i], p);

//			int type = (int)RZlist_flat[j * 3 + 2];  // 0 or 1
			float32x4_t type = vld1q_z_f32(&RZlist_flat_p[NUM_CELLS*(j*3+2)+i], p);
			float *ix_values_p = (float*)cell->ix_values;
			vst1q_p(&ix_values_p[NUM_CELLS*j+i], vdupq_n_f32(0.0f), p);
			int32x4_t ii = vcvtq_s32_f32(type);
			int *indices_p = (int*)cell->indices;
			vst1q_p(&indices_p[NUM_CELLS*j+i], ii, p);

			float32x4_t mul_vec = vmulq_m_f32(vdupq_n_f32(0.0f), Z_vec, R_vec, p);
			float32x4_t div_vec = vdupq_n_f32(0.0f);
//			for(int zz=0; zz<4; zz++){
//				((float*)&div_vec)[zz] = ((float*)&Z_vec)[i]*((float*)&R_vec)[i];
//			}

			// Create mask for type == 1
			mve_pred16_t mask0 = vcmpeqq_m_f32(type, vdupq_n_f32(0.0f), p);
			mve_pred16_t mask1 = vcmpeqq_m_f32(type, vdupq_n_f32(1.0f), p);
//			mask = vandq_m(mask, p);                   // combine with tail


//			float32x4_t tau_vec = vmovq_m	(mask, div_vec, mul_vec);
//			if (type == 0) {
				// tau = R * Z, and add R to esrDC
//				cell->tau_values[j] = R * Z;

			// Store into cell->tau_values[j] for 4 cells
			float* tau_values_p = (float*)cell->tau_values;
			vst1q_p(&tau_values_p[NUM_CELLS*j+i], mul_vec, mask0);
			vst1q_p(&tau_values_p[NUM_CELLS*j+i], div_vec, mask1);

			vaddq_m(vdupq_n_f32(0.0f), esrDC_vec, R_vec, mask0);

//			float32_t s = 0.0f;
//			int lanes = (26-i)/4?4:(26-i)%4;
//			for(int k=0;i<lanes;k++){
//				s+=vgetq_lane_f32(R_vec, i);
//			}


			// Select tau per lane
//				cell->esrDC += R;
			// For esrDC, only lanes with type==0 add R
//			float32x4_t esrDC_contrib = vbslq_f32(mask, vdupq_n_f32(0.0f), R_vec);
//			float esrDC_sum = vaddvq_m(R_vec, mask0);   // horizontal sum across vector
//			cell->esrDC += esrDC_sum;
//			} else {
				// tau = Z / R
//				cell->tau_values[j] = Z / R;
//			}
		}
		vst1q_p(&cell->esrDC[i], esrDC_vec, p);

	}
	cell->time = 0;
    return;
}

// Calculate delta V using ESR model
float32x4_t cellESR_calculateDeltaV_x4(volatile cellESR_t* cell, float32x4_t i_new, float dt, mve_pred16_t p, int32_t i) {
//    if (!cell) return 0.0f;
	float32x4_t dV = vdupq_n_f32(0.0f);
	cell->time += dt;

	float32x4_t ESRout = vld1q_z_f32(&cell->R0[i], p);
	float32x4_t multiplier = vld1q_z_f32(&cell->R0multiplier[i], p);
	ESRout = ESRout * multiplier;
	dV = i_new * ESRout;

		// Process each RZ element
	for (int j = 0; j < cell->num_elements; j++) {
			// Update ix_values using: ix_next = (ix_previous + i_new * dt / tau) / (1 + dt / tau)
//			cell->ix_values[i] = (cell->ix_values[i] + i_new * dt / cell->tau_values[i]) /
//								 (1.0f + dt / cell->tau_values[i]);
		float32_t *R_values_p = (float32_t*)cell->R_values;
		float *ix_values_p = (float*)cell->ix_values;
		float *tau_values_p = (float*)cell->tau_values;
		int *indices_p = (int*)cell->indices;
		float32x4_t R_values_vec = vld1q_z_f32(&R_values_p[NUM_CELLS*j+i], p);
		float32x4_t ix_values_vec = vld1q_z_f32(&ix_values_p[NUM_CELLS*j+i], p);
		float32x4_t tau_values_vec = vld1q_z_f32(&tau_values_p[NUM_CELLS*j+i], p);
		float32x4_t inv_tau_values_vec = 1.0f/tau_values_vec;
		float32x4_t new_ix_values_vec_num = vfmaq_m_f32(ix_values_vec, i_new, inv_tau_values_vec, p);
		float32x4_t new_ix_values_vec_den = vfmaq_m_n_f32(vdupq_n_f32(1.0f), inv_tau_values_vec, dt, p);
		int32x4_t cell_ind_vec = vld1q_z_s32(&indices_p[NUM_CELLS*j+i], p);
		ix_values_vec = new_ix_values_vec_num/new_ix_values_vec_den;

			// Clamp values to prevent overshoot // unimplemented
//			if (fabs(cell->ix_values[i]) > fabs(i_new)) {
//			mve_pred16_t clamp = vcmpgeq_m_f32(ix_values_vec, i_new, p);
//				cell->ix_values[i] = i_new;
//			}
//
		mve_pred16_t mask1 = vcmpeqq_m_s32(cell_ind_vec, vdupq_n_s32(1), p);
//			mve_pred16_t mask0 = vcmpeqq_m_s32(vdupq_n_s32(0), cell_ind_vec, vdupq_n_s32(0), p);

			// Calculate deltaV contribution based on type
//			if (cell->indices[i] == 1) {
				// For indices == 1: dV += R * (i_new - ix_values)
//				dV += cell->R_values[i] * (i_new - cell->ix_values[i]);
		ix_values_vec = vsubq_m_f32(ix_values_vec, ix_values_vec, i_new, mask1);
		vst1q_p_f32(&ix_values_p[NUM_CELLS*j+i], ix_values_vec, p);
		dV = dV + R_values_vec * ix_values_vec;

//			} else {
				// For indices == 0: dV += R * ix_values
//				dV += cell->R_values[i] * cell->ix_values[i];
//			}
	}

	vst1q_p_f32(&cell->deltaV[i], dV, p);
	dV;
    return dV;
}


// Temperature array: -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
static const float temperatures[TEMP_POINTS] = {
    -40.0f, -30.0f, -20.0f, -10.0f, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f
};

// Global objects
//static cellESR_t* esr_model ;
//static cellESR_t* esr_model_list[26];

//static int ocv_initialized = 0;

void init_ekf_system(void) {

    // Initialize ESR model
	const int numRZ = 2;
    float RZlist_flat_cell[6] = {1e-3f, 10.0f, 0.0f, 1.5e-3f, 1e4f, 0.0f};  // [R1,tau1,idx1, R2,tau2,idx2]
    float RZlist_flat[6][NUM_CELLS];
//    float two_d_arr[2][2] = {{1.0f,2.0f},{3.0f,4.0f}};
    float *RZlist_flat_p = (float*)(RZlist_flat);
    float R0_cell = 0.015f;
    float R0[NUM_CELLS] = {0};
    for(int i=0; i < numRZ; i++){
    	for(int j=0; j < NUM_CELLS; j+=4){
//    		RZlist_flat_p[(3*i+0)*26+j]=RZlist_flat_cell[3*i+0];
//    		RZlist_flat_p[(3*i+1)*26+j]=RZlist_flat_cell[3*i+1];
//    		RZlist_flat_p[(3*i+2)*26+j]=RZlist_flat_cell[3*i+2];
    		mve_pred16_t p = vctp32q(NUM_CELLS-j);
//    		RZlist_flat[i*3][j] =RZlist_flat_cell[i*3]; // for R1
    		vst1q_p(&RZlist_flat_p[(3*i+0)*NUM_CELLS+j], vdupq_n_f32(RZlist_flat_cell[i*3+0]), p);
//    		RZlist_flat[i*3+1][j] =RZlist_flat_cell[i*3+1]; // for Z1
    		vst1q_p(&RZlist_flat_p[(3*i+1)*NUM_CELLS+j], vdupq_n_f32(RZlist_flat_cell[i*3+1]), p);
//    		RZlist_flat[i*3+2][j] =RZlist_flat_cell[i*3+2]; // for type
    		vst1q_p(&RZlist_flat_p[(3*i+2)*NUM_CELLS+j], vdupq_n_f32(RZlist_flat_cell[i*3+2]), p);

    	}
    }
    for(int j=0; j<NUM_CELLS; j+=4){
    	mve_pred16_t p = vctp32q(NUM_CELLS-j);
    	vst1q_p(&R0[j], vdupq_n_f32(R0_cell), p);

    }
    cellESR_create(R0, RZlist_flat_p, numRZ);
//    int x=0;
//    x=4;
    return;
}

// Global variable to store the last computed OCV derivative
//static float last_ocv_derivative[26] = {0};


float32x4_t get_ocv_bilinear(float32x4_t soc,volatile  float32x4_t temp, float ocv_d_SoC[]) {
    // Clamp inputs for corner cases
//    if (soc < 0.0f) soc = 0.0f;
//    if (soc > 100.0f) soc = 100.0f;
//    if (temp < -40.0f) temp = -40.0f;
//    if (temp > 90.0f) temp = 90.0f;

    // Use int soc as index
    int32x4_t soc_lower = vcvtq_s32_f32(soc);
//    int soc_a[] = {99,100,98,100};
//    soc_lower = vld1q_s32(soc_a);
    mve_pred16_t clamp = vcmpneq_s32(soc_lower, vdupq_n_s32(100));
    int32x4_t soc_upper = soc_lower;
    soc_upper = vaddq_m(soc_upper, soc_lower, vdupq_n_s32(1), clamp);
    float32x4_t soc_weight = vdupq_n_f32(0.0f);

    // Handle SoC corner cases
//    if (soc_lower < 100) {
//        soc_upper = soc_lower + 1;
    soc_weight = soc - vcvtq_f32_s32(soc_lower);
//    }

    int32x4_t temp_lower = vcvtq_s32_f32(temp / vdupq_n_f32(10.0f))+vdupq_n_s32(4);
    int32x4_t temp_upper = temp_lower + vdupq_n_s32(1);
//    int32x4_t temp_lower = temp_index;
//    int32x4_t temp_upper = temp_index;
    float32x4_t temp_weight = vdupq_n_f32(0.0f);

    // Handle temperature corner cases
//    if (temp_index < 0) {
        // temp < -40, already clamped above
//        temp_lower = temp_upper = 0;
//    } else if (temp_index >= TEMP_POINTS) {
        // temp > 90, already clamped above
//        temp_lower = temp_upper = TEMP_POINTS - 1;
//    } else if (temp_index < TEMP_POINTS - 1) {
        // Normal case: interpolate between two temperature points
//        temp_upper = temp_index + 1;
    float32x4_t temp_at_lower = vldrwq_gather_shifted_offset_f32(temperatures, temp_lower);
//        float32x4_t temp_at_lower = (float32x4_t)(((int32x4_t)(temp/10.0f))*10.0f);
//        float32x4_t temp_at_upper = vldrwq_gather_offset_f32(temperatures, temp_upper);
    temp_weight = (temp - temp_at_lower) *0.1f; // 1/(temp_at_upper - temp_at_lower)
//    }
    // If temp_index == TEMP_POINTS-1, temp_upper stays same as temp_lower

    // Get the four corner points for bilinear interpolation
    float32x4_t q11 = vldrwq_gather_shifted_offset_f32(ocv_table_p, (100-soc_lower)*14+temp_lower);
    float32x4_t q12 = vldrwq_gather_shifted_offset_f32(ocv_table_p, (100-soc_lower)*14+temp_upper);
    float32x4_t q21 = vldrwq_gather_shifted_offset_f32(ocv_table_p, (100-soc_upper)*14+temp_lower);
    float32x4_t q22 = vldrwq_gather_shifted_offset_f32(ocv_table_p, (100-soc_upper)*14+temp_upper);
    // printf("socl--%d socu--%d templ--%d tempu--%d\n",soc_lower,soc_upper,temp_lower,temp_upper);
    // printf("q11-%.4f q12-%.4f q21-%.4f q22-%.4f ",q11,q12,q21,q22);

    // Perform bilinear interpolation
    float32x4_t temp_weight_inv = 1.0f - temp_weight;
    float32x4_t soc_weight_inv = 1.0f - soc_weight;

    float32x4_t r1 = q11 * temp_weight_inv + q12 * temp_weight;
    float32x4_t r2 = q21 * temp_weight_inv + q22 * temp_weight;

    // Compute derivative analytically: d(OCV)/d(soc) = r2 - r1
    // This is the derivative per 1% SOC change
    vst1q_f32(ocv_d_SoC,(r2 - r1));
//    printf("H: %f",ocv_d_SoC[0]);

    // int y=get_cycle_count();
    // bilinearCounts+=y-x;
    return r1 * soc_weight_inv + r2 * soc_weight;
}

// Returns the OCV derivative computed during the last call to get_ocv_bilinear
//static float32x4_t ocv_derivative_soc(float32x4_t soc, float32x4_t temp) {
//    // Simply return the derivative that was computed in get_ocv_bilinear
//    // Note: get_ocv_bilinear must be called first with the same soc and temp
//    return last_ocv_derivative;
//}

static const float ocvSoH_gain_factor = 1.0f;

// Process and measurement noise covariances
static const float SigmaW = 1.0f / 6.0f;  // Process noise covariance (matches Python)
static const float SigmaV = 0.01f;         // Measurement noise covariance (matches Python)

// Battery parameters
static float MaxCapacity = 65000.0f;  // mAh
static float coulombEfficiency = 0.97f;  // Coulombic efficiency
static float SoH = 100.0f;  // State of Health (%)

// Global objects
//static cellESR_t* esr_model ;
//static int ocv_initialized = 1;


float32x4_t ekf_step_x4(float32x4_t current, float32x4_t voltage_measurement, float32x4_t temperature,
		 float32_t dt, float32_t soc_estimate[], float32_t state_covariance[], mve_pred16_t p, int32_t i);


void ekf_step_vec(float32_t current[], float32_t voltage_measurement[], float32_t temperature[],
		float32_t dt, float32_t soc_estimate[], float32_t state_covariance[]){
	for(int32_t i=0; i<NUM_CELLS; i+=4){
		mve_pred16_t p = vctp32q(NUM_CELLS-i);
		float32x4_t current_vec = vld1q_z_f32(&current[i],p);
		float32x4_t voltage_measurement_vec = vld1q_z_f32(&voltage_measurement[i], p);
		float32x4_t temperature_vec = vld1q_z_f32(&temperature[i],p);
//		float32x4_t soc_estimate = vld1q_z(&soc_estimate[i],p);
//		float32x4_t state_covariance = vld1q_z(&state_covariance[i],p);
		(void)ekf_step_x4(current_vec, voltage_measurement_vec, temperature_vec, dt, soc_estimate+i, state_covariance+i, p, i);
	}

}
float32x4_t ekf_step_x4(float32x4_t current, float32x4_t voltage_measurement, float32x4_t temperature,
		 float32_t dt, float32_t soc_estimate[], float32_t state_covariance[], mve_pred16_t p, int32_t i) {
//    if (!ocv_initialized || !esr_model) {
//        return *soc_estimate;
//    }

    // Current estimates

//    float xhat = *soc_estimate;
	float32x4_t xhat = vld1q_z(&soc_estimate[0],p);
//    float SigmaX = *state_covariance;
	float32x4_t SigmaX = vld1q_z(&state_covariance[0],p);
//    if (!isfinite(xhat)) {
//        xhat = 0.0f;
//    }
//    if (!isfinite(SigmaX) || SigmaX < 1e-6f) {
//        SigmaX = 1e-6f;
//    }

    // === TIME UPDATE ===
//    float discharge = current * dt; // As
	float32x4_t discharge = vmulq_m(vdupq_n_f32(0.0f), current, vdupq_n_f32(dt), p);

//    float dSoC = 100.0f * discharge / (3.6f * coulombEfficiency * MaxCapacity * SoH / 100.0f);
	float32x4_t dSoC_n = vmulq_m(vdupq_n_f32(0.0f), discharge, vdupq_n_f32(100.0f*100.0f), p);
	float dSoC_d =  (3.6f * coulombEfficiency * MaxCapacity * SoH);
	float dSoC_d_inv = 1.0f/dSoC_d;
	float32x4_t dSoC = vmulq_m(vdupq_n_f32(0.0f), dSoC_n, vdupq_n_f32(dSoC_d_inv), p);

//    float xhat_minus = xhat - dSoC;
	float32x4_t xhat_minus = vdupq_n_f32(0.0f);
    xhat_minus = vsubq_m(vdupq_n_f32(0.0f), xhat, dSoC, p);

//    float SigmaX_minus = SigmaX + SigmaW;
    float32x4_t SigmaX_minus = vaddq_m(vdupq_n_f32(0.0f), SigmaX, SigmaW, p);

    // === MEASUREMENT UPDATE ===
//    float soc_bounded = xhat_minus;
    float32x4_t soc_bounded = xhat_minus;

//    if (soc_bounded < 0.0f) soc_bounded = 0.0f;
    mve_pred16_t p_low = vcmpleq_m_f32(soc_bounded, vdupq_n_f32(0.0f), p); // true where x < 0

    soc_bounded = vdupq_m_n_f32(soc_bounded, 0.0f, p_low);            // clamp lower

    //    if (soc_bounded > 100.0f) soc_bounded = 100.0f;
    mve_pred16_t p_high = vcmpgeq_m_f32(soc_bounded, vdupq_n_f32(100.0f), p); // true where x < 0
    soc_bounded = vdupq_m_n_f32(soc_bounded, 100.0f, p_high);            // clamp lower
//    p_low = vandq_m(p_low, p);                           // mask out inactive lanes
//    soc_bounded = vbslq_f32(p_low, vdupq_n_f32(100.0f), soc_bounded);            // clamp lower

    float32x4_t deltaV = cellESR_calculateDeltaV_x4(cell, current, dt, p, i);
//    float32x4_t deltaV = vld1q_f32(&cell->)
    float32_t ocv_d_SoC[4] = {0.0f};
    float32x4_t ocv = get_ocv_bilinear(soc_bounded, temperature, ocv_d_SoC);
    float32x4_t yhat = ocv * ocvSoH_gain_factor + deltaV;
    float32x4_t ocv_d_SoC_vec = vld1q_f32(ocv_d_SoC);
    float32x4_t H = ocv_d_SoC_vec * ocvSoH_gain_factor;

    float32x4_t S = H * SigmaX_minus * H + SigmaV;
//    float32x4_t K = (S > 1e-12f) ? (SigmaX_minus * H / S) : 0.0f;
    float32x4_t K = SigmaX_minus * H / S;

    float32x4_t innovation = voltage_measurement - yhat;
    float32x4_t xhat_plus = xhat_minus + K * innovation;
    float32x4_t SigmaX_plus = (1.0f - K * H) * SigmaX_minus;

	mve_pred16_t clamp = vcmpgeq_m_f32(xhat_plus, vdupq_n_f32(100.0f), p);

	xhat_plus = vdupq_m_n_f32(xhat_plus, 100.0f, clamp);


//    if (SigmaX_plus < 1e-6f) {
//        SigmaX_plus = 1e-6f;
//    }

//    if (xhat_plus < 0.0f) xhat_plus = 0.0f;
//    if (xhat_plus > 100.0f) xhat_plus = 100.0f;
//    *state_covariance = SigmaX_plus;
    vst1q_p_f32(state_covariance, SigmaX_plus, p);
//    *soc_estimate = xhat_plus;
    vst1q_p_f32(soc_estimate, xhat_plus, p);

    return xhat_plus;
}

// #if 0
__attribute__((noreturn)) int main(void) {
    init_ekf_system();
    float soc_estimate[NUM_CELLS] = {100.0f};
    float state_covariance[NUM_CELLS] = {1.0f};
    float current[NUM_CELLS] = {12.49f};
    float temperature[NUM_CELLS] = {25.0f};
    float dt = 1.0f;

    int max_iterations = 1000;
    for(int i=1;i<NUM_CELLS;i+=4){
    	mve_pred16_t p = vctp32q(NUM_CELLS-i);
    	vst1q_p_f32((soc_estimate+i),vdupq_n_f32(*soc_estimate),p);
    	vst1q_p_f32((state_covariance+i),vdupq_n_f32(*state_covariance),p);
    	vst1q_p_f32((current+i),vdupq_n_f32(*current),p);
    	vst1q_p_f32((temperature+i),vdupq_n_f32(*temperature),p);
    }
    for (int k = 0; k < max_iterations && *soc_estimate > 1.1f; k++) {
    	float voltage_true[NUM_CELLS] = {0.0f};
    	float voltage_measurement_k = voltage_measurement[k];
    	 for(int i=0;i<NUM_CELLS;i+=4){
    	    mve_pred16_t p = vctp32q(NUM_CELLS-i);
    	    vst1q_p_f32((voltage_true+i),vdupq_n_f32(voltage_measurement_k),p);
    	 }
        ekf_step_vec(current, voltage_true, temperature, dt, soc_estimate, state_covariance);
        #ifdef DEBUGG
//        for(int i=0;i<NUM_CELLS;i++)
        printf("%.4f ",soc_estimate[0]);
        printf("\n");
        #endif
    }
    #ifdef ARMCM55
    while( 1 )
       {
       	__asm volatile("nop");
       }
    #endif
}

