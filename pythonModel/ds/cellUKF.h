#ifndef CELLUKF_H
#define CELLUKF_H

/**
 * Cell UKF (Unscented Kalman Filter) for Battery SoC Estimation
 * 
 * This library implements an Unscented Kalman Filter for estimating
 * Battery State of Charge (SoC) using terminal voltage measurements.
 */

/**
 * Initialize the UKF system with OCV data and ESR model
 * @param ocv_file: Path to CSV file containing OCV-SoC lookup data
 * @return: 0 on success, -1 on failure
 */
int init_ukf_system(const char* ocv_file);

/**
 * Perform one UKF estimation step
 * @param current: Battery current in Amperes (positive for discharge)
 * @param voltage_measurement: Measured terminal voltage in Volts
 * @param temperature: Temperature in Celsius
 * @param dt: Time step in seconds
 * @param soc_estimate: Current SoC estimate (%) - updated by this function
 * @param state_covariance: State covariance - updated by this function
 * @return: Updated SoC estimate in %
 */
float ukf_step(float current, float voltage_measurement, float temperature,
                float dt, float* soc_estimate, float* state_covariance);

/**
 * Clean up UKF system resources
 */
void cleanup_ukf_system(void);

/**
 * Initialize UKF weights (internal function)
 */
void init_ukf_weights(void);

/**
 * 3x3 Cholesky decomposition (internal function)
 * @param A: Input symmetric positive definite matrix
 * @param L: Output lower triangular matrix
 */
void cholesky_3x3(float A[3][3], float L[3][3]);

#endif /* CELLUKF_H */
