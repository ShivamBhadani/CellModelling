#ifndef OCV_SOC_H
#define OCV_SOC_H

/**
 * OCV (Open Circuit Voltage) Lookup with Bilinear Interpolation
 * 
 * This library provides functions to interpolate OCV values from a CSV lookup table
 * based on State of Charge (SoC) and Temperature using bilinear interpolation.
 */

/**
 * Load OCV data from CSV file into the lookup table
 * @param filename: Path to the CSV file containing OCV data
 * @return: 0 on success, -1 on failure
 */
int load_ocv_data(const char* filename);

/**
 * Get OCV value using bilinear interpolation
 * @param soc: State of Charge (0.0 to 100.0 %)
 * @param temp: Temperature in Celsius (-40.0 to 90.0 °C)
 * @return: OCV value in Volts, or -1.0 on error
 * 
 * Note: Values outside the supported ranges will be clamped to the boundaries.
 *       Call load_ocv_data() before using this function.
 */
double get_ocv_bilinear(double soc, double temp);

/**
 * Print the OCV table for debugging purposes
 */
void print_ocv_table(void);

#endif /* OCV_SOC_H */