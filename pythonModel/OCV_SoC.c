#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "OCV_SoC.h"

#define SOC_POINTS 101  // SoC from 0 to 100 in 1% increments
#define TEMP_POINTS 14  // Temperature from -40 to 90 in 10°C steps

// Temperature array: -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90
static const double temperatures[TEMP_POINTS] = {
    -40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0
};

// OCV lookup table [SoC][Temperature]
static double ocv_table[SOC_POINTS][TEMP_POINTS];
static int table_initialized = 0;

int load_ocv_data(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open file %s\n", filename);
        return -1;
    }

    char line[2048];
    int row = 0;

    // Skip header line
    if (fgets(line, sizeof(line), file) == NULL) {
        fclose(file);
        return -1;
    }

    // Read data lines
    while (fgets(line, sizeof(line), file) && row < SOC_POINTS) {
        char *token = strtok(line, ",");
        int col = 0;

        // Skip SoC column (first column)
        token = strtok(NULL, ",");
        
        // Read OCV values for each temperature
        while (token && col < TEMP_POINTS) {
            ocv_table[100-row][col] = atof(token);  // SoC decreases in CSV (100 to 0)
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }

    fclose(file);
    table_initialized = 1;
    return 0;
}

void find_interpolation_params(double value, const double* array, int size, 
                              int* lower_idx, int* upper_idx, double* weight) {
    // Handle boundary cases
    if (value <= array[0]) {
        *lower_idx = 0;
        *upper_idx = 0;
        *weight = 0.0;
        return;
    }
    
    if (value >= array[size-1]) {
        *lower_idx = size-1;
        *upper_idx = size-1;
        *weight = 0.0;
        return;
    }

    // Find the bracketing indices
    for (int i = 0; i < size-1; i++) {
        if (value >= array[i] && value <= array[i+1]) {
            *lower_idx = i;
            *upper_idx = i+1;
            *weight = (value - array[i]) / (array[i+1] - array[i]);
            return;
        }
    }
}

double get_ocv_bilinear(double soc, double temp) {
    if (!table_initialized) {
        printf("Error: OCV table not initialized. Call load_ocv_data() first.\n");
        return -1.0;
    }

    // Bound checking
    if (soc < 0.0) soc = 0.0;
    if (soc > 100.0) soc = 100.0;
    if (temp < -40.0) temp = -40.0;
    if (temp > 90.0) temp = 90.0;

    // Create SoC array for interpolation (0 to 100 in 1% steps)
    double soc_array[SOC_POINTS];
    for (int i = 0; i < SOC_POINTS; i++) {
        soc_array[i] = (double)i;
    }

    // Find interpolation parameters for SoC
    int soc_lower, soc_upper;
    double soc_weight;
    find_interpolation_params(soc, soc_array, SOC_POINTS, &soc_lower, &soc_upper, &soc_weight);

    // Find interpolation parameters for Temperature
    int temp_lower, temp_upper;
    double temp_weight;
    find_interpolation_params(temp, temperatures, TEMP_POINTS, &temp_lower, &temp_upper, &temp_weight);

    // Get the four corner points for bilinear interpolation
    double q11 = ocv_table[soc_lower][temp_lower];   // (soc_lower, temp_lower)
    double q12 = ocv_table[soc_lower][temp_upper];   // (soc_lower, temp_upper)
    double q21 = ocv_table[soc_upper][temp_lower];   // (soc_upper, temp_lower)
    double q22 = ocv_table[soc_upper][temp_upper];   // (soc_upper, temp_upper)

    // Perform bilinear interpolation
    double r1 = q11 * (1.0 - temp_weight) + q12 * temp_weight;  // Interpolate along temperature for soc_lower
    double r2 = q21 * (1.0 - temp_weight) + q22 * temp_weight;  // Interpolate along temperature for soc_upper
    double result = r1 * (1.0 - soc_weight) + r2 * soc_weight;  // Interpolate along SoC

    return result;
}

void print_ocv_table() {
    if (!table_initialized) {
        printf("Table not initialized\n");
        return;
    }

    printf("SoC\\Temp");
    for (int j = 0; j < TEMP_POINTS; j++) {
        printf("\t%.0f°C", temperatures[j]);
    }
    printf("\n");

    for (int i = 0; i < SOC_POINTS; i += 10) {  // Print every 10% for brevity
        printf("%d%%", i);
        for (int j = 0; j < TEMP_POINTS; j++) {
            printf("\t%.3f", ocv_table[i][j]);
        }
        printf("\n");
    }
}

#ifdef TEST
int main() {
    // Load the OCV data
    if (load_ocv_data("Sample_OCV_SoC_1.0pct_10deg.csv") != 0) {
        printf("Failed to load OCV data\n");
        return -1;
    }

    // Test cases
    double test_socs[] = {0.0, 25.5, 50.0, 75.3, 100.0};
    double test_temps[] = {-40.0, -15.0, 0.0, 25.0, 60.0, 90.0};

    printf("Test Results:\n");
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 6; j++) {
            double ocv = get_ocv_bilinear(test_socs[i], test_temps[j]);
            printf("SoC: %6.1f%%, Temp: %6.1f°C -> OCV: %7.4f V\n", 
                   test_socs[i], test_temps[j], ocv);
        }
    }

    printf("\nInterpolation Examples:\n");
    printf("SoC: 33.7%%, Temp: 22.5°C -> OCV: %.4f V\n", get_ocv_bilinear(33.7, 22.5));
    printf("SoC: 67.2%%, Temp: -5.5°C -> OCV: %.4f V\n", get_ocv_bilinear(67.2, -5.5));
    printf("SoC: 89.4%%, Temp: 55.3°C -> OCV: %.4f V\n", get_ocv_bilinear(89.4, 55.3));
    printf("20-11: %.4f\n",get_ocv_bilinear(20,11));

    return 0;
}
#endif 