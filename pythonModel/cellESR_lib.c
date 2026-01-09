#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Maximum number of RZ elements for embedded systems
#define MAX_RZ_ELEMENTS 20

// Structure for cellESR class
struct cellESR {
    // Basic parameters
    double R0;
    double R0multiplier;
    double esrDC;
    double deltaV;
    int time;
    
    // Arrays for RZ elements
    double R_values[MAX_RZ_ELEMENTS];
    double tau_values[MAX_RZ_ELEMENTS];
    int indices[MAX_RZ_ELEMENTS];
    double ix_values[MAX_RZ_ELEMENTS];
    int num_elements;
    
};

// Export functions for shared library
#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

// Function prototypes with EXPORT
EXPORT cellESR_t* cellESR_create(double R0, double* RZlist_flat, int num_RZ);
EXPORT void cellESR_destroy(cellESR_t* cell);
EXPORT double cellESR_calculateDeltaV(cellESR_t* cell, double i_new, double dt);
EXPORT double cellESR_get_R0(cellESR_t* cell);
EXPORT double cellESR_get_esrDC(cellESR_t* cell);
EXPORT int cellESR_get_time(cellESR_t* cell);
EXPORT double cellESR_get_deltaV(cellESR_t* cell);
EXPORT int cellESR_get_num_elements(cellESR_t* cell);
EXPORT void cellESR_get_R_values(cellESR_t* cell, double* output);
EXPORT void cellESR_get_tau_values(cellESR_t* cell, double* output);
EXPORT void cellESR_get_indices(cellESR_t* cell, int* output);
EXPORT void cellESR_get_ix_values(cellESR_t* cell, double* output);

// Create and initialize cellESR structure
EXPORT cellESR_t* cellESR_create(double R0, double* RZlist_flat, int num_RZ) {
    cellESR_t* cell = (cellESR_t*)malloc(sizeof(cellESR_t));
    if (!cell) return NULL;
    
    // Initialize basic parameters
    cell->R0 = R0;
    cell->R0multiplier = 1.0;
    cell->esrDC = R0;
    cell->deltaV = 0.0;
    cell->time = 0;
    cell->num_elements = (num_RZ < MAX_RZ_ELEMENTS) ? num_RZ : MAX_RZ_ELEMENTS;
    
    // Process RZlist (flattened array: [R0, Z0, type0, R1, Z1, type1, ...])
    for (int i = 0; i < cell->num_elements; i++) {
        double R = RZlist_flat[i * 3];
        double Z = RZlist_flat[i * 3 + 1];  
        int type = (int)RZlist_flat[i * 3 + 2];  // 0 or 1
        
        cell->R_values[i] = R;
        cell->ix_values[i] = 0.0;
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
        free(cell);
    }
}

// Calculate delta V using ESR model
EXPORT double cellESR_calculateDeltaV(cellESR_t* cell, double i_new, double dt) {
    if (!cell) return 0.0;
    
    // Update time
    cell->time += dt;
    
    // Calculate ESR contribution
    double ESRout = cell->R0 * cell->R0multiplier;
    double dV = i_new * ESRout;
    
    // Process each RZ element
    for (int i = 0; i < cell->num_elements; i++) {
        // Update ix_values using: ix_next = (ix_previous + i_new * dt / tau) / (1 + dt / tau)
        cell->ix_values[i] = (cell->ix_values[i] + i_new * dt / cell->tau_values[i]) / 
                             (1.0 + dt / cell->tau_values[i]);
        
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
    }
    
    cell->deltaV = dV;
    return dV;
}

// Getter functions for accessing private data
EXPORT double cellESR_get_R0(cellESR_t* cell) {
    return cell ? cell->R0 : 0.0;
}

EXPORT double cellESR_get_esrDC(cellESR_t* cell) {
    return cell ? cell->esrDC : 0.0;
}

EXPORT int cellESR_get_time(cellESR_t* cell) {
    return cell ? cell->time : 0;
}

EXPORT double cellESR_get_deltaV(cellESR_t* cell) {
    return cell ? cell->deltaV : 0.0;
}

// Additional getter functions for arrays
EXPORT int cellESR_get_num_elements(cellESR_t* cell) {
    return cell ? cell->num_elements : 0;
}

EXPORT void cellESR_get_R_values(cellESR_t* cell, double* output) {
    if (cell && output) {
        for (int i = 0; i < cell->num_elements; i++) {
            output[i] = cell->R_values[i];
        }
    }
}

EXPORT void cellESR_get_tau_values(cellESR_t* cell, double* output) {
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

EXPORT void cellESR_get_ix_values(cellESR_t* cell, double* output) {
    if (cell && output) {
        for (int i = 0; i < cell->num_elements; i++) {
            output[i] = cell->ix_values[i];
        }
    }
}