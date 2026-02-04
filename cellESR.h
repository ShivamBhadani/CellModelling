#ifndef CELLESR_H
#define CELLESR_H

#ifdef __cplusplus
extern "C" {
#endif

/* Platform export macro */
#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

/* Forward declaration (opaque type for users) */
typedef struct cellESR cellESR_t;

/* Constructor / destructor */
EXPORT cellESR_t* cellESR_create(float R0, float* RZlist_flat, int num_RZ);
EXPORT void cellESR_destroy(cellESR_t* cell);

/* Core functionality */
EXPORT float cellESR_calculateDeltaV(cellESR_t* cell, float i_new, float dt);

/* Getters */
EXPORT float cellESR_get_R0(cellESR_t* cell);
EXPORT float cellESR_get_esrDC(cellESR_t* cell);
EXPORT int    cellESR_get_time(cellESR_t* cell);
EXPORT float cellESR_get_deltaV(cellESR_t* cell);

/* RZ element access */
EXPORT int  cellESR_get_num_elements(cellESR_t* cell);
EXPORT void cellESR_get_R_values(cellESR_t* cell, float* output);
EXPORT void cellESR_get_tau_values(cellESR_t* cell, float* output);
EXPORT void cellESR_get_indices(cellESR_t* cell, int* output);
EXPORT void cellESR_get_ix_values(cellESR_t* cell, float* output);

#ifdef __cplusplus
}
#endif

#endif /* CELLESR_H */
