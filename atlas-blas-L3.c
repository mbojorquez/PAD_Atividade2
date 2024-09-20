#include <stdio.h>
#include <stdlib.h>
#include <papi.h>
#include <cblas.h>

#define N 1024 // Size of matrices

// Function to initialize the matrices
void initialize_matrices(double* A, double* B, double* C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = rand() % 5; // Row-major order
            B[i * n + j] = rand() % 5; // Row-major order
            C[i * n + j] = 0.0;        // Initialize result matrix
        }
    }
}

// Function to print the matrix
void print_matrix(const char* name, double* M, int n) {
    printf("%s:\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%6.2f ", M[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}


int main() {
    int EventSet = PAPI_NULL;
    long long values[6]; // Array to store event counts
    long long s, e;
    int retval;

    // Allocate memory for matrices
    double* A = (double*)malloc(N * N * sizeof(double));
    double* B = (double*)malloc(N * N * sizeof(double));
    double* C = (double*)malloc(N * N * sizeof(double));

    if (A == NULL || B == NULL || C == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    // Initialize PAPI library
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        printf("PAPI_library_init failed: %d\n", retval);
        return 1;
    }

    // Create an event set
    retval = PAPI_create_eventset(&EventSet);
    if (retval != PAPI_OK) {
        printf("PAPI_create_eventset failed: %d\n", retval);
        return 1;
    }

    // Add events to the event set
    if (PAPI_add_event(EventSet, PAPI_DP_OPS) != PAPI_OK ||
        PAPI_add_event(EventSet, PAPI_L1_TCM) != PAPI_OK ||
        PAPI_add_event(EventSet, PAPI_L2_TCM) != PAPI_OK ||
        PAPI_add_event(EventSet, PAPI_L3_TCM) != PAPI_OK) {
        printf("PAPI_add_event failed\n");
        return 1;
    }

    // Start counting events
    retval = PAPI_start(EventSet);
    if (retval != PAPI_OK) {
        printf("PAPI_start failed\n");
        return 1;
    }

    // Initialize matrices

    initialize_matrices(A, B, C, N);
    //print_matrix("Matrix A", A, N);
    //print_matrix("Matrix B", B, N);
    s = PAPI_get_real_usec();

    // Perform matrix multiplication using ATLAS/BLAS
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N,
                 1.0, A, N, B, N, 0.0, C, N);

    e = PAPI_get_real_usec();

    //print_matrix("Result Matrix C", C, N);


    // Stop counting events
    retval = PAPI_stop(EventSet, values);
    if (retval != PAPI_OK) {
        printf("Error in PAPI_stop\n");
        return 1;
    }

    // Calculate CPI
    double cycles = (double)values[0];
    double instructions = (double)values[1];
    double cpi = (instructions != 0) ? (cycles / instructions) : 0;


    // Print event counts
    printf("Wallclock_Time: %lld ms\n", ((e-s)/1000));
    printf("PAPI_L1_TCM: %lld\n", values[1]);
    printf("PAPI_L2_TCM: %lld\n", values[2]);
    printf("PAPI_L3_TCM: %lld\n", values[3]);
    printf("+++\n---\n");  

    // Clean up PAPI
    PAPI_cleanup_eventset(EventSet);
    PAPI_destroy_eventset(&EventSet);
    PAPI_shutdown();

    // Free allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}
