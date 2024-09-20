#include <stdio.h>
#include <stdlib.h>
#include <papi.h>

#define N 1024  // Size of matrices (must be a power of 2 for simplicity)

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

// Helper functions for Strassen's algorithm
void add_matrices(double* A, double* B, double* C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = A[i * n + j] + B[i * n + j];
        }
    }
}

void subtract_matrices(double* A, double* B, double* C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = A[i * n + j] - B[i * n + j];
        }
    }
}

void strassen(double* A, double* B, double* C, int n) {
    if (n <= 2) {
        // Base case: standard matrix multiplication
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                C[i * n + j] = 0;
                for (int k = 0; k < n; k++) {
                    C[i * n + j] += A[i * n + k] * B[k * n + j];
                }
            }
        }
        return;
    }

    int newSize = n / 2;
    double* A11 = (double*)malloc(newSize * newSize * sizeof(double));
    double* A12 = (double*)malloc(newSize * newSize * sizeof(double));
    double* A21 = (double*)malloc(newSize * newSize * sizeof(double));
    double* A22 = (double*)malloc(newSize * newSize * sizeof(double));
    double* B11 = (double*)malloc(newSize * newSize * sizeof(double));
    double* B12 = (double*)malloc(newSize * newSize * sizeof(double));
    double* B21 = (double*)malloc(newSize * newSize * sizeof(double));
    double* B22 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M1 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M2 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M3 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M4 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M5 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M6 = (double*)malloc(newSize * newSize * sizeof(double));
    double* M7 = (double*)malloc(newSize * newSize * sizeof(double));
    double* temp1 = (double*)malloc(newSize * newSize * sizeof(double));
    double* temp2 = (double*)malloc(newSize * newSize * sizeof(double));

    if (!A11 || !A12 || !A21 || !A22 || !B11 || !B12 || !B21 || !B22 || !M1 || !M2 || !M3 || !M4 || !M5 || !M6 || !M7 || !temp1 || !temp2) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    // Divide matrices into submatrices
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            A11[i * newSize + j] = A[i * n + j];
            A12[i * newSize + j] = A[i * n + (j + newSize)];
            A21[i * newSize + j] = A[(i + newSize) * n + j];
            A22[i * newSize + j] = A[(i + newSize) * n + (j + newSize)];
            B11[i * newSize + j] = B[i * n + j];
            B12[i * newSize + j] = B[i * n + (j + newSize)];
            B21[i * newSize + j] = B[(i + newSize) * n + j];
            B22[i * newSize + j] = B[(i + newSize) * n + (j + newSize)];
        }
    }

    // M1 = (A11 + A22) * (B11 + B22)
    add_matrices(A11, A22, temp1, newSize);
    add_matrices(B11, B22, temp2, newSize);
    strassen(temp1, temp2, M1, newSize);

    // M2 = (A21 + A22) * B11
    add_matrices(A21, A22, temp1, newSize);
    strassen(temp1, B11, M2, newSize);

    // M3 = A11 * (B12 - B22)
    subtract_matrices(B12, B22, temp2, newSize);
    strassen(A11, temp2, M3, newSize);

    // M4 = A22 * (B21 - B11)
    subtract_matrices(B21, B11, temp2, newSize);
    strassen(A22, temp2, M4, newSize);

    // M5 = (A11 + A12) * B22
    add_matrices(A11, A12, temp1, newSize);
    strassen(temp1, B22, M5, newSize);

    // M6 = (A21 - A11) * (B11 + B12)
    subtract_matrices(A21, A11, temp1, newSize);
    add_matrices(B11, B12, temp2, newSize);
    strassen(temp1, temp2, M6, newSize);

    // M7 = (A12 - A22) * (B21 + B22)
    subtract_matrices(A12, A22, temp1, newSize);
    add_matrices(B21, B22, temp2, newSize);
    strassen(temp1, temp2, M7, newSize);

    // C11 = M1 + M4 - M5 + M7
    add_matrices(M1, M4, temp1, newSize);
    subtract_matrices(temp1, M5, temp2, newSize);
    add_matrices(temp2, M7, C, newSize);

    // C12 = M3 + M5
    add_matrices(M3, M5, C + newSize, newSize);

    // C21 = M2 + M4
    add_matrices(M2, M4, C + newSize * n, newSize);

    // C22 = M1 - M2 + M3 + M6
    subtract_matrices(M1, M2, temp1, newSize);
    add_matrices(temp1, M3, temp2, newSize);
    add_matrices(temp2, M6, C + newSize * (n + 1), newSize);

    // Free allocated memory
    free(A11);
    free(A12);
    free(A21);
    free(A22);
    free(B11);
    free(B12);
    free(B21);
    free(B22);
    free(M1);
    free(M2);
    free(M3);
    free(M4);
    free(M5);
    free(M6);
    free(M7);
    free(temp1);
    free(temp2);
}

int main() {
    int EventSet = PAPI_NULL;
    long long values[5]; // Array to store event counts
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
    if (PAPI_add_event(EventSet, PAPI_DP_OPS) != PAPI_OK) {
        printf("PAPI_DP_OPS failed");
        return 1;
    }

    if (PAPI_add_event(EventSet, PAPI_L1_TCM) != PAPI_OK) {
        printf("PAPI_L1_TCM failed");
        return 1;
    }

    if (PAPI_add_event(EventSet, PAPI_L2_TCM) != PAPI_OK) {
        printf("PAPI_L2_TCM failed");
        return 1;
    }

    if (PAPI_add_event(EventSet, PAPI_L3_TCM) != PAPI_OK) {
        printf("PAPI_L3_TCM failed");
        return 1;
    }

    // Start counting events
    retval = PAPI_start(EventSet);
    if (retval != PAPI_OK) {
        printf("PAPI_start failed");
        return 1;
    }

    // Initialize matrices
    initialize_matrices(A, B, C, N);

    // Measure time and perform matrix multiplication
    long long start_time = PAPI_get_real_usec();
    strassen(A, B, C, N);
    long long end_time = PAPI_get_real_usec();

    // Stop counting events
    retval = PAPI_stop(EventSet, values);
    if (retval != PAPI_OK) {
        printf("Error in PAPI_stop");
        return 1;
    }

    // Calculate CPI
    double cycles = (double)values[0]; // Total Cycles
    double instructions = (double)values[1]; // Total Instructions
    double cpi = (instructions != 0) ? (cycles / instructions) : 0;

    // Print event counts
    printf("Wallclock_Time: %lld ms\n", ((end_time - start_time) / 1000));
    printf("PAPI_L1_TCM: %lld\n", values[1]);
    printf("PAPI_L2_TCM: %lld\n", values[2]);
    printf("PAPI_L3_TCM: %lld\n", values[3]);

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
