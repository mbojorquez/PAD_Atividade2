#include <stdio.h>
#include <stdlib.h>
#include <papi.h>

#define N 1024  // Size of matrices
#define BLOCK_SIZE 64  // Tamanho do bloco

// Function to initialize the matrices
void initialize_matrices(double** A, double** B, double** C) {
    for (long long int i = 0; i < N; i++) {
        for (long long int j = 0; j < N; j++) {
            A[i][j] = rand() % 5;
            B[i][j] = rand() % 5;
            C[i][j] = 0.0;
        }
    }
}

// Function to print the matrix
void print_matrix(double** M) {
    for (long long int i = 0; i < N; i++) {
        for (long long int j = 0; j < N; j++) {
            printf("%f ", M[i][j]);
        }
        printf("\n");
    }
}


// Multiplicação de matrizes ijk com blocagem 
void matmul_blocked_ijk(double** A, double** B, double** C) {
    int i, j, k, ii, jj, kk;
    for (ii = 0; ii < N; ii += BLOCK_SIZE) {
        for (jj = 0; jj < N; jj += BLOCK_SIZE) {
            for (kk = 0; kk < N; kk += BLOCK_SIZE) {
                // Multiplica blocos de tamanho BLOCK_SIZE x BLOCK_SIZE
                for (i = ii; i < ii + BLOCK_SIZE && i < N; i++) {
                    for (j = jj; j < jj + BLOCK_SIZE && j < N; j++) {
                        for (k = kk; k < kk + BLOCK_SIZE && k < N; k++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

// Multiplicação de matrizes ikj com blocagem 
void matmul_blocked_ikj(double** A, double** B, double** C) {
    int i, j, k, ii, jj, kk;
    for (ii = 0; ii < N; ii += BLOCK_SIZE) {
        for (jj = 0; jj < N; jj += BLOCK_SIZE) {
            for (kk = 0; kk < N; kk += BLOCK_SIZE) {
                // Multiplica blocos de tamanho BLOCK_SIZE x BLOCK_SIZE
                for (i = ii; i < ii + BLOCK_SIZE && i < N; i++) {
                    for (k = kk; k < kk + BLOCK_SIZE && k < N; k++) {
                        for (j= jj; j < jj + BLOCK_SIZE && j < N; j++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

// Multiplicação de matrizes jik com blocagem 
void matmul_blocked_jik(double** A, double** B, double** C) {
    int i, j, k, ii, jj, kk;
    for (ii = 0; ii < N; ii += BLOCK_SIZE) {
        for (jj = 0; jj < N; jj += BLOCK_SIZE) {
            for (kk = 0; kk < N; kk += BLOCK_SIZE) {
                // Process blocks of size BLOCK_SIZE x BLOCK_SIZE
                for (j = jj; j < jj + BLOCK_SIZE && j < N; j++) {
                    for (i = ii; i < ii + BLOCK_SIZE && i < N; i++) {
                        for (k = kk; k < kk + BLOCK_SIZE && k < N; k++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

// Multiplicação de matrizes jki com blocagem 
void matmul_blocked_jki(double** A, double** B, double** C) {
    int i, j, k, ii, jj, kk;
    for (ii = 0; ii < N; ii += BLOCK_SIZE) {
        for (jj = 0; jj < N; jj += BLOCK_SIZE) {
            for (kk = 0; kk < N; kk += BLOCK_SIZE) {
                // Process blocks of size BLOCK_SIZE x BLOCK_SIZE
                for (j = jj; j < jj + BLOCK_SIZE && j < N; j++) {
                    for (k = kk; k < kk + BLOCK_SIZE && k < N; k++) {
                        for (i = ii; i < ii + BLOCK_SIZE && i < N; i++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

// Multiplicação de matrizes kij com blocagem 
void matmul_blocked_kij(double** A, double** B, double** C) {
    int i, j, k, ii, jj, kk;
    for (ii = 0; ii < N; ii += BLOCK_SIZE) {
        for (jj = 0; jj < N; jj += BLOCK_SIZE) {
            for (kk = 0; kk < N; kk += BLOCK_SIZE) {
                // Process blocks of size BLOCK_SIZE x BLOCK_SIZE
                for (k = kk; k < kk + BLOCK_SIZE && k < N; k++) {
                    for (i = ii; i < ii + BLOCK_SIZE && i < N; i++) {
                        for (j = jj; j < jj + BLOCK_SIZE && j < N; j++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

// Multiplicação de matrizes kji com blocagem 
void matmul_blocked_kji(double** A, double** B, double** C) {
    int i, j, k, ii, jj, kk;
    for (ii = 0; ii < N; ii += BLOCK_SIZE) {
        for (jj = 0; jj < N; jj += BLOCK_SIZE) {
            for (kk = 0; kk < N; kk += BLOCK_SIZE) {
                // Process blocks of size BLOCK_SIZE x BLOCK_SIZE
                for (k = kk; k < kk + BLOCK_SIZE && k < N; k++) {
                    for (j = jj; j < jj + BLOCK_SIZE && j < N; j++) {
                        for (i = ii; i < ii + BLOCK_SIZE && i < N; i++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}



int main() {
    int EventSet = PAPI_NULL;
    long long values[6]; // Array to store event counts
    long long s,e;
    int retval;

    // Allocate memory for matrices
    double** A = (double**)malloc(N * sizeof(double*));
    double** B = (double**)malloc(N * sizeof(double*));
    double** C = (double**)malloc(N * sizeof(double*));
    
    for (long long int i = 0; i < N; i++) {
        A[i] = (double*)malloc(N * sizeof(double));
        B[i] = (double*)malloc(N * sizeof(double));
        C[i] = (double*)malloc(N * sizeof(double));
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
    /*if (PAPI_add_event(EventSet, PAPI_TOT_INS) != PAPI_OK) {
        printf("PAPI_TOT_INS event addition failed\n");
        return 1;
    }*/

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
    initialize_matrices(A, B, C);

    s = PAPI_get_real_usec();

    // Perform matrix multiplication
    // You can swap this function call with others as needed
    matmul_blocked_ikj(A, B, C);

    e = PAPI_get_real_usec();

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

    // Print the results
    /*printf("Matrix A:\n");
    print_matrix(A);
    
    printf("Matrix B:\n");
    print_matrix(B);
    
    printf("Result (C = A * B):\n");
    print_matrix(C);
    */
   
    // Print event counts
    printf("Wallclock_Time: %lld ms\n", ((e-s)/1000));
    printf("PAPI_L1_TCM: %lld\n", values[1]);
    printf("PAPI_L2_TCM: %lld\n", values[2]);
    printf("PAPI_L3_TCM: %lld\n", values[3]);

    // Clean up PAPI
    PAPI_cleanup_eventset(EventSet);
    PAPI_destroy_eventset(&EventSet);
    PAPI_shutdown();

    // Free allocated memory
    for (long long int i = 0; i < N; i++) {
        free(A[i]);
        free(B[i]);
        free(C[i]);
    }
    free(A);
    free(B);
    free(C);
    
    return 0;
}
