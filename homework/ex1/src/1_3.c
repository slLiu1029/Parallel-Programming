#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4

/* Parallel computing */
void parallel(int *A, int *B, int *C) {
    int i;
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(i)
    for (i = 2; i <= 20; i++)
        if (A[i] > 0)
            B[i] = C[i - 1] + 1;
        else
            C[i] = B[i] - 1;
}

/* Serial computing */
void serial(int *A, int *B, int *C) {
    int i;
    for (i = 2; i <= 20; i++)
        if (A[i] > 0)
            B[i] = C[i - 1] + 1;
        else
            C[i] = B[i] - 1;
}

/* Check if there is any difference between two arrays */
int check_correctness(int *A, int *A1) {
    int i, raise_errors = 0;
    for (i = 0; i < 21; i++)
        if (A[i] != A1[i]) raise_errors = 1;
    return raise_errors;
}

int main(int argc, char const *argv[]) {
    int A[21], B[21], C[21], A1[21], B1[21], C1[21], i;
    /* Initialize the arrays */
    for (i = 0; i < 21; i++) {
        if (i < 13)
            A[i] = A1[i] = 1;
        else
            A[i] = A1[i] = -1;
        B[i] = B1[i] = i * i;
        C[i] = C1[i] = 3 * i;
    }

    serial(A1, B1, C1);
    parallel(A, B, C);

    int B_errors = check_correctness(B, B1);
    int C_errors = check_correctness(C, C1);
    if (B_errors || C_errors)
        printf("Errors!");
    else
        printf("No errors!");

    return 0;
}