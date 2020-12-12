#include <math.h>
#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4

int main(int argc, char const *argv[]) {
    /* Initialize the arrays */
    float A[101], B[101], C[101], D[101];
    float A1[101], B1[101], C1[101], D1[101];
    int i;
    for (i = 0; i < 101; i++) {
        A[i] = A1[i] = i;
        B[i] = B1[i] = (i + 1) / 2;
        C[i] = C1[i] = i * i * i;
        D[i] = D1[i] = i * i + 3 * i + 1;
    }

    /* Serial computing */
    for (i = 1; i <= 100; i++) {
        A1[i] = A1[i] + B1[i - 1];
        B1[i] = C1[i - 1] * 2;
        C1[i] = 1 / B1[i];
        D1[i] = C1[i] * C1[i];
    }

    /* Parallel computing */
    for (i = 1; i <= 100; i++) {
        B[i] = C[i - 1] * 2;
        C[i] = 1 / B[i];
        if (B[i] != B1[i]) printf("ERROR in B[%d]\n", i);
        if (C[i] != C1[i]) printf("ERROR in C[%d]\n", i);
    }

    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel sections
    {
#pragma omp section
#pragma omp parallel for private(i)
        for (i = 1; i <= 100; i++) {
            A[i] = A[i] + B[i - 1];
            if (A[i] != A1[i]) printf("ERROR in A[%d]\n", i);
        }

#pragma omp section
#pragma omp parallel for private(i)
        for (i = 1; i <= 100; i++) {
            D[i] = C[i] * C[i];
            if (D[i] != D1[i]) printf("ERROR in D[%d]\n", i);
        }
    }

    return 0;
}