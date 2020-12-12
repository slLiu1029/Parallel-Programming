#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4

int main(int argc, char const *argv[]) {
    /* Initialize the arrays */
    double A[1000], B[1000], C[1000], D[1000];
    double A1[1000], B1[1000], C1[1000], D1[1000];
    int i;
    for (i = 0; i < 1000; i++) {
        A[i] = A1[i] = i * i + 1;
        B[i] = B1[i] = i + 10;
        C[i] = C1[i] = 3 * i + 5;
        D[i] = D1[i] = i;
    }

    /* Serial computing */
    for (i = 1; i <= 999; i++) {
        A1[i] = B1[i] + C1[i];
        D1[i] = (A1[i] + A1[999 - i + 1]) / 2;
    }

    /* Parallel computing */
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(i)
    for (i = 1; i <= 499; i++) {
        A[i] = B[i] + C[i];
        if (A[i] != A1[i]) printf("ERROR in A[%d]\n", i);
    }
#pragma omp parallel for private(i)
    for (i = 1; i <= 499; i++) {
        D[i] = (A[i] + A[999 - i + 1]) / 2;
        if (D[i] != D1[i]) printf("ERROR in D[%d]\n", i);
    }
#pragma omp parallel for private(i)
    for (i = 500; i <= 999; i++) {
        A[i] = B[i] + C[i];
        if (A[i] != A1[i]) printf("ERROR in A[%d]\n", i);
    }
#pragma omp parallel for private(i)
    for (i = 500; i <= 999; i++) {
        D[i] = (A[i] + A[999 - i + 1]) / 2;
        if (D[i] != D1[i]) printf("ERROR in D[%d]\n", i);
    }

    return 0;
}