#include <math.h>
#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 3
#define min(a, b) (a < b ? a : b)

int main(int argc, char const *argv[]) {
    /* Initialize the arrays */
    int A[20], B[20];
    int A1[20], B1[20];
    int i, j, k;
    for (i = 0; i < 20; i++) {
        A[i] = A1[i] = i * i;
        B[i] = B1[i] = i + 1;
    }

    /* Serial computing */
    for (k = 1; k <= 16; k += 5)
        for (i = k; i <= min(16, i + 4); i++) A1[i + 3] = A1[i] + B1[i];

    /* Parallel computing */
    omp_set_num_threads(NUM_THREADS);
    for (i = 1; i <= 16; i += 3)
#pragma omp parallel for private(j)
        for (j = i; j <= i + 2; j++)
            if (j <= 16) {
                A[j + 3] = A[j] + B[j];
                if (A[j + 3] != A1[j + 3]) printf("ERROR in %d\n", j);
            }

    return 0;
}