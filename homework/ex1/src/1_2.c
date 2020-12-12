#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 2

int main(int argc, char const *argv[]) {
    /* Initialize the array */
    int A[43], A_copy[43], B[43], i, j, raise_errors = 0;
    for (i = 0; i < 43; i++) {
        A[i] = i * i * i;
        A_copy[i] = i * i * i;
        B[i] = i + i * i;
    }
    for (i = 2; i <= 20; i++) A_copy[2 * i + 2] = A_copy[2 * i - 2] + B[i];

    /* Parallel program */
    omp_set_num_threads(NUM_THREADS);
    for (i = 2; i <= 20; i += 2) {
#pragma omp parallel for private(j)
        for (j = i; j <= i + 1; j++) {
            if (j > 20) continue;
            A[2 * i + 2] = A[2 * i - 2] + B[i];
            if (A[2 * i + 2] != A_copy[2 * i + 2]) raise_errors = 1;
        }
    }

    /* Check if there exists any error */
    if (raise_errors)
        printf("Errors!");
    else
        printf("No errors!");

    return 0;
}
