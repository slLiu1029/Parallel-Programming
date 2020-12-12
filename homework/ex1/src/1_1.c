#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4

int main(int argc, char const *argv[]) {
    /* Initialize the array */
    int A[12][12], A_copy[12][12], i, j, raise_errors = 0;
    for (i = 0; i <= 11; i++)
        for (j = 0; j <= 11; j++) {
            A[i][j] = i * i + j * j;
            A_copy[i][j] = i * i + j * j;
        }

    for (i = 2; i <= 10; i++)
        for (j = 2; j <= 10; j++)
            A_copy[i][j] = (A_copy[i - 1][j - 1] + A_copy[i + 1][j + 1]) / 2;

    /* Parallel program */
    for (i = 2; i <= 10; i++) {
        omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(j)
        for (j = 2; j <= 10; j++) {
            A[i][j] = (A[i - 1][j - 1] + A[i + 1][j + 1]) / 2;
            if (A[i][j] != A_copy[i][j]) raise_errors = 1;
        }
    }

    /* Check if there exists any error */
    if (raise_errors)
        printf("Errors!");
    else
        printf("No errors!");

    return 0;
}
