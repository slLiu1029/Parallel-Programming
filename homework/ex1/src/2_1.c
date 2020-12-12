#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4
#define M 15
#define N 20
#define C 2

int main(int argc, char const *argv[]) {
    /* Initialize the array */
    int A[M + 2][N + 2], A1[M + 2][N + 2], i, j, raise_errors = 0;
    for (i = 0; i <= M + 1; i++)
        for (j = 0; j <= N + 1; j++) {
            A[i][j] = i * i + j * j;
            A1[i][j] = i * i + j * j;
        }

    for (i = 1; i <= M; i++)
        for (j = 1; j <= N; j++) A1[i + 1][j + 1] = A1[i][j] + C;

    /* Parallel program */
    for (i = 1; i <= M; i++) {
        omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(j)
        for (j = 1; j <= N; j++) {
            A[i + 1][j + 1] = A[i][j] + C;
            if (A[i + 1][j + 1] != A1[i + 1][j + 1]) raise_errors = 1;
        }
    }

    /* Check if there exists any error */
    if (raise_errors)
        printf("Errors!");
    else
        printf("No errors!");

    return 0;
}