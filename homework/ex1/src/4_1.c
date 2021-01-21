#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4
#define MAX(i, j) (i >= j ? i : j)
#define MIN(i, j) (i < j ? i : j)

int main(int argc, char const *argv[]) {
    double A[11][11], _A[11][11];
    int i, j, k;
    for (i = 0; i <= 10; i++)
        for (j = 0; j <= 10; j++) A[i][j] = _A[i][j] = i * i + j * j;

    for (i = 2; i <= 10; i++)
        for (j = i; j <= 10; j++) _A[i][j] = (_A[i][j - 1] + _A[i - 1][j]) / 2;

    /* Parallelize */
    omp_set_num_threads(NUM_THREADS);
    for (k = 4; k <= 20; k++) {
#pragma omp parallel for private(i, j)
        for (i = MAX(2, k - 10); i <= MIN(k / 2, 10); i++) {
            j = k - i;
            A[i][j] = (A[i][j - 1] + A[i - 1][j]) / 2;
            if (A[i][j] != _A[i][j]) printf("ERROR in (%d, %d)", i, j);
        }
    }

    printf("Matrix _A:\n");
    for (i = 2; i <= 10; i++) {
        for (j = 2; j <= 10; j++) printf("%lf\t", _A[i][j]);
        printf("\n");
    }

    printf("Matrix A:\n");
    for (i = 2; i <= 10; i++) {
        for (j = 2; j <= 10; j++) printf("%lf\t", A[i][j]);
        printf("\n");
    }

    return 0;
}
