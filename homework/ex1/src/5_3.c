#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4
#define N 1000

int A[N][N], C[N][N], D[N][N];
int A1[N][N], C1[N][N], D1[N][N];

int main(int argc, char const *argv[]) {
    /* Initialize the arrays */
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            A[i][j] = A1[i][j] = i + j;
            C[i][j] = C1[i][j] = i - j;
            D[i][j] = D1[i][j] = i * j;
        }

    /* Serial computing */
    for (i = 1; i <= 100; i++)
        for (j = 1; j <= 100; j++) {
            A1[3 * i + 2 * j][2 * j] = C1[i][j] * 2;
            if (i - j + 6 >= 0) D1[i][j] = A1[i - j + 6][i + j];
        }

    /* Parallel computing */
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(i, j)
    for (i = 1; i <= 100; i++)
#pragma omp parallel for private(j)
        for (j = 1; j <= 100; j++) {
            A[3 * i + 2 * j][2 * j] = C[i][j] * 2;
            if (A[3 * i + 2 * j][2 * j] != A1[3 * i + 2 * j][2 * j])
                printf("ERROR in A[%d][%d]\n", i, j);
            if (i - j + 6 >= 0) {
                D[i][j] = A[i - j + 6][i + j];
                if (D[i][j] != D1[i][j]) printf("ERROR in D[%d][%d]\n", i, j);
            }
        }

    return 0;
}
