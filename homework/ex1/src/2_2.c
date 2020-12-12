#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4
#define N 50

/* Serial computing */
void serial(int* X, int* Y, int A[][101], int* B, int C[][101]) {
    int i, j, k;
    for (i = 1; i <= 100; i++) {
        X[i] = Y[i] + 10;
        for (j = 1; j <= 100; j++) {
            B[j] = A[j][N];
            for (k = 1; k <= 100; k++) A[j + 1][k] = B[j] + C[j][k];
            Y[i + j] = A[j + 1][N];
        }
    }
}

/* Parallel computing */
void parallel(int* X, int* Y, int A[][101], int* B, int C[][101]) {
    int i, j, k;
    for (i = 1; i <= 100; i++) {
        for (j = 1; j <= 100; j++) {
            B[j] = A[j][N];
            omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(k)
            for (k = 1; k <= 100; k++) A[j + 1][k] = B[j] + C[j][k];
        }

        omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(j)
        for (j = 1; j <= 100; j++) Y[i + j] = A[j + 1][N];
    }

    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private(i)
    for (i = 1; i <= 100; i++) X[i] = Y[i] + 10;
}

/* Check if there is any difference between two one-dimensional arrays */
int check_correctness1(int* A, int* A1, int n) {
    int i, raise_error = 0;
    for (i = 0; i < n; i++)
        if (A[i] != A1[i]) raise_error = 1;

    return raise_error;
}

/* Check if there is any difference between two two-dimensional arrays */
int check_correctness2(int A[][101], int A1[][101], int n) {
    int i, j, raise_error = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < 101; j++)
            if (A[i][j] != A1[i][j]) raise_error = 1;

    return raise_error;
}

int main(int argc, char const* argv[]) {
    /* Initialize the arrays */
    int i, j;
    int X[101], Y[201], A[102][101], B[101], C[101][101];
    int X1[101], Y1[201], A1[102][101], B1[101], C1[101][101];
    for (i = 0; i <= 100; i++) {
        X[i] = X1[i] = i;
        B[i] = B1[i] = i * i;
        for (j = 0; j <= 100; j++) C[i][j] = C1[i][j] = i + j + i * j;
    }
    for (i = 0; i <= 200; i++) Y[i] = Y1[i] = i - 2;
    for (i = 0; i <= 101; i++)
        for (j = 0; j <= 100; j++) A[i][j] = A1[i][j] = i * i + j;

    /* Calculate arrays with serial and parallel computing */
    serial(X1, Y1, A1, B1, C1);
    parallel(X, Y, A, B, C);

    /* Check if there is any error */
    int X_error = check_correctness1(X, X1, 101);
    int Y_error = check_correctness1(Y, Y1, 201);
    int B_error = check_correctness1(B, B1, 101);
    int A_error = check_correctness2(A, A1, 102);

    if (X_error) printf("ERRORS in X!");
    if (Y_error) printf("ERRORS in Y!");
    if (A_error) printf("ERRORS in A!");
    if (B_error) printf("ERRORS in B!");

    return 0;
}