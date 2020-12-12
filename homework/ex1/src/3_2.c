#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 4

/* Parallel computing */
void parallel(int x, int y, int z, int* A, int* B, int* C, int D[][51]) {
    int i, j;
    x = y * 2;
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel sections private(i, j)
    {
#pragma omp section
        for (i = 1; i <= 100; i++) {
            C[i] = B[i] + x;
            A[i] = C[i - 1] + z;
            C[i + 1] = A[i] * B[i];
        }

#pragma omp section
#pragma omp parallel for private(i, j)
        for (i = 1; i <= 100; i++)
            for (j = 1; j <= 50; j++) D[i][j] = D[i][j - 1] + x;
    }

    z = y + 4;
}

/* Serial computing */
void serial(int x, int y, int z, int* A, int* B, int* C, int D[][51]) {
    int i, j;
    x = y * 2;
    for (i = 1; i <= 100; i++) {
        C[i] = B[i] + x;
        A[i] = C[i - 1] + z;
        C[i + 1] = A[i] * B[i];
        for (j = 1; j <= 50; j++) D[i][j] = D[i][j - 1] + x;
    }
    z = y + 4;
}

/* Check if there is any difference between two one-dimensional arrays */
int check_correctness1(int* A, int* A1, int n) {
    int i, raise_error = 0;
    for (i = 0; i < n; i++)
        if (A[i] != A1[i]) raise_error = 1;

    return raise_error;
}

/* Check if there is any difference between two two-dimensional arrays */
int check_correctness2(int A[][51], int A1[][51], int n) {
    int i, j, raise_error = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < 101; j++)
            if (A[i][j] != A1[i][j]) raise_error = 1;

    return raise_error;
}

int main(int argc, char const* argv[]) {
    /* Initialize the array */
    int x, y = 3, z;
    int A[101], B[101], C[102], D[101][51];
    int A1[101], B1[101], C1[102], D1[101][51];
    int i, j;
    for (i = 0; i < 101; i++) {
        A[i] = A1[i] = i;
        B[i] = B1[i] = i * i;
        for (j = 0; j < 51; j++) D[i][j] = D1[i][j] = i * i + j * j + i * j;
    }
    for (i = 0; i < 102; i++) C[i] = C1[i] = i * i * i;

    /* Calculate arrays with serial and parallel computing */
    serial(x, y, z, A1, B1, C1, D1);
    parallel(x, y, z, A, B, C, D);

    /* Check if there is any error */
    int A_error = check_correctness1(A, A1, 101);
    int B_error = check_correctness1(B, B1, 101);
    int C_error = check_correctness1(C, C1, 102);
    int D_error = check_correctness2(D, D1, 101);

    if (A_error) printf("ERRORS in A!");
    if (B_error) printf("ERRORS in B!");
    if (C_error) printf("ERRORS in C!");
    if (D_error) printf("ERRORS in D!");

    return 0;
}