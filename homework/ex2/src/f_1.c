#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define INDEX(i, j, n) ((i) * (n) + (j))
#define N 8

void initialize(int rank, double* A, double* B, double* A1, double* B1,
                int rows) {
    int row, col;
    int i, j;
    for (i = 0; i < rows; i++)
        for (j = 0; j < N; j++) {
            row = rank * rows + i;
            col = j;

            A[INDEX(i, j, N)] = A1[INDEX(row, col, N)];
            B[INDEX(i, j, N)] = 0;
        }
}

void initialize1(double* A, double* B) {
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            A[INDEX(i, j, N)] = i * i + j;
            B[INDEX(i, j, N)] = 0;
        }

    for (i = 1; i < N - 1; i++)
        for (j = 1; j < N - 1; j++)
            B[INDEX(i, j, N)] =
                (A[INDEX(i, j + 1, N)] + A[INDEX(i, j - 1, N)] +
                 A[INDEX(i + 1, j, N)] + A[INDEX(i - 1, j, N)]) /
                4;
}

void communicate_with_neighbor(int rank, int size, MPI_Status status, double* A,
                               int rows, double* up, double* low) {
    MPI_Sendrecv(A, N, MPI_DOUBLE, (rank - 1 + size) % size, 0, low, N,
                 MPI_DOUBLE, (rank + 1 + size) % size, 0, MPI_COMM_WORLD,
                 &status);
    MPI_Sendrecv(A + (rows - 1) * N, N, MPI_DOUBLE, (rank + 1 + size) % size, 0,
                 up, N, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD,
                 &status);
}

void calculate(int rank, double* A, double* B, int rows, int size, double* up,
               double* low) {
    int i, j;
    for (j = 1; j < N - 1; j++) {
        if (rank)
            B[INDEX(0, j, N)] = (A[INDEX(0, j + 1, N)] + A[INDEX(0, j - 1, N)] +
                                 A[INDEX(1, j, N)] + up[j]) /
                                4;
        if (rank != size - 1)
            B[INDEX(rows - 1, j, N)] =
                (A[INDEX(rows - 1, j + 1, N)] + A[INDEX(rows - 1, j - 1, N)] +
                 A[INDEX(rows - 2, j, N)] + low[j]) /
                4;
        for (i = 1; i < rows - 1; i++)
            B[INDEX(i, j, N)] =
                (A[INDEX(i, j + 1, N)] + A[INDEX(i, j - 1, N)] +
                 A[INDEX(i + 1, j, N)] + A[INDEX(i - 1, j, N)]) /
                4;
    }
}

void print_matrix(double* A) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) printf("%lf\t", A[INDEX(i, j, N)]);
        printf("\n");
    }
}

void print_row_matrix(double* A, int rows) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < N; j++) printf("%lf\t", A[INDEX(i, j, N)]);
        printf("\n");
    }
}

void print_all_blocks_in_order(int rank, int size, double* B, int rows) {
    int print_permission = 0;
    int block_row = rank * rows;

    if (!rank) {
        printf("Partition with rows calculate result:\n");
        print_permission = 1;
        print_row_matrix(B, rows);
    } else {
        MPI_Recv(&print_permission, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        print_row_matrix(B, rows);
    }
    if (rank + 1 < size)
        MPI_Send(&print_permission, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
}

int check_correctness(int rank, int rows, double* B, double* B1) {
    int is_right = 1;
    int i, j;
    for (i = 0; i < rows; i++)
        for (j = 0; j < N; j++)
            if (B[INDEX(i, j, N)] != B1[INDEX(i + rank * rows, j, N)])
                is_right = 0;

    return is_right;
}

int main(int argc, char const* argv[]) {
    double A1[N * N], B1[N * N];
    initialize1(A1, B1);

    int rank, size;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows_per_block = N / size;
    double A[rows_per_block * N], B[rows_per_block * N];
    initialize(rank, A, B, A1, B1, rows_per_block);

    double up[N], low[N];
    communicate_with_neighbor(rank, size, status, A, rows_per_block, up, low);
    calculate(rank, A, B, rows_per_block, size, up, low);

    if (!rank) {
        printf("Matrix A:\n");
        print_matrix(A1);
        printf("Matrix B should be:\n");
        print_matrix(B1);
    }

    print_all_blocks_in_order(rank, size, B, rows_per_block);

    int is_right = check_correctness(rank, rows_per_block, B, B1);
    int is_all_right;
    MPI_Reduce(&is_right, &is_all_right, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (!rank) {
        if (is_all_right == size)
            printf("Right!\n");
        else
            printf("ERROR\n");
    }

    MPI_Finalize();

    return 0;
}
