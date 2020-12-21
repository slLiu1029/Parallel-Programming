#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define INDEX(i, j, n) ((i) * (n) + (j))
#define N 8

void initialize(int rank, double* A, double* B, double* A1, double* B1,
                int block_dim, int dim) {
    int block_row = rank / block_dim, block_col = rank % block_dim;
    int row, col;

    int i, j;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            row = block_row * dim + i;
            col = block_col * dim + j;

            A[INDEX(i, j, dim)] = A1[INDEX(row, col, N)];
            B[INDEX(i, j, dim)] = 0;
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

int block_row_col_to_rank(int block_row, int block_col, int blocks_dim) {
    return block_row * blocks_dim + block_col;
}

void communicate_with_neighbor(int rank, int size, MPI_Status status, double* A,
                               int block_dim, int dim, double* up, double* down,
                               double* left, double* right) {
    MPI_Sendrecv(A, dim, MPI_DOUBLE, (rank - block_dim + size) % size, 0, down,
                 dim, MPI_DOUBLE, (rank + block_dim) % size, 0, MPI_COMM_WORLD,
                 &status);
    MPI_Sendrecv(A + (dim - 1) * dim, dim, MPI_DOUBLE,
                 (rank + block_dim) % size, 0, up, dim, MPI_DOUBLE,
                 (rank - block_dim + size) % size, 0, MPI_COMM_WORLD, &status);
    int i;
    double left_most[dim], right_most[dim];
    for (i = 0; i < dim; i++) {
        left_most[i] = A[INDEX(i, 0, dim)];
        right_most[i] = A[INDEX(i, dim - 1, dim)];
    }
    MPI_Sendrecv(left_most, dim, MPI_DOUBLE, (rank - 1 + size) % size, 0, right,
                 dim, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD,
                 &status);
    MPI_Sendrecv(right_most, dim, MPI_DOUBLE, (rank + 1) % size, 0, left, dim,
                 MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD,
                 &status);
}

void calculate(int rank, double* A, double* B, int dim, int block_dim,
               double* up, double* down, double* left, double* right) {
    int block_row = rank / block_dim, block_col = rank % block_dim;
    int row, col;
    double up_num, down_num, left_num, right_num;

    int i, j;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            row = block_row * dim + i;
            col = block_col * dim + j;
            if (row <= 0 || row >= N - 1 || col <= 0 || col >= N - 1) continue;

            up_num = !i ? up[j] : A[INDEX(i - 1, j, dim)];
            down_num = i == dim - 1 ? down[j] : A[INDEX(i + 1, j, dim)];
            left_num = !j ? left[i] : A[INDEX(i, j - 1, dim)];
            right_num = j == dim - 1 ? right[i] : A[INDEX(i, j + 1, dim)];

            B[INDEX(i, j, dim)] =
                (up_num + down_num + left_num + right_num) / 4;
        }
}

void print_matrix(double* A) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) printf("%lf\t", A[INDEX(i, j, N)]);
        printf("\n");
    }
}

void print_row_matrix(double* row_matrix, int dim) {
    int i, j, k, start;
    for (i = 0; i < dim; i++) {
        start = i * dim;
        for (j = start; j < dim * N; j += dim * dim) {
            for (k = j; k < j + dim; k++) printf("%lf\t", row_matrix[k]);
        }
        printf("\n");
    }
}

void print_all_blocks_in_order(int rank, int size, double* C, int dim,
                               int blocks_dim) {
    int print_permission = 0;
    int block_row = rank / blocks_dim, block_col = rank % blocks_dim;
    double row_matrix[dim * N];

    if (block_col != 0) {
        int dest = block_row_col_to_rank(block_row, 0, blocks_dim);
        MPI_Send(C, dim * dim, MPI_INT, dest, 0, MPI_COMM_WORLD);
    } else {
        int i;
        for (i = 0; i < dim * dim; i++) row_matrix[i] = C[i];
        for (i = 1; i < blocks_dim; i++) {
            int source = block_row_col_to_rank(block_row, i, blocks_dim);
            MPI_Recv(row_matrix + i * dim * dim, dim * dim, MPI_INT, source, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (!rank) {
        printf("A times B calculated by fox:\n");
        print_permission = 1;
        print_row_matrix(row_matrix, dim);
        int dest = block_row_col_to_rank(1, 0, blocks_dim);
        if (dest < size)
            MPI_Send(&print_permission, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
    } else if (!block_col) {
        int source = block_row_col_to_rank(block_row - 1, 0, blocks_dim);
        MPI_Recv(&print_permission, 1, MPI_INT, source, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        print_row_matrix(row_matrix, dim);
        int dest = block_row_col_to_rank(block_row + 1, 0, blocks_dim);
        if (dest < size)
            MPI_Send(&print_permission, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
    }
}

int check_correctness(int rank, int block_dim, int dim, double* B, double* B1) {
    int is_right = 1;
    int block_row = rank / block_dim, block_col = rank % block_dim;
    int row, col;

    int i, j;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            row = block_row * dim + i;
            col = block_col * dim + j;
            if (B[INDEX(i, j, dim)] != B1[INDEX(row, col, N)]) is_right = 0;
        }

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

    int block_dim = (int)sqrt((double)size), dim = N / block_dim;
    double A[dim * dim], B[dim * dim];
    initialize(rank, A, B, A1, B1, block_dim, dim);

    double up[dim], down[dim], left[dim], right[dim];
    communicate_with_neighbor(rank, size, status, A, block_dim, dim, up, down,
                              left, right);
    calculate(rank, A, B, dim, block_dim, up, down, left, right);

    if (!rank) {
        printf("Matrix A:\n");
        print_matrix(A1);
        printf("Matrix B should be:\n");
        print_matrix(B1);
    }

    print_all_blocks_in_order(rank, size, B, dim, block_dim);

    int is_right = check_correctness(rank, block_dim, dim, B, B1);
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
