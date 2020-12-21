#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define INDEX(i, j, n) ((i) * (n) + (j))
#define N 8

void matrix_mul(int* A, int* B, int* C, int dim) {
    int i, j, k, sum;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            C[INDEX(i, j, dim)] = 0;
            for (k = 0; k < dim; k++)
                C[INDEX(i, j, dim)] +=
                    A[INDEX(i, k, dim)] * B[INDEX(k, j, dim)];
        }
}

void initialize_matrix_AB(int* A, int* B) {
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) {
            A[INDEX(i, j, N)] = i + j;
            B[INDEX(i, j, N)] = i * j - 2;
        }
}

void initialize_block_matrix(int* A, int* B, int* A1, int* B1, int blocks_dim,
                             int dim, int rank) {
    int row, col;
    int i, j;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            row = (rank / blocks_dim) * dim + i;
            col = (rank % blocks_dim) * dim + j;

            A[INDEX(i, j, dim)] = A1[INDEX(row, col, N)];
            B[INDEX(i, j, dim)] = B1[INDEX(row, col, N)];
        }
}

int block_row_col_to_rank(int block_row, int block_col, int blocks_dim) {
    return block_row * blocks_dim + block_col;
}

int check_size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (sqrt((double)size) - (int)sqrt((double)size) != 0.0) {
        printf("Number of processes must be perfect squre~\n");
        return 0;
    } else
        return 1;
}

void copy(int* A, int* A1, int dim) {
    int i, j;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) A[INDEX(i, j, dim)] = A1[INDEX(i, j, dim)];
}

void fox(int* A1, int* B1, int* C, int rank, int size) {
    int blocks_dim = (int)sqrt((double)size), dim = N / blocks_dim;
    int A[dim * dim], B[dim * dim];
    initialize_block_matrix(A, B, A1, B1, blocks_dim, dim, rank);
    int block_row = rank / blocks_dim, block_col = rank % blocks_dim;

    MPI_Barrier(MPI_COMM_WORLD);
    int loop;
    int i, j;
    int recv_buff[dim * dim], C_temp[dim * dim];
    for (loop = 0; loop < blocks_dim; loop++) {
        if (block_col == (block_row + loop) % blocks_dim) {
            for (j = 0; j < blocks_dim; j++)
                if (j != block_col)
                    MPI_Send(A, dim * dim, MPI_INT,
                             block_row_col_to_rank(block_row, j, blocks_dim), 0,
                             MPI_COMM_WORLD);

            matrix_mul(A, B, C_temp, dim);
        } else {
            MPI_Recv(
                recv_buff, dim * dim, MPI_INT,
                block_row_col_to_rank(
                    block_row, (block_row + loop) % blocks_dim, blocks_dim),
                0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            matrix_mul(recv_buff, B, C_temp, dim);
        }
        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                C[INDEX(i, j, dim)] += C_temp[INDEX(i, j, dim)];

        MPI_Sendrecv(
            B, dim * dim, MPI_INT,
            block_row_col_to_rank((block_row - 1 + blocks_dim) % blocks_dim,
                                  block_col, blocks_dim),
            0, recv_buff, dim * dim, MPI_INT,
            block_row_col_to_rank((block_row + 1) % blocks_dim, block_col,
                                  blocks_dim),
            0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        copy(B, recv_buff, dim);

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int check_correctness(int* C, int* C1, int dim, int rank, int blocks_dim) {
    int row, col;
    int is_right = 1;

    int i, j;
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            row = (rank / blocks_dim) * dim + i;
            col = (rank % blocks_dim) * dim + j;

            if (C[INDEX(i, j, dim)] != C1[INDEX(row, col, N)]) is_right = 0;
        }

    return is_right;
}

void print_matrix(int* A) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) printf("%d\t", A[INDEX(i, j, N)]);
        printf("\n");
    }
}

void print_row_matrix(int* row_matrix, int dim) {
    int i, j, k, start;
    for (i = 0; i < dim; i++) {
        start = i * dim;
        for (j = start; j < dim * N; j += dim * dim) {
            for (k = j; k < j + dim; k++) printf("%d\t", row_matrix[k]);
        }
        printf("\n");
    }
}

void print_all_blocks_in_order(int rank, int size, int* C, int dim,
                               int blocks_dim) {
    int print_permission = 0;
    int block_row = rank / blocks_dim, block_col = rank % blocks_dim;
    int row_matrix[dim * N];

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

int main(int argc, char const* argv[]) {
    int A1[N * N], B1[N * N], C1[N * N];
    initialize_matrix_AB(A1, B1);
    matrix_mul(A1, B1, C1, N);

    int rank, size;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!check_size()) return 0;

    int blocks_dim = (int)sqrt((double)size), dim = N / blocks_dim;
    int C[dim * dim];
    int i;
    for (i = 0; i < dim * dim; i++) C[i] = 0;
    fox(A1, B1, C, rank, size);

    if (!rank) {
        printf("Matrix A:\n");
        print_matrix(A1);
        printf("Matrix B:\n");
        print_matrix(B1);
        printf("A times B should be:\n");
        print_matrix(C1);
    }

    print_all_blocks_in_order(rank, size, C, dim, blocks_dim);

    int is_right = check_correctness(C, C1, dim, rank, blocks_dim);
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
