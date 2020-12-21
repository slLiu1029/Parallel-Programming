#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void my_MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                     void *recvbuf, int recvcount, MPI_Datatype recvtype,
                     MPI_Comm comm) {
    int i, j, index, size, rank;

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (i = 0; i < size; i++) {
        index = i * sendcount * sizeof(sendtype);
        if (rank == i) {
            MPI_Sendrecv(sendbuf + index, sendcount, sendtype, rank, 1,
                         recvbuf + index, recvcount, recvtype, rank, 1,
                         MPI_COMM_WORLD, &status);
        } else {
            MPI_Send(sendbuf + index, sendcount, sendtype, i, 1,
                     MPI_COMM_WORLD);
        }
    }

    for (i = 0; i < size; i++) {
        index = i * sendcount * sizeof(recvtype);
        if (rank != i) {
            MPI_Recv(recvbuf + index, recvcount, recvtype, i, 1, MPI_COMM_WORLD,
                     &status);
        }
    }
}

int main(int argc, char const *argv[]) {
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *sendbuf = (int *)malloc(sizeof(int) * size * 3);
    int *recvbuf = (int *)malloc(sizeof(int) * size * 3);
    int *sendbuf1 = (int *)malloc(sizeof(int) * size * 3);
    int *recvbuf1 = (int *)malloc(sizeof(int) * size * 3);
    int i;
    for (i = 0; i < 3 * size; i++) sendbuf1[i] = sendbuf[i] = i + rank;

    MPI_Barrier(MPI_COMM_WORLD);
    double begin_time = MPI_Wtime();
    my_MPI_Alltoall(sendbuf1, 3, MPI_INT, recvbuf1, 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double time_of_my_func = end_time - begin_time;

    MPI_Barrier(MPI_COMM_WORLD);
    begin_time = MPI_Wtime();
    MPI_Alltoall(sendbuf, 3, MPI_INT, recvbuf, 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    double time_of_MPI_func = end_time - begin_time;

    int is_right = 1;
    for (i = 0; i < 3 * size; i++)
        if (recvbuf1[i] != recvbuf[i]) is_right = 0;
    if (is_right == 0) printf("rank %d has errors\n", rank);

    MPI_Barrier(MPI_COMM_WORLD);

    int is_all_right;
    MPI_Reduce(&is_right, &is_all_right, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
    if (!rank) {
        if (is_all_right == size) {
            printf("Right!\n");
            printf("My function time: %lf\n", time_of_my_func);
            printf("MPI function time: %lf\n", time_of_MPI_func);
        } else
            printf("ERROR!\n");
    }

    MPI_Finalize();
    return 0;
}
