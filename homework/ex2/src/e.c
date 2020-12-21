#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define P 3
#define MAX_LOOP 2

int main(int argc, char const *argv[]) {
    int rank, size;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int Q = size - P, servers[P], i;
    for (i = 0; i < P; i++) servers[i] = i;

    MPI_Comm server_world;
    MPI_Group world_group, server_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, P, servers, &server_group);
    MPI_Comm_create(MPI_COMM_WORLD, server_group, &server_world);

    srand(rank + time(0));
    int loop = 0;
    while (loop < MAX_LOOP) {
        if (rank >= P) {
            int msg = rand() % 1000, served_by = rank % P;
            printf("rank: %d\tserver: %d\tmessage: %d\n", rank, served_by, msg);
            MPI_Send(&msg, 1, MPI_INT, served_by, 0, MPI_COMM_WORLD);
            double mean_result;
            MPI_Recv(&mean_result, 1, MPI_DOUBLE, served_by, 0, MPI_COMM_WORLD,
                     &status);
            printf("rank: %d\tserver: %d\tmean message: %f\tloop: %d\n", rank,
                   served_by, mean_result, loop);
        } else {
            int recv_msg, msg_sum = 0, msg_num = 0;
            for (i = rank + P; i < size; i += P) {
                MPI_Recv(&recv_msg, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                msg_num++;
                msg_sum += recv_msg;
            }

            int all_msg_sum, all_msg_num;
            MPI_Barrier(server_world);
            MPI_Allreduce(&msg_sum, &all_msg_sum, 1, MPI_INT, MPI_SUM,
                          server_world);
            MPI_Allreduce(&msg_num, &all_msg_num, 1, MPI_INT, MPI_SUM,
                          server_world);
            double msgs_mean = all_msg_sum * 1.0 / all_msg_num;

            for (i = rank + P; i < size; i += P)
                MPI_Send(&msgs_mean, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        loop++;
    }

    MPI_Finalize();

    return 0;
}
