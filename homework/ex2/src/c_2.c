#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) {
    int size, rank;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int group = 2, sum = rank, process_comm, temp;
    while (group <= size) {
        if (rank % group == 0) {
            process_comm = rank + group / 2;
            MPI_Recv(&temp, 1, MPI_INT, process_comm, 0, MPI_COMM_WORLD,
                     &status);
            sum += temp;
        } else if (rank % group == group / 2) {
            process_comm = rank - group / 2;
            MPI_Send(&sum, 1, MPI_INT, process_comm, 0, MPI_COMM_WORLD);
        }

        group *= 2;
    }

    group = size;
    while (group >= 2) {
        if (rank % group == 0) {
            process_comm = rank + group / 2;
            MPI_Send(&sum, 1, MPI_INT, process_comm, 0, MPI_COMM_WORLD);
        } else if (rank % group == group / 2) {
            process_comm = rank - group / 2;
            MPI_Recv(&sum, 1, MPI_INT, process_comm, 0, MPI_COMM_WORLD,
                     &status);
        }

        group /= 2;
    }

    printf("rank: %d\tsum: %d\n", rank, sum);

    MPI_Finalize();

    return 0;
}
