#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char const *argv[]) {
    int size, rank;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int group = 2, sum = rank, temp, process_comm;
    while (group <= size) {
        if (rank % group < group / 2)
            process_comm = rank + group / 2;
        else
            process_comm = rank - group / 2;

        MPI_Sendrecv(&sum, 1, MPI_INT, process_comm, 0, &temp, 1, MPI_INT,
                     process_comm, 0, MPI_COMM_WORLD, &status);

        sum += temp;
        group *= 2;
    }

    printf("rank: %d\tsum:%d\n", rank, sum);

    MPI_Finalize();
    return 0;
}
