#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char const* argv[]) {
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char node_name[20];
    int len;
    int* node_nums = (int*)malloc(sizeof(int) * size);
    MPI_Get_processor_name(node_name, &len);
    int node_num = (int)node_name[len - 1] - '0';
    MPI_Gather(&node_num, 1, MPI_INT, node_nums, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* new_ranks = (int*)malloc(sizeof(int) * size);
    int process_nums_in_node[100];
    int i;
    for (i = 0; i < 100; i++) process_nums_in_node[i] = 0;
    if (!rank)
        for (i = 0; i < size; i++) {
            int node_index = node_nums[i];
            new_ranks[i] = process_nums_in_node[node_index];
            process_nums_in_node[node_index]++;
        }

    MPI_Bcast(new_ranks, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Comm my_world, split_world;
    MPI_Comm_dup(MPI_COMM_WORLD, &my_world);
    MPI_Comm_split(my_world, node_num, new_ranks[rank], &split_world);

    char message[15];
    if (!rank) {
        strcpy(message, "This is Sili!");
        for (i = 1; i < size; i++)
            if (!new_ranks[i])
                MPI_Send(message, 13, MPI_CHAR, i, 0, MPI_COMM_WORLD);
    }

    if (rank && !new_ranks[rank])
        MPI_Recv(&message, 13, MPI_CHAR, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

    int group_size, group_rank;
    MPI_Comm_rank(split_world, &group_rank);
    MPI_Comm_size(split_world, &group_size);
    MPI_Bcast(message, 13, MPI_CHAR, 0, split_world);

    printf("old rank: %d\tnode No. : %d\tnew rank: %d\tmessage: %s\n", rank,
           node_num, new_ranks[rank], message);

    free(node_nums);
    free(new_ranks);

    MPI_Finalize();

    return 0;
}
