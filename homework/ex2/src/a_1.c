#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char const* argv[]) {
    // 初始化
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
    MPI_Comm myworld, split_world;
    MPI_Comm_dup(MPI_COMM_WORLD, &myworld);
    MPI_Comm_split(myworld, node_num, new_ranks[rank], &split_world);

    int new_rank;
    MPI_Comm_rank(split_world, &new_rank);
    printf("old rank: %d\tnode No. : %d\tnew rank: %d\n", rank, node_num,
           new_rank);

    free(node_nums);
    free(new_ranks);

    MPI_Finalize();

    return 0;
}
