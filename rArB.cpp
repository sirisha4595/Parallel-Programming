#include <iostream>
#include <cmath>
using namespace std;

int *
init_block (int myrank, int **matrix,
int blocks, int block_size)
{
	int *block_M;
	int row, col;

	block_M = (int *) malloc(pow(block_size, 2) * sizeof(int));
	row = myrank / blocks;
	col = myrank % blocks;

	for (int i = row * block_size; i < (row * block_size) + block_size; i++) {
		for (int j = col * block_size; j < (col * block_size) + block_size; j++) {
			*block_M = matrix[i][j];
			block_M++;
       		}
    	}

    	return block_M;
}

static void 
MM_rArB(int n, int p)
{
	int blocks = sqrt(p);
	int block_size = (n / sqrt_p);
	int **sub_matrix_X, int **sub_matrix_Y, int **sub_matrix_Z;	

	MPI_Status x_status[p];
	MPI_Status y_status[p];
	MPI_Status z_status[p];
	
	sub_matrix_X = (int *) malloc(sizeof(int) * pow(block_size, 2));
	sub_matrix_Y = (int *) malloc(sizeof(int) * pow(block_size, 2));
	sub_matrix_Z = (int *) malloc(sizeof(int) * pow(block_size, 2));
	
	if (myrank == 0) {
		for (int i = 1; i < blocks * blocks; i++) {
			int *block_X = init_block(i, X, blocks, block_size);
			MPI_Send(&block_X, block_size * block_size, MPI_INT, i, x_tag, MPI_COMM_WORLD);
		
			int *block_Y = init_block(i, Y, blocks, block_size)
			MPI_Send(&block_Y, block_size * block_size, MPI_INT, i, y_tag, MPI_COMM_WORLD);
		
			int *block_Z = init_block(i, Z, blocks, block_size)
			MPI_Send(&block_Z, block_size * block_size, MPI_INT, i, z_tag, MPI_COMM_WORLD);

			free(block_X);
			free(block_Y);
			free(block_Z);
		}
	} else if (myrank != 0) {
		MPI_Recv(&(sub_matrix_X[0][0]), pow(block_size, 2), MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_status[myrank]);
		MPI_Recv(&(sub_matrix_Y[0][0]), pow(block_size, 2), MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_status[myrank]);
		MPI_Recv(&(sub_matrix_Z[0][0]), pow(block_size, 2), MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_status[myrank]);
	}
	
