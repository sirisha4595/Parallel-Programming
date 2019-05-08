#include<stdio.h>
#include<mpi.h>

int main(int argc, char *argv[])
{
	int p, myrank;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	printf("This is process %d out of process %d\n", p, myrank);
	
	MPI_Finalize();
	
	return 0;
}
