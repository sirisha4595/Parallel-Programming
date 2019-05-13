#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char *argv[ ])
{
	int p, myrank, v = 121;
	int s_tag = 10;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank == 0) {
		MPI_Send(&v, 1, MPI_INT, 1, s_tag, MPI_COMM_WORLD);
		cout << "Process " << myrank << " sent " << v << endl;
	} 
	if (myrank == 1) {
		MPI_Recv(&v, 1, MPI_INT, 0, s_tag, MPI_COMM_WORLD, &status);
		cout << "Process " << myrank << " receieved " << v << endl;
	}

	MPI_Finalize();

	return 0;
}
