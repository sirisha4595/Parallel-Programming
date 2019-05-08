#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char *argv[ ])
{
	int p, myrank, v = 121;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	if (myrank == 0) {
		MPI_Send(&v, 1, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD);
		cout << "Process " << myrank << " sent " << v << endl;
	} else if (myrank == 1) {
		MPI_Recv(&v, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		cout << "Process " << myrank << " receieved " << v << endl;
	}

	MPI_Finalize();

	return 0;
}
		
	


