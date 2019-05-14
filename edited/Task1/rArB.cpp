#include "rArB.h"
using namespace std;

static int fastrand()
{
    int ret;
    
    g_seed = (343243 * g_seed + 31231233);
    ret = (g_seed >> 16) & 0x7FFF; 

    return ret;
}

static void
new_matrix(int ***array, int n, int flag)
{
    int *new_matrix = (int *) malloc(pow(n, 2) * sizeof(int));
    
    (*array) = (int **) malloc(n * sizeof(int *));
    
    for (int i = 0; i < n; i++)
       (*array)[i] = &(new_matrix[i * n]);

    if (flag == 1)
	init_matrix(*array, n, 0);
    return;
}

static void
free_matrix(int ***arr)
{
    free(&(*(arr)[0][0]));
    free(*arr);
    return;
}

static void
init_matrix(int **arr, int n, int flag)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
	    if (flag == 1)
            	arr[i][j] = fastrand();
	    else
		arr[i][j] = 0;
}

static void
MMultiply(int **X, int **Y, int **Z,
int x_row, int x_col, int y_row, int y_col,
int z_row, int z_col, int n)
{
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                Z[z_row + i][z_col + j] += X[x_row + i][x_col + k] \
					   * Y[y_row + k][y_col + j];
            
}

static void
updateZ(int *temp, int rank, int **Z,
int block_size, int blocks)
{
    int row = rank / blocks;
    int col = rank % blocks;
    int count = 0;

    for (int i = row * block_size; i < (row * block_size) + block_size; i++)
        for (int j = col * block_size; j < (col * block_size) + block_size; j++)
            Z[i][j] += *(temp + count++);
}

static int *
clone(int rank, int  **arr,
int block_size, int blocks)
{
    int row = rank / blocks;
    int col = rank % blocks;
    int count = 0;

    int *ret = (int *) malloc(pow(block_size, 2) * sizeof(int));

    for (int i = block_size * row; i < (block_size * row) + block_size; i++)
        for (int j = block_size * col; j < (block_size * col) + block_size; j++)
		*(ret + count++) = arr[i][j];

    return ret;
}

static void
MM_rArB(int rank, int n, int p) 
{
    int blocks = std::sqrt(p);
    int block_size = n / blocks;

    int **sub_X = NULL;
    int **sub_Y = NULL;
    int **sub_Z = NULL;

    MPI_Status lu_status;
    MPI_Status merge_status[p];
    MPI_Status x_status[p];
    MPI_Status y_status[p];
    MPI_Status z_status[p];

    int x_tag = 01;
    int y_tag = 02;
    int z_tag = 02;
    int rs_tag = 04;
    int merge_tag = 05;

    new_matrix(&sub_X, block_size, 0);
    new_matrix(&sub_Y, block_size, 0);
    new_matrix(&sub_Z, block_size, 1);

    if (rank == 0) {
	for (int i = 0; i < block_size; ++i) {
	   for (int j = 0; j < block_size; ++j) {
               sub_X[i][j] = *(*X + i * n + j);
	       sub_Y[i][j] = *(*Y + i * n + j);
	       sub_Z[i][j] = *(*Z + i * block_size + j);
	   }
	}

        for (int i = 1; i < blocks * blocks; ++i) {
            int *sub_X = clone(i, X, block_size, blocks);
            int *sub_Y = clone(i, Y, block_size, blocks);
            int *sub_Z = clone(i, Z, block_size, blocks);

            MPI_Send(sub_X, block_size * block_size, MPI_INT, i, x_tag, MPI_COMM_WORLD);
	    MPI_Send(sub_Y, block_size * block_size, MPI_INT, i, y_tag, MPI_COMM_WORLD);
	    MPI_Send(sub_Z, block_size * block_size, MPI_INT, i, z_tag, MPI_COMM_WORLD);

            free(sub_X);
            free(sub_Y);
            free(sub_Z);
        }
    } 

    if (rank > 0) {
        MPI_Recv(&(sub_X[0][0]), block_size * block_size, MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_status[rank]);
        MPI_Recv(&(sub_Y[0][0]), block_size * block_size, MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_status[rank]);
        MPI_Recv(&(sub_Z[0][0]), block_size * block_size, MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_status[rank]);
    }
    if (rank > - 1) {
    	int xy_row = rank / blocks; 
    	int xy_col = rank % blocks;

    	int send_xrow = xy_row; 
    	int send_xcol = (blocks + xy_col - xy_row) % blocks; 
    	int send_xproc = (blocks * send_xrow ) + send_xcol;

    	int rec_xrow = xy_row;
    	int rec_xcol = (blocks + xy_col + xy_row) % blocks; 
    	int rec_xproc = (blocks * rec_xrow ) + rec_xcol;

    	int send_yrow = (blocks + xy_row - xy_col) % blocks; 
    	int send_ycol = xy_col; 
    	int send_yproc = (blocks * send_yrow ) + send_ycol;
	
	int rec_yrow = (blocks + xy_row + xy_col) % blocks;
   	int rec_ycol = xy_col; 
   	int rec_yproc = (blocks * rec_yrow ) + rec_ycol;

    	MPI_Sendrecv_replace(&(sub_X[0][0]), block_size * block_size, MPI_INT, send_xproc, rs_tag, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);
    	MPI_Sendrecv_replace(&(sub_Y[0][0]), block_size * block_size, MPI_INT, send_yproc, rs_tag, rec_yproc, rs_tag, MPI_COMM_WORLD, &lu_status);

    	for (int l = 0; l < blocks; l++) {
        	MMultiply(sub_X, sub_Y, sub_Z, 0, 0, 0, 0, 0, 0, block_size);

        	send_xcol = (blocks + xy_col - 1) % blocks; 
        	send_xproc = (blocks * send_xrow ) + send_xcol;

        	rec_xcol = (blocks + xy_col + 1) % blocks; 
        	rec_xproc = (blocks * rec_xrow ) + rec_xcol;

        	send_yrow = (blocks + xy_row - 1) % blocks; 
        	send_ycol = xy_col; 
        	send_yproc = (blocks * send_yrow ) + send_ycol;

        	rec_yrow = (blocks + xy_row + 1) % blocks;
        	rec_ycol = xy_col; 
        	rec_yproc = (blocks * rec_yrow ) + rec_ycol;

       	 	MPI_Sendrecv_replace(&(sub_X[0][0]), block_size * block_size, MPI_INT,
                		     send_xproc, rs_tag, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);
        	MPI_Sendrecv_replace(&(sub_Y[0][0]), block_size * block_size, MPI_INT,
                                     send_yproc, rs_tag, rec_yproc, rs_tag, MPI_COMM_WORLD, &lu_status);
    	}

    	if (rank != 0)
		MPI_Send(&(sub_Z[0][0]), block_size * block_size, MPI_INT, 0, merge_tag, MPI_COMM_WORLD);
    	if (rank == 0) {
		updateZ(&(sub_Z[0][0]), 0, Z, block_size, blocks);
	        for (int i = 1; i < pow(blocks, 2); ++i) {
	            int *tempZ = (int *) malloc(pow(block_size, 2) * sizeof(int));
		    MPI_Recv(tempZ, block_size * block_size, MPI_INT, i, merge_tag, MPI_COMM_WORLD, &merge_status[i]);
		    updateZ(tempZ, i, Z, block_size, blocks);
		    free(tempZ);
        	}
  	  }
    }

    free_matrix(&sub_X);
    free_matrix(&sub_Y);
    free_matrix(&sub_Z);
    
    return;
}

static void 
init_var(int argc, char **argv,
int *n, int &myrank, int &p)
{
    if (argc < 2)
	goto out;

    *n = pow(2, atoi(argv[1]));
    srand(time(NULL));
    g_seed = rand();
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    if (myrank == 0) {
        new_matrix(&X, *n, 0);
        new_matrix(&Y, *n, 0);
        new_matrix(&Z, *n, 1);
    }
    return;
out:
    cout << "Invalid Input: mpirun -n #processors ./rArB matrix_size" << endl;
    exit(0);
}

int main(int argc, char *argv[])
{
    int n, myrank, p;
    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point finish;

    init_var(argc, argv, &n, myrank, p);
    if (myrank == 0)
	start = std::chrono::high_resolution_clock::now();

    MM_rArB(myrank, n, p);

    if (myrank == 0) {
	finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Elapsed time: " << elapsed.count() << " s"<< endl;

	free_matrix(&X);
	free_matrix(&Y);
	free_matrix(&Z);
    }

    MPI_Finalize();
    
    return 0;
}
