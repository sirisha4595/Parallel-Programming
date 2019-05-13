#include "rArB.h"
using namespace std;

static int fastrand()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed>>16) & 0x7FFF; 
}

static void
create_2d_array(int ***array, int n, int m)
{
    int *p = (int *)malloc(n*m*sizeof(int));
    if(!p) {
        std::cout << "Malloc Failed !!!\n";
        return;
    }

    (*array) = (int**)malloc(n*sizeof(int*));
    if (!(*array)) {
       free(p);
       std::cout << "Malloc Failed !!!\n";
       return;
    }

    for (int i=0; i<n; i++) {
       (*array)[i] = &(p[i*m]);
    }
    return;
}

static void
delete_2d_array(int ***array)
{
    free(&((*array)[0][0]));
    free(*array);
    return;
}

static void
fillMatrix(int **arr, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            arr[i][j] = fastrand();        
}

static void
init_sub_matrix(int **arr, int m, int n)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j< m; j++)
            arr[i][j] = 0;
}

static void
Matrix_Multiply(int **X, int **Y, int **Z, int x_row, int x_col, 
int y_row, int y_col, int z_row, int z_col, int n)
{
    for (int i = 0; i < n; i++)
        for (int k = 0; k<n; k++)
            for (int j = 0; j < n; j++)
                Z[z_row + i][z_col + j] += X[x_row + i][x_col + k] * Y[y_row + k][y_col + j];
            
}

static void
update_result(int *temp, int rank, int **Z, int n , int p)
{
    int blocks = std::sqrt(p);
    int block_size = n / blocks;
 
    int row = rank / blocks;
    int col = rank % blocks;
    int *trav = temp;

    for (int i = row * block_size; i < (row * block_size) + block_size; i++) {
        for (int j = col * block_size; j < (col * block_size) + block_size; j++) {
            Z[i][j] += *trav;
            trav++;
        }
    }
}

int * get_buff_copy(int rank,int  **arr, int size, int procs) {
    int blocks = std::sqrt(procs);
    int block_size = size / blocks;

    int * res = (int *)malloc(block_size*block_size*sizeof(int));
    int *trav = res;
    int row = rank / blocks;
    int col = rank % blocks;

    for (int i = row * block_size; i < row * block_size + block_size; ++i) {
        for (int j = col * block_size; j < col * block_size + block_size; ++j) {
            *trav = arr[i][j];
            trav++;
        }
    }

    return res;
}

void
MM_rArB(int n, int p, int rank) 
{
    int blocks = std::sqrt(p);
    int block_size = n / blocks;

    int **sub_matrix_X = NULL;
    int **sub_matrix_Y = NULL;
    int **sub_matrix_Z = NULL;

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

    create_2d_array(&sub_matrix_X, block_size, block_size);
    create_2d_array(&sub_matrix_Y, block_size, block_size);
    create_2d_array(&sub_matrix_Z, block_size, block_size);
    init_sub_matrix(sub_matrix_Z, block_size, block_size);

    if (rank == 0) {
	for (int i = 0; i < block_size; ++i) {
	   for (int j = 0; j < block_size; ++j) {
               sub_matrix_X[i][j] = *(*X + i * n + j);
	       sub_matrix_Y[i][j] = *(*Y + i * n + j);
	       sub_matrix_Z[i][j] = *(*Z + i * block_size + j);
	   }
	}

        for (int i = 1; i < blocks * blocks; ++i) {
            int *sub_X = get_buff_copy(i, X, n, p);
            int *sub_Y = get_buff_copy(i, Y, n, p);
            int *sub_Z = get_buff_copy(i, Z, n, p);

            MPI_Send(sub_X, block_size * block_size, MPI_INT, i, x_tag, MPI_COMM_WORLD);
	    MPI_Send(sub_Y, block_size * block_size, MPI_INT, i, y_tag, MPI_COMM_WORLD);
	    MPI_Send(sub_Z, block_size * block_size, MPI_INT, i, z_tag, MPI_COMM_WORLD);

            free(sub_X);
            free(sub_Y);
            free(sub_Z);
        }
    } 

    if (rank > 0) {
        MPI_Recv(&(sub_matrix_X[0][0]), block_size * block_size, MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_status[rank]);
        MPI_Recv(&(sub_matrix_Y[0][0]), block_size * block_size, MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_status[rank]);
        MPI_Recv(&(sub_matrix_Z[0][0]), block_size * block_size, MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_status[rank]);
    }

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

    MPI_Sendrecv_replace(&(sub_matrix_X[0][0]), block_size * block_size, MPI_INT, send_xproc, rs_tag, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);

    MPI_Sendrecv_replace(&(sub_matrix_Y[0][0]), block_size * block_size, MPI_INT, send_yproc, rs_tag, rec_yproc, rs_tag, MPI_COMM_WORLD, &lu_status);

    for (int l = 0; l < blocks; l++) {
        Matrix_Multiply(sub_matrix_X, sub_matrix_Y, sub_matrix_Z,
                        0, 0, 0, 0, 0, 0, block_size);

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

        MPI_Sendrecv_replace(&(sub_matrix_X[0][0]), block_size * block_size, MPI_INT,
                                 send_xproc, rs_tag, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);
        MPI_Sendrecv_replace(&(sub_matrix_Y[0][0]), block_size * block_size, MPI_INT,
                                 send_yproc, rs_tag, rec_yproc, rs_tag, MPI_COMM_WORLD, &lu_status);
    }

    if (rank != 0)
	MPI_Send(&(sub_matrix_Z[0][0]), block_size * block_size, MPI_INT, 0, merge_tag, MPI_COMM_WORLD);
    if (rank == 0) {
	update_result(&(sub_matrix_Z[0][0]), 0, Z, n , p);
        for (int i = 1; i < blocks * blocks; ++i) {
            int *temp = (int *)malloc(block_size*block_size*sizeof(int));
            MPI_Recv(temp, block_size * block_size, MPI_INT, i, merge_tag, MPI_COMM_WORLD, &merge_status[i]);
            update_result(temp, i, Z, n , p);
            free(temp);
        }
    }

    delete_2d_array(&sub_matrix_X);
    delete_2d_array(&sub_matrix_Y);
    delete_2d_array(&sub_matrix_Z);
    
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
        create_2d_array(&X, *n, *n);
        fillMatrix(X, *n);
        create_2d_array(&Y, *n, *n);
        fillMatrix(Y, *n);
        create_2d_array(&Z, *n, *n);
        init_sub_matrix(Z, *n, *n);
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

    init_var(argc, argv, &n, myrank, processors);
    if (myrank == 0)
	start = std::chrono::high_resolution_clock::now();

    MM_rArB(n, processors, myrank);

    if (myrank == 0) {
	finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Elapsed time: " << elapsed.count() << " s"<< endl;

	delete_2d_array(&X);
	delete_2d_array(&Y);
	delete_2d_array(&Z);
    }

    MPI_Finalize();
    
    return 0;
}
