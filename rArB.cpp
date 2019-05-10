#include "rArB.h"

static int fastrand()
{
	g_seed = (214013 * g_seed + 2531011); 
	return (g_seed>>16) & 0x7FFF; 
}

static void
fillMatrix(int **arr, int n, int flag) 
{
	int count = 0;

	for(int i = 0; i < n; i++) 
        	for(int j = 0; j < n; j++) 
			if (flag)
        	    		arr[i][j] = fastrand();
			else
				arr[i][j] = 0;
}

static void 
MMultiply(int **X, int **Y, int **Z, int x_row, int x_col, 
int y_row, int y_col, int z_row, int z_col, int n)
{
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) 
            for (int j = 0; j < n; j++) 
                Z[z_row + i][z_col + j] += X[x_row + i][x_col + k] * Y[y_row + k][y_col + j];
}

static void
updateZ(int *updateZ, int myrank,
int **Z, int blocks , int block_size)
{
    int row = myrank / blocks;
    int col = myrank % blocks;

    for (int i = row * block_size; i < (row * block_size) + block_size; i++) {
        for (int j = col * block_size; j < (col * block_size) + block_size; j++) {
            Z[i][j] += *updateZ;
            updateZ++;
        }
    }
}

static int *
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

static void delete_2d_array(int **X){
	free(X);
}
static void create_2d_array(int ***array, int block_size) {

    int *p = (int *)malloc(pow(block_size, 2)*sizeof(int));
    if(!p) {
        std::cout << "Malloc Failed !!!\n";
        return;
    }

    (*array) = (int**)malloc(block_size *sizeof(int*));
    if (!(*array)) {
       free(p);
       std::cout << "Malloc Failed !!!\n";
       return;
    }

    for (int i=0; i<block_size; i++) {
       (*array)[i] = &(p[i*block_size]);
    }
    return;
}
static void 
MM_rArB(int myrank,  int n, int p)
{
	int blocks = sqrt(p);
	int block_size = (n / sqrt(p));
	int **sub_matrix_X = NULL;
	int **sub_matrix_Y = NULL;
	int **sub_matrix_Z = NULL;	
	cout<<"blocks"<<blocks<<endl;
	cout<<"block size"<<block_size<<endl;
	MPI_Status lu_status;
	MPI_Status merge_status[p];
	MPI_Status x_status[p];
	MPI_Status y_status[p];
	MPI_Status z_status[p];
	
	int x_tag  = 100;
	int y_tag  = 101;
	int z_tag  = 102;
	int rs_tag = 103;
	int merge_tag = 104;
	cout<<"initialized tags"<<endl;
	create_2d_array(&sub_matrix_X, block_size);
	create_2d_array(&sub_matrix_Y, block_size);
	create_2d_array(&sub_matrix_Z, block_size);
	cout <<"created successfully"<<endl;
	for(int i= 0;i<block_size;i++){
		for(int j=0;j<block_size;j++){
			sub_matrix_Z[i][j] = 0;
		}
	}
	cout<<"Z initialized"<<endl;
	/*sub_matrix_X = (int**) malloc(sizeof(int*) * block_size);
	for(int i = 0; i< block_size;i++)
		sub_matrix_X[i] = (int*)malloc(block_size * sizeof(int));
	sub_matrix_Y = (int**) malloc(sizeof(int*) * block_size);
        for(int i = 0; i< block_size;i++)
                sub_matrix_Y[i] = (int*)malloc(block_size * sizeof(int));
	sub_matrix_Z = (int**) malloc(sizeof(int*) * block_size);
        for(int i = 0; i< block_size;i++)
                sub_matrix_Z[i] = (int*)malloc(block_size * sizeof(int));

	cout <<"malloced sub matrix X"<<sub_matrix_X<<endl;
	//&sub_matrix_Y = (int***) malloc(sizeof(int*) * pow(block_size, 2));
	cout <<"malloced sub_matrixY"<<sub_matrix_Y<<endl;
	//&sub_matrix_Z = (int***) malloc(sizeof(int*) * pow(block_size, 2));*/
	cout <<"malloced submatrix Z"<<sub_matrix_Z<<endl;
	cout<<"my rank"<<myrank<<endl;	
	if (myrank == 0) {
		for (int i = 0; i < block_size; i++) {
	   		for (int j = 0; j < block_size; j++) {
           		    cout <<"X"<<X<<endl;
			    sub_matrix_X[i][j] = *(*X + (i * n) + j);
			    cout <<"sub_matX[i][j]"<<sub_matrix_X[i][j]<<endl;
	   		    sub_matrix_Y[i][j] = *(*Y + (i * n)+ j);
			    cout<<"sub matY[i][j]"<<sub_matrix_Y[i][j]<<endl;
	   		    sub_matrix_Z[i][j] = *(*Z + (i * block_size) + j);
			    cout<<"submatZ[i][j]"<<sub_matrix_Z[i][j]<<endl;
	  		 }
	 	}
		cout<<"done with first for loop"<<endl; 
		for (int i = 1; i < blocks * blocks; i++) {
			int *block_X = init_block(i, X, blocks, block_size);
			MPI_Send(&block_X, block_size * block_size, MPI_INT, i, x_tag, MPI_COMM_WORLD);
		
			int *block_Y = init_block(i, Y, blocks, block_size);
			MPI_Send(&block_Y, block_size * block_size, MPI_INT, i, y_tag, MPI_COMM_WORLD);
			
			int *block_Z = init_block(i, Z, blocks, block_size);
			MPI_Send(&block_Z, block_size * block_size, MPI_INT, i, z_tag, MPI_COMM_WORLD);

			free(block_X);
			free(block_Y);
			free(block_Z);
		}
	} else if (myrank != 0) {
		cout<<"before 1st MPIREc"<<endl;
		MPI_Recv(&(sub_matrix_X[0][0]), pow(block_size, 2), MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_status[myrank]);
		cout<<"1st MPIRecv"<<endl;
		MPI_Recv(&(sub_matrix_Y[0][0]), pow(block_size, 2), MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_status[myrank]);
		cout<<"2nd MPIRecv"<<endl;	
		MPI_Recv(&(sub_matrix_Z[0][0]), pow(block_size, 2), MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_status[myrank]);
		cout<<"3rd MPI REcv "<<endl;
	}
	
	if (myrank > -1) {
	    int xy_row = myrank / blocks;
	    int xy_col = myrank % blocks;
	    
	    int send_xrow  = xy_row;
	    int send_xcol  = (blocks + xy_col - xy_row) % blocks;
	    int send_xproc = (blocks * send_xrow) + send_xcol;
	    
	    int send_yrow  = (blocks + xy_col - xy_row) % blocks;
	    int send_ycol  = xy_col;
	    int send_yproc = (blocks * send_yrow) + send_ycol;

	    int rec_xrow  = xy_row;
	    int rec_xcol  = (blocks + xy_col + xy_row) % blocks;
	    int rec_xproc = (blocks * rec_xrow) + rec_xcol;

	    int rec_yrow  = (blocks + xy_col + xy_row) % blocks;
	    int rec_ycol  = xy_col;
	    int rec_yproc = (blocks * rec_yrow) + rec_ycol;

	    MPI_Send(&(sub_matrix_X[0][0]), pow(block_size, 2), MPI_INT, send_xproc, rs_tag, MPI_COMM_WORLD);
	    MPI_Recv(&(sub_matrix_X[0][0]), pow(block_size, 2), MPI_INT, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);

	    MPI_Send(&(sub_matrix_Y[0][0]), pow(block_size, 2), MPI_INT, send_xproc, rs_tag, MPI_COMM_WORLD);
	    MPI_Recv(&(sub_matrix_Y[0][0]), pow(block_size, 2), MPI_INT, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);

	   for (int l = 1; l < blocks; l++) {
		MMultiply(sub_matrix_X, sub_matrix_Y, sub_matrix_Z,
				0, 0, 0, 0, 0, 0, block_size);
		send_xrow = xy_row;
		send_xcol  = (blocks + xy_col - 1) % blocks;
		send_xproc = (blocks * send_xrow) + send_xcol;

		send_yrow  = (blocks + xy_col - 1) % blocks;
		send_ycol  = xy_col;
		send_yproc = (blocks * send_yrow) + send_ycol;

		rec_xrow  = xy_row;
		rec_xcol  = (blocks + xy_col + 1) % blocks;
		rec_xproc = (blocks * rec_xrow) + rec_xcol;

		rec_yrow  = (blocks + xy_col + 1) % blocks;
		rec_ycol  = xy_col;
		rec_yproc = (blocks * rec_yrow) + rec_ycol;

		MPI_Send(&(sub_matrix_X[0][0]), pow(block_size, 2), MPI_INT, send_xproc, rs_tag, MPI_COMM_WORLD);
		MPI_Recv(&(sub_matrix_X[0][0]), pow(block_size, 2), MPI_INT, rec_xproc,rs_tag, MPI_COMM_WORLD, &lu_status);

		MPI_Send(&(sub_matrix_Y[0][0]), pow(block_size, 2), MPI_INT, send_xproc, rs_tag, MPI_COMM_WORLD);
		MPI_Recv(&(sub_matrix_Y[0][0]), pow(block_size, 2), MPI_INT, rec_xproc, rs_tag, MPI_COMM_WORLD, &lu_status);
	   }

	   if (myrank != 0)
	       MPI_Send(&(sub_matrix_Z[0][0]), pow(block_size, 2), MPI_INT, 0, merge_tag, MPI_COMM_WORLD);
	   else if (myrank == 0) {
		updateZ(&(sub_matrix_Z[0][0]), 0, Z, n , p);
		int *update_Z = (int *) malloc(sizeof(int) * pow(block_size, 2));
	       
		for (int i = 1; i < pow(blocks, 2); i++) {
			MPI_Recv(update_Z, pow(block_size, 2) , MPI_INT, i, merge_tag, MPI_COMM_WORLD, &merge_status[i]);
			updateZ(update_Z, i, Z, n , p);
	       }
	       free(update_Z);
	   }

	   free(sub_matrix_X);
	   free(sub_matrix_Y);
	   free(sub_matrix_Z);

	}
}
static int input_valid(int argc){
	if(argc != 2)
		return 1;
	else 
		return 0;
}
int main(int argc, char *argv[])
{	cout<<"argc"<<argc<<endl;
	int ret = input_valid(argc);
	cout<<"ret is:"<<ret<<endl;
	if (ret) {
		std::cout << "Missing Input params - ibrun -n <processors> <binary> <matrix_size>\n";
		exit(1);
	}
	cout <<"outisede if"<<endl;
	int n = atoi(argv[1]);
	cout <<"n is"<<n <<endl;
	int processors;
	int myrank; 
	srand(time(NULL));
	cout <<"srand"<<endl;
	g_seed=rand();
	cout <<"gseed set"<<endl;
	std::chrono::system_clock::time_point start;
	std::chrono::system_clock::time_point finish;
	cout<<"Initialized variabled"<<endl;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &processors);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	cout<<"set up done"<<endl;
	if (myrank == 0) {
		X = (int**) malloc(sizeof(int*) * pow(n, 2));
		cout<<"mallced X"<<X<<endl;
		fillMatrix(X, n, 1);
		cout <<"filled Matrix X"<<endl;
	} else if (myrank == 1) {
		Y = (int**) malloc(sizeof(int*) * pow(n, 2));
		cout<<"malloced Y"<<Y<<endl;
		fillMatrix(Y, n, 1);
		cout<<"Filled Matrix Y"<<endl;
	} else if (myrank == 2) {
		Z = (int**) malloc(sizeof(int*) * pow(n, 2));
		cout<<"malloced Z"<<Z<<endl;
		fillMatrix(Z, n, 0);
		cout <<"filled matrix Z"<<endl;
	}

	if (myrank == 0) 
	    start = std::chrono::high_resolution_clock::now();
	cout<<"before calling rotatte function"<<endl;
	MM_rArB(n, processors, myrank);
	cout <<"after calling rotate function"<<endl;
	
	if (myrank == 0) {
		finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time: " << elapsed.count() << " s\n";
		delete_2d_array(X);
		delete_2d_array(Y);
		delete_2d_array(Z);
	}
	
	MPI_Finalize();

	return 0;
}
