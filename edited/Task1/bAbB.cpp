#include "bAbB.h"


static int fastrand() {
  g_seed = (214013 * g_seed + 2531011); 
  return (g_seed>>16) & 0x7FFF; 
}

static void create_2d_array(int ***array, int n, int m) {

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

static void delete_2d_array(int ***array) {
    free(&((*array)[0][0]));
    free(*array);
    return;
}

static void fillMatrix(int **arr, int n) {
    int count = 0;
    for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            arr[i][j] = /*++count;*/ fastrand();
        }
    }
}

static void fillValuesInMatrix(int **arr, int n, int val) {
    for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            arr[i][j] = val; //fastrand();
        }
    }
}


static void printMatrix(int **arr, int n) {
     std::cout << "\n";
     for(int i = 0; i<n; ++i) {
        for(int j = 0; j<n; ++j) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }   
}

static void init_block(int **arr, int m, int n) {
    for(int i = 0; i< m; ++i) {
        for(int j = 0; j< m; ++j) {
            arr[i][j] = 0;
        }
    }
}


static void MMultiply(int **X, int **Y, int **Z, int x_row, int x_col, 
                     int y_row, int y_col, int z_row, int z_col, int n) {
    for(int i = 0; i<n; ++i){
        for(int k = 0; k<n; ++k) {
            for(int j = 0; j<n; ++j) {
                Z[z_row + i][z_col + j] += X[x_row + i][x_col + k] * Y[y_row + k][y_col + j];
            }
        }
    }
}


static void updateZ(int *temp, int rank, int **Z, int n , int p) {
    int sqrt_p = std::sqrt(p);
    int proc_mat_size = n / sqrt_p;

    int row = rank / sqrt_p;
    int col = rank % sqrt_p;
    int *trav = temp;

    for (int i = row * proc_mat_size; i < row * proc_mat_size + proc_mat_size; ++i) {
        for (int j = col * proc_mat_size; j < col * proc_mat_size + proc_mat_size; ++j) {

            //std::cout << " i = " << i << " j = " << j << " val = " << *trav;
            Z[i][j] += *trav;
            trav++;
        }
        //std::cout << "\n";
    }
}

static int * get_buffer(int rank,int  **arr, int size, int procs) {
    int sqrt_p = std::sqrt(procs);
    int proc_mat_size = size / sqrt_p;

    int * res = (int *)malloc(proc_mat_size*proc_mat_size*sizeof(int));
    int *trav = res;
    int row = rank / sqrt_p;
    int col = rank % sqrt_p;

    for (int i = row * proc_mat_size; i < row * proc_mat_size + proc_mat_size; ++i) {
        for (int j = col * proc_mat_size; j < col * proc_mat_size + proc_mat_size; ++j) {
            *trav = arr[i][j];
            trav++;
        }
    }

    return res;
}

static void MM_bAbB(int n, int p, int rank) {
    int sqrt_p = std::sqrt(p);
    int proc_mat_size = n / sqrt_p;
    int mat_size[2] = {n , n};
    int sub_mat_size[2] = {proc_mat_size, proc_mat_size};
    int mat_start[2] = {0, 0};

    int **sub_matrix_X;
    create_2d_array(&sub_matrix_X, proc_mat_size, proc_mat_size);

    int **sub_matrix_Y;
    create_2d_array(&sub_matrix_Y, proc_mat_size, proc_mat_size);

    int **sub_matrix_Z;
    create_2d_array(&sub_matrix_Z, proc_mat_size, proc_mat_size);
    init_block(sub_matrix_Z, proc_mat_size, proc_mat_size);

    MPI_Status status;
    MPI_Datatype mat_type, sub_mat_type;
    MPI_Type_create_subarray(2, mat_size, sub_mat_size, mat_start, MPI_ORDER_C, MPI_INT, &mat_type);
    MPI_Type_create_resized(mat_type, 0, proc_mat_size * sizeof(int), &sub_mat_type);
    MPI_Type_commit(&sub_mat_type);

    int color = rank % sqrt_p;
    MPI_Comm col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &col_comm);

    int row_color = sqrt_p + (rank / sqrt_p);
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, row_color, rank, &row_comm);

    int send_mat_count[p]; 
    int send_mat_start[p];


    MPI_Status x_status[n];
    MPI_Status y_status[n];
    MPI_Status z_status[n];

    int x_tag = 40;
    int y_tag = 41;
    int z_tag = 42;



    if(rank == 0) {

    for (int i = 0; i < proc_mat_size; ++i) {
       for (int j = 0; j < proc_mat_size; ++j) {
               sub_matrix_X[i][j] = *(*X + i * n + j);
           sub_matrix_Y[i][j] = *(*Y + i * n + j);
           sub_matrix_Z[i][j] = *(*Z + i * proc_mat_size + j);
       }
    }
        for (int i = 1; i < sqrt_p * sqrt_p; ++i) {
            int *sub_X = get_buffer(i, X, n, p);
            MPI_Send(sub_X, proc_mat_size * proc_mat_size, MPI_INT, i, x_tag, MPI_COMM_WORLD);

            int *sub_Y = get_buffer(i, Y, n, p);
            MPI_Send(sub_Y, proc_mat_size * proc_mat_size, MPI_INT, i, y_tag, MPI_COMM_WORLD);

            int *sub_Z = get_buffer(i, Z, n, p);
            MPI_Send(sub_Z, proc_mat_size * proc_mat_size, MPI_INT, i, z_tag, MPI_COMM_WORLD);

            free(sub_X);
            free(sub_Y);
            free(sub_Z);
        }

    }

    if (rank != 0) {

        MPI_Recv(&(sub_matrix_X[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, x_tag, MPI_COMM_WORLD, &x_status[rank]);
        MPI_Recv(&(sub_matrix_Y[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, y_tag, MPI_COMM_WORLD, &y_status[rank]);
        MPI_Recv(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, z_tag, MPI_COMM_WORLD, &z_status[rank]);
    }



    int tag = 99;
    int p_src_row = rank / sqrt_p; 
    int p_src_col = rank % sqrt_p;

    int x_send_to_p_dst_row = p_src_row;     
    int x_recv_from_p_src_row = p_src_row;

    int **a_bcast_ptr;
    create_2d_array(&a_bcast_ptr, proc_mat_size, proc_mat_size);
    init_block(a_bcast_ptr, proc_mat_size, proc_mat_size);
    int *a_g_bcast_ptr = &(a_bcast_ptr[0][0]);

    int **b_bcast_ptr;
    create_2d_array(&b_bcast_ptr, proc_mat_size, proc_mat_size);
    init_block(b_bcast_ptr, proc_mat_size, proc_mat_size);
    int *b_g_bcast_ptr = &(b_bcast_ptr[0][0]);
    
    for(int l = 1; l <= sqrt_p; l++) {
        int k = (l - 1);
        if (k == p_src_row) {
           b_bcast_ptr = sub_matrix_Y;
        }

        if (k == p_src_col) {
        a_bcast_ptr = sub_matrix_X;
        }

        MPI_Bcast(&(b_bcast_ptr[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 
                  k, col_comm);

        MPI_Bcast(&(a_bcast_ptr[0][0]), proc_mat_size * proc_mat_size, MPI_INT,
                  k, row_comm);

        MMultiply(a_bcast_ptr, b_bcast_ptr, sub_matrix_Z,
                        0, 0, 0, 0, 0, 0, proc_mat_size);

    }


    int rcv_tag;
    MPI_Status rcv_status[n];

    if (rank != 0) {
        MPI_Send(&(sub_matrix_Z[0][0]), proc_mat_size * proc_mat_size, MPI_INT, 0, rcv_tag, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        updateZ(&(sub_matrix_Z[0][0]), 0, Z, n , p);
        for (int i = 1; i < sqrt_p * sqrt_p; ++i) {
            int *temp = (int *)malloc(proc_mat_size*proc_mat_size*sizeof(int));
            MPI_Recv(temp, proc_mat_size * proc_mat_size, MPI_INT, i, rcv_tag, MPI_COMM_WORLD, &rcv_status[i]);
            updateZ(temp, i, Z, n , p);
            free(temp);
        }
    }


    return;
}

int main(int argc, char *argv[]) {
    if(argc < 2)  {
        std::cout << "Missing Input params - ibrun -n <procs> <binary> <matrix_size>\n";
        return 1;
    }
    int n = atoi(argv[1]);

    srand(time(NULL));
    g_seed=rand();

    int myrank, v = 121;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    
    if(myrank == 0) {
        create_2d_array(&X, n, n);
        fillMatrix(X, n);
        //printMatrix(X, n);

        create_2d_array(&Y, n, n);
        fillMatrix(Y, n);

        create_2d_array(&Z, n, n);
        init_block(Z, n, n);
    }


    using namespace std::chrono;    
    high_resolution_clock::time_point start_time = high_resolution_clock::now();

    MM_bAbB(n, procs, myrank);

    high_resolution_clock::time_point end_time = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end_time - start_time);

    if(myrank == 0) {
        std::cout << "Exectution Time: " << time_span.count() << " seconds.";
        std::cout << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 1;
}
