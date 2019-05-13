#include <iostream>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <chrono>

using namespace std;

// Input Matrices
int **X; 
int **Y;
int **Z;

int g_seed;

int procs = 1;

static int
fastrand();

static void
create_2d_array(int ***array, int n, int m);

static void 
delete_2d_array(int ***array);

static void 
fillMatrix(int **arr, int n);

static void 
fillValuesInMatrix(int **arr, int n, int val);

static void 
printMatrix(int **arr, int n);

static void 
init_block(int **arr, int m, int n);

static void 
MMultiply(int **X, int **Y, int **Z, int x_row, int x_col, 
                     int y_row, int y_col, int z_row, int z_col, int n);

static void 
parallel_MMultiply(int **X, int **Y, int **Z, int x_row, int x_col,
				 int y_row, int y_col, int z_row, int z_col, int n, int base);

static void 
updateZ(int *temp, int rank, int **Z, int n , int p);

static int * 
get_buffer(int rank,int  **arr, int size, int procs);

static void 
MM_bAbB(int n, int p, int myrank);