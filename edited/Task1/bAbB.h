#include <iostream>
#include <cmath>
#include <mpi.h>
#include <algorithm>
#include <chrono>
using namespace std;

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

static int * 
get_buffer(int myrank,int  **arr, int size, int procs);

static void 
updateZ(int *temp, int myrank, int **Z, int n , int p);

static void 
MM_bAbB(int n, int p, int myrank);