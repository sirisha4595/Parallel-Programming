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

static int fastrand();

static void
new_matrix(int ***array, int n, int flag);

static void
free_matrix(int ***arr);

static void
init_matrix(int **arr, int n, int flag);

static void 
MMultiply(int **X, int **Y, int **Z,
int x_row, int x_col, int y_row, int y_col,
int z_row, int z_col, int n);

static void
updateZ(int *updateZ, int myrank,
int **Z, int blocks , int block_size);

static int *
clone(int rank, int  **arr,
int blocks, int block_size);

static void
MM_rArB(int myrank, int n, int p);

static void 
init_var(int argc, char **argv,
int *n, int &myrank, int &p);
