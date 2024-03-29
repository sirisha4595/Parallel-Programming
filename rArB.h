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

static int
fastrand();

static int
fillMatrix(int **arr, int n, int flag);

static void 
MMultiply(int **X, int **Y, int **Z, int x_row, int x_col,
int y_row, int y_col, int z_row, int z_col, int n);

static void
updateZ(int *updateZ, int myrank,
int **Z, int blocks , int block_size);

static int *
init_block (int myrank, int **matrix,
int blocks, int block_size);

static void
MM_rArB(int myrank, int n, int p);

static void delete_2d_array(int **X);

static int input_valid(int argc);
static void create_2d_array(int ***array, int block_size); 
