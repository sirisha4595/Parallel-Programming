#include <iostream>
#include <cmath>
using namespace std;
    
static void 
MM_rArB(int n, int p)
{
    int sqrt_p = sqrt(p);
    int block_size = (n / sqrt_p);

