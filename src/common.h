#ifndef __COMMON_H__
#define  __COMMON_H__

#include <iostream>
#include <algorithm>
typedef double Float;

const Float alpha = 9.8;
const Float beta =  15.0;
const Float T_ambient = 100;
const Float epsilon = 0.25;
const int CHUNKSIZE = 1;
const Float eps = 1e-8;
#define ijk2idx(i,j,k) ((i) + grid.grid_num_x * (j) + grid.grid_num_x * grid.grid_num_y * (k))

using std::cout;
using std::endl;


#endif