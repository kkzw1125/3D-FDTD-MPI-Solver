#pragma once
#ifndef BASIC_PARAMETERS_H
#define BASIC_PARAMETERS_H
#include"math.h"
#include"Basic_functions.h"

// 常数
const double PI = 3.141592654;
const double eps0 = 8.85e-12;
const double epsr = 1;
const double pi = 3.1415926;
const double mu0 = PI * 4e-7;
const double c0 = 1 / sqrt(eps0 * mu0);

// 网格数
extern int MAX_X;
extern int MAX_Y;
extern int MAX_Z;
/*
extern int i0, j_0 , k0 ;
extern int ia , jb , kc ;
*/
// 时间步长和空间步长
extern double dx;
extern double dy;
extern double dz;
extern double dt;
extern const int timesteps;

extern double CA;
extern double CP;

// 场分量和系数数组（改为 extern 声明，不在这里定义）
extern double*** Ex, *** Ey, *** Ez, *** Hx, *** Hy, *** Hz;
extern double*** CBx, *** CBy, *** CBz;
extern double*** CQx, *** CQy, *** CQz;
/*
extern double* ez_inc, * hx_inc;
*/
void Initialization_grid();

#endif