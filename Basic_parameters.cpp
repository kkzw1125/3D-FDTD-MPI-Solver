#include "Basic_parameters.h"


// 在这里定义所有 extern 变量
int MAX_X = 100;
int MAX_Y = 80;
int MAX_Z = 50;

double dx = 1e-2;
double dy = 1e-2;
double dz = 1e-2;
double dt = 0.9 / sqrt(1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz)) / c0;
const int timesteps = 300;

double CA = 1.0;
double CP = 1.0;

// 场分量和系数数组定义
double*** Ex, *** Ey, *** Ez, *** Hx, *** Hy, *** Hz;
double*** CBx, *** CBy, *** CBz;
double*** CQx, *** CQy, *** CQz;
/*
double* ez_inc, * hx_inc;
// 网格尺寸
int i0 = 11, j_0 = 11, k0 = 11;
int ia = MAX_X - 11, jb = MAX_Y - 11, kc =  MAX_Z - 11;
*/
void Initialization_grid() {
    Ex = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z + 1, 0);
    Ey = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z + 1, 0);
    Ez = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z + 1, 0);
    Hx = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z + 1, 0);
    Hy = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z + 1, 0);
    Hz = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z + 1, 0);

    CBx = create3DArray(MAX_X, MAX_Y + 1, MAX_Z + 1, dt / eps0);
    CBy = create3DArray(MAX_X + 1, MAX_Y, MAX_Z + 1, dt / eps0);
    CBz = create3DArray(MAX_X + 1, MAX_Y + 1, MAX_Z, dt / eps0);
    CQx = create3DArray(MAX_X + 1, MAX_Y, MAX_Z, dt / mu0);
    CQy = create3DArray(MAX_X, MAX_Y + 1, MAX_Z, dt / mu0);
    CQz = create3DArray(MAX_X, MAX_Y, MAX_Z + 1, dt / mu0);

    /*// 入射波
    ez_inc = create1DArray( MAX_Y + 1, 0);
    hx_inc = create1DArray( MAX_Y , 0);
    */
}