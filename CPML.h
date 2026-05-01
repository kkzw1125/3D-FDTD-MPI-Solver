#pragma once
#ifndef CPML_H
#define CPML_H

#include "Basic_functions.h"
#include "Basic_parameters.h"
#include "math.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
using namespace std;

extern int ihT;
extern int jhT;
extern int khT;
extern int ieT;
extern int jeT;
extern int keT;

extern int iBC, jBC, kBC;
// CPML参数设置
extern int deep;
extern int nxPML1, nxPML2;
extern int nyPML1, nyPML2;
extern int nzPML1, nzPML2;

extern int m;
extern double sig_x_max, sig_y_max, sig_z_max;
extern double kap_x_max, kap_y_max, kap_z_max;
extern double alp_x_max, alp_y_max, alp_z_max;

// CPML辅助变量（改为指针形式）
extern double*** psi_Exy, *** psi_Exz, *** psi_Eyx, *** psi_Eyz, *** psi_Ezx, *** psi_Ezy;
extern double*** psi_Hxy, *** psi_Hxz, *** psi_Hyx, *** psi_Hyz, *** psi_Hzx, *** psi_Hzy;
// CPML系数
extern double* be_x1, * ce_x1, * sige_xPML1, * kape_xPML1, * alpe_xPML1;
extern double* be_x2, * ce_x2, * sige_xPML2, * kape_xPML2, * alpe_xPML2;

extern double* be_y1, * ce_y1, * sige_yPML1, * kape_yPML1, * alpe_yPML1;
extern double* be_y2, * ce_y2, * sige_yPML2, * kape_yPML2, * alpe_yPML2;

extern double* be_z1, * ce_z1, * sige_zPML1, * kape_zPML1, * alpe_zPML1;
extern double* be_z2, * ce_z2, * sige_zPML2, * kape_zPML2, * alpe_zPML2;

extern double* bh_x1, * ch_x1, * sigh_xPML1, * kaph_xPML1, * alph_xPML1;
extern double* bh_x2, * ch_x2, * sigh_xPML2, * kaph_xPML2, * alph_xPML2;

extern double* bh_y1, * ch_y1, * sigh_yPML1, * kaph_yPML1, * alph_yPML1;
extern double* bh_y2, * ch_y2, * sigh_yPML2, * kaph_yPML2, * alph_yPML2;

extern double* bh_z1, * ch_z1, * sigh_zPML1, * kaph_zPML1, * alph_zPML1;
extern double* bh_z2, * ch_z2, * sigh_zPML2, * kaph_zPML2, * alph_zPML2;

// 创建完整的系数数组
extern double* be_x_tmp, * ce_x_tmp;
extern double* be_y_tmp, * ce_y_tmp;
extern double* be_z_tmp, * ce_z_tmp;

extern double* bh_x_tmp, * ch_x_tmp;
extern double* bh_y_tmp, * ch_y_tmp;
extern double* bh_z_tmp, * ch_z_tmp;

// 函数声明
void init_cpml_variables();
void calculate_cpml();
void free_cpml_variables();

#endif