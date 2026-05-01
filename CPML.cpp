#include "CPML.h"

// 定义全局变量
int ihT;
int jhT;
int khT;
int ieT;
int jeT;
int keT;

int iBC = 9, jBC = 9, kBC = 9;
// CPML参数设置
int deep = iBC;
int nxPML1 = deep, nxPML2 = deep;
int nyPML1 = deep, nyPML2 = deep;
int nzPML1 = deep, nzPML2 = deep;

int m = 4;
double sig_x_max, sig_y_max, sig_z_max;
double kap_x_max = 1.0, kap_y_max = 1.0, kap_z_max = 1.0;
double alp_x_max = 0.0, alp_y_max = 0.0, alp_z_max = 0.0;

// CPML辅助变量
double*** psi_Exy = NULL, *** psi_Exz = NULL, *** psi_Eyx = NULL, *** psi_Eyz = NULL, *** psi_Ezx = NULL, *** psi_Ezy = NULL;
double*** psi_Hxy = NULL, *** psi_Hxz = NULL, *** psi_Hyx = NULL, *** psi_Hyz = NULL, *** psi_Hzx = NULL, *** psi_Hzy = NULL;

// CPML系数
double* be_x1 = NULL, * ce_x1 = NULL, * sige_xPML1 = NULL, * kape_xPML1 = NULL, * alpe_xPML1 = NULL;
double* be_x2 = NULL, * ce_x2 = NULL, * sige_xPML2 = NULL, * kape_xPML2 = NULL, * alpe_xPML2 = NULL;

double* be_y1 = NULL, * ce_y1 = NULL, * sige_yPML1 = NULL, * kape_yPML1 = NULL, * alpe_yPML1 = NULL;
double* be_y2 = NULL, * ce_y2 = NULL, * sige_yPML2 = NULL, * kape_yPML2 = NULL, * alpe_yPML2 = NULL;

double* be_z1 = NULL, * ce_z1 = NULL, * sige_zPML1 = NULL, * kape_zPML1 = NULL, * alpe_zPML1 = NULL;
double* be_z2 = NULL, * ce_z2 = NULL, * sige_zPML2 = NULL, * kape_zPML2 = NULL, * alpe_zPML2 = NULL;

double* bh_x1 = NULL, * ch_x1 = NULL, * sigh_xPML1 = NULL, * kaph_xPML1 = NULL, * alph_xPML1 = NULL;
double* bh_x2 = NULL, * ch_x2 = NULL, * sigh_xPML2 = NULL, * kaph_xPML2 = NULL, * alph_xPML2 = NULL;

double* bh_y1 = NULL, * ch_y1 = NULL, * sigh_yPML1 = NULL, * kaph_yPML1 = NULL, * alph_yPML1 = NULL;
double* bh_y2 = NULL, * ch_y2 = NULL, * sigh_yPML2 = NULL, * kaph_yPML2 = NULL, * alph_yPML2 = NULL;

double* bh_z1 = NULL, * ch_z1 = NULL, * sigh_zPML1 = NULL, * kaph_zPML1 = NULL, * alph_zPML1 = NULL;
double* bh_z2 = NULL, * ch_z2 = NULL, * sigh_zPML2 = NULL, * kaph_zPML2 = NULL, * alph_zPML2 = NULL;

// 创建完整的系数数组
double* be_x_tmp = NULL, * ce_x_tmp = NULL;
double* be_y_tmp = NULL, * ce_y_tmp = NULL;
double* be_z_tmp = NULL, * ce_z_tmp = NULL;

double* bh_x_tmp = NULL, * ch_x_tmp = NULL;
double* bh_y_tmp = NULL, * ch_y_tmp = NULL;
double* bh_z_tmp = NULL, * ch_z_tmp = NULL;

// 初始化CPML变量
void init_cpml_variables() {
    // 计算网格大小
    ieT = MAX_X;
    jeT = MAX_Y;
    keT = MAX_Z;
    ihT = ieT+1;
    jhT = jeT+1;
    khT = keT+1;

    // 计算最大sigma值
    sig_x_max = (m + 1) / (sqrt(epsr) * 150 * pi * dx);
    sig_y_max = (m + 1) / (sqrt(epsr) * 150 * pi * dy);
    sig_z_max = (m + 1) / (sqrt(epsr) * 150 * pi * dz);

    // 分配CPML辅助变量内存
   // 分配CPML辅助变量内存（使用 create3DArray）
    psi_Exy = create3DArray(ieT, jeT, keT, 0.0);
    psi_Exz = create3DArray(ieT, jeT, keT, 0.0);
    psi_Eyx = create3DArray(ieT, jeT, keT, 0.0);
    psi_Eyz = create3DArray(ieT, jeT, keT, 0.0);
    psi_Ezx = create3DArray(ieT, jeT, keT, 0.0);
    psi_Ezy = create3DArray(ieT, jeT, keT, 0.0);

    psi_Hxy = create3DArray(ihT, jhT, khT, 0.0);
    psi_Hxz = create3DArray(ihT, jhT, khT, 0.0);
    psi_Hyx = create3DArray(ihT, jhT, khT, 0.0);
    psi_Hyz = create3DArray(ihT, jhT, khT, 0.0);
    psi_Hzx = create3DArray(ihT, jhT, khT, 0.0);
    psi_Hzy = create3DArray(ihT, jhT, khT, 0.0);

    // 分配系数数组
    be_x1 = create1DArray(nxPML1, 0.0);
    ce_x1 = create1DArray(nxPML1, 0.0);
    sige_xPML1 = create1DArray(nxPML1, 0.0);
    kape_xPML1 = create1DArray(nxPML1, 0.0);
    alpe_xPML1 = create1DArray(nxPML1, 0.0);

    be_x2 = create1DArray(nxPML2, 0.0);
    ce_x2 = create1DArray(nxPML2, 0.0);
    sige_xPML2 = create1DArray(nxPML2, 0.0);
    kape_xPML2 = create1DArray(nxPML2, 0.0);
    alpe_xPML2 = create1DArray(nxPML2, 0.0);

    be_y1 = create1DArray(nyPML1, 0.0);
    ce_y1 = create1DArray(nyPML1, 0.0);
    sige_yPML1 = create1DArray(nyPML1, 0.0);
    kape_yPML1 = create1DArray(nyPML1, 0.0);
    alpe_yPML1 = create1DArray(nyPML1, 0.0);

    be_y2 = create1DArray(nyPML2, 0.0);
    ce_y2 = create1DArray(nyPML2, 0.0);
    sige_yPML2 = create1DArray(nyPML2, 0.0);
    kape_yPML2 = create1DArray(nyPML2, 0.0);
    alpe_yPML2 = create1DArray(nyPML2, 0.0);

    be_z1 = create1DArray(nzPML1, 0.0);
    ce_z1 = create1DArray(nzPML1, 0.0);
    sige_zPML1 = create1DArray(nzPML1, 0.0);
    kape_zPML1 = create1DArray(nzPML1, 0.0);
    alpe_zPML1 = create1DArray(nzPML1, 0.0);

    be_z2 = create1DArray(nzPML2, 0.0);
    ce_z2 = create1DArray(nzPML2, 0.0);
    sige_zPML2 = create1DArray(nzPML2, 0.0);
    kape_zPML2 = create1DArray(nzPML2, 0.0);
    alpe_zPML2 = create1DArray(nzPML2, 0.0);

    bh_x1 = create1DArray(nxPML1 - 1, 0.0);
    ch_x1 = create1DArray(nxPML1 - 1, 0.0);
    sigh_xPML1 = create1DArray(nxPML1 - 1, 0.0);
    kaph_xPML1 = create1DArray(nxPML1 - 1, 0.0);
    alph_xPML1 = create1DArray(nxPML1 - 1, 0.0);

    bh_x2 = create1DArray(nxPML2 - 1, 0.0);
    ch_x2 = create1DArray(nxPML2 - 1, 0.0);
    sigh_xPML2 = create1DArray(nxPML2 - 1, 0.0);
    kaph_xPML2 = create1DArray(nxPML2 - 1, 0.0);
    alph_xPML2 = create1DArray(nxPML2 - 1, 0.0);

    bh_y1 = create1DArray(nyPML1 - 1, 0.0);
    ch_y1 = create1DArray(nyPML1 - 1, 0.0);
    sigh_yPML1 = create1DArray(nyPML1 - 1, 0.0);
    kaph_yPML1 = create1DArray(nyPML1 - 1, 0.0);
    alph_yPML1 = create1DArray(nyPML1 - 1, 0.0);

    bh_y2 = create1DArray(nyPML2 - 1, 0.0);
    ch_y2 = create1DArray(nyPML2 - 1, 0.0);
    sigh_yPML2 = create1DArray(nyPML2 - 1, 0.0);
    kaph_yPML2 = create1DArray(nyPML2 - 1, 0.0);
    alph_yPML2 = create1DArray(nyPML2 - 1, 0.0);

    bh_z1 = create1DArray(nzPML1 - 1, 0.0);
    ch_z1 = create1DArray(nzPML1 - 1, 0.0);
    sigh_zPML1 = create1DArray(nzPML1 - 1, 0.0);
    kaph_zPML1 = create1DArray(nzPML1 - 1, 0.0);
    alph_zPML1 = create1DArray(nzPML1 - 1, 0.0);

    bh_z2 = create1DArray(nzPML2 - 1, 0.0);
    ch_z2 = create1DArray(nzPML2 - 1, 0.0);
    sigh_zPML2 = create1DArray(nzPML2 - 1, 0.0);
    kaph_zPML2 = create1DArray(nzPML2 - 1, 0.0);
    alph_zPML2 = create1DArray(nzPML2 - 1, 0.0);

    // 创建完整的系数数组
    be_x_tmp = create1DArray(ihT, 0.0);
    ce_x_tmp = create1DArray(ihT, 0.0);
    be_y_tmp = create1DArray(jhT, 0.0);
    ce_y_tmp = create1DArray(jhT, 0.0);
    be_z_tmp = create1DArray(khT, 0.0);
    ce_z_tmp = create1DArray(khT, 0.0);

    bh_x_tmp = create1DArray(ieT, 0.0);
    ch_x_tmp = create1DArray(ieT, 0.0);
    bh_y_tmp = create1DArray(jeT, 0.0);
    ch_y_tmp = create1DArray(jeT, 0.0);
    bh_z_tmp = create1DArray(keT, 0.0);
    ch_z_tmp = create1DArray(keT, 0.0);
}

// 计算CPML系数
void calculate_cpml() {
    // 计算电场CPML系数
    for (int i = 0; i < nxPML1; i++) {
        sige_xPML1[i] = sig_x_max * pow((nxPML1 - 1 - i) / (nxPML1 - 1.0), m);
        kape_xPML1[i] = 1 + (kap_x_max - 1) * pow((nxPML1 - 1 - i) / (nxPML1 - 1.0), m);
        alpe_xPML1[i] = alp_x_max * ((nxPML1 - 1 - i) / (nxPML1 - 1.0));
        be_x1[i] = exp(-(sige_xPML1[i] / (eps0 * kape_xPML1[i]) + alpe_xPML1[i] / eps0) * dt);

        if (sige_xPML1[i] == 0 && alpe_xPML1[i] == 0 && i == nxPML1 - 1) {
            ce_x1[i] = 0;
        }
        else {
            ce_x1[i] = sige_xPML1[i] / (sige_xPML1[i] * kape_xPML1[i] + kape_xPML1[i] * kape_xPML1[i] * alpe_xPML1[i]) * (be_x1[i] - 1);
        }
    }

    for (int i = 0; i < nxPML2; i++) {
        sige_xPML2[i] = sig_x_max * pow((nxPML2 - 1 - i) / (nxPML2 - 1.0), m);
        kape_xPML2[i] = 1 + (kap_x_max - 1) * pow((nxPML2 - 1 - i) / (nxPML2 - 1.0), m);
        alpe_xPML2[i] = alp_x_max * ((nxPML2 - 1 - i) / (nxPML2 - 1.0));
        be_x2[i] = exp(-(sige_xPML2[i] / (eps0 * kape_xPML2[i]) + alpe_xPML2[i] / eps0) * dt);

        if (sige_xPML2[i] == 0 && alpe_xPML2[i] == 0 && i == nxPML2 - 1) {
            ce_x2[i] = 0;
        }
        else {
            ce_x2[i] = sige_xPML2[i] / (sige_xPML2[i] * kape_xPML2[i] + kape_xPML2[i] * kape_xPML2[i] * alpe_xPML2[i]) * (be_x2[i] - 1);
        }
    }

    for (int i = 0; i < nyPML1; i++) {
        sige_yPML1[i] = sig_y_max * pow((nyPML1 - 1 - i) / (nyPML1 - 1.0), m);
        kape_yPML1[i] = 1 + (kap_y_max - 1) * pow((nyPML1 - 1 - i) / (nyPML1 - 1.0), m);
        alpe_yPML1[i] = alp_y_max * ((nyPML1 - 1 - i) / (nyPML1 - 1.0));
        be_y1[i] = exp(-(sige_yPML1[i] / (eps0 * kape_yPML1[i]) + alpe_yPML1[i] / eps0) * dt);

        if (sige_yPML1[i] == 0 && alpe_yPML1[i] == 0 && i == nyPML1 - 1) {
            ce_y1[i] = 0;
        }
        else {
            ce_y1[i] = sige_yPML1[i] / (sige_yPML1[i] * kape_yPML1[i] + kape_yPML1[i] * kape_yPML1[i] * alpe_yPML1[i]) * (be_y1[i] - 1);
        }
    }

    for (int i = 0; i < nyPML2; i++) {
        sige_yPML2[i] = sig_y_max * pow((nyPML2 - 1 - i) / (nyPML2 - 1.0), m);
        kape_yPML2[i] = 1 + (kap_y_max - 1) * pow((nyPML2 - 1 - i) / (nyPML2 - 1.0), m);
        alpe_yPML2[i] = alp_y_max * ((nyPML2 - 1 - i) / (nyPML2 - 1.0));
        be_y2[i] = exp(-(sige_yPML2[i] / (eps0 * kape_yPML2[i]) + alpe_yPML2[i] / eps0) * dt);

        if (sige_yPML2[i] == 0 && alpe_yPML2[i] == 0 && i == nyPML2 - 1) {
            ce_y2[i] = 0;
        }
        else {
            ce_y2[i] = sige_yPML2[i] / (sige_yPML2[i] * kape_yPML2[i] + kape_yPML2[i] * kape_yPML2[i] * alpe_yPML2[i]) * (be_y2[i] - 1);
        }
    }

    for (int i = 0; i < nzPML1; i++) {
        sige_zPML1[i] = sig_z_max * pow((nzPML1 - 1 - i) / (nzPML1 - 1.0), m);
        kape_zPML1[i] = 1 + (kap_z_max - 1) * pow((nzPML1 - 1 - i) / (nzPML1 - 1.0), m);
        alpe_zPML1[i] = alp_z_max * ((nzPML1 - 1 - i) / (nzPML1 - 1.0));
        be_z1[i] = exp(-(sige_zPML1[i] / (eps0 * kape_zPML1[i]) + alpe_zPML1[i] / eps0) * dt);

        if (sige_zPML1[i] == 0 && alpe_zPML1[i] == 0 && i == nzPML1 - 1) {
            ce_z1[i] = 0;
        }
        else {
            ce_z1[i] = sige_zPML1[i] / (sige_zPML1[i] * kape_zPML1[i] + kape_zPML1[i] * kape_zPML1[i] * alpe_zPML1[i]) * (be_z1[i] - 1);
        }
    }

    for (int i = 0; i < nzPML2; i++) {
        sige_zPML2[i] = sig_z_max * pow((nzPML2 - 1 - i) / (nzPML2 - 1.0), m);
        kape_zPML2[i] = 1 + (kap_z_max - 1) * pow((nzPML2 - 1 - i) / (nzPML2 - 1.0), m);
        alpe_zPML2[i] = alp_z_max * ((nzPML2 - 1 - i) / (nzPML2 - 1.0));
        be_z2[i] = exp(-(sige_zPML2[i] / (eps0 * kape_zPML2[i]) + alpe_zPML2[i] / eps0) * dt);

        if (sige_zPML2[i] == 0 && alpe_zPML2[i] == 0 && i == nzPML2 - 1) {
            ce_z2[i] = 0;
        }
        else {
            ce_z2[i] = sige_zPML2[i] / (sige_zPML2[i] * kape_zPML2[i] + kape_zPML2[i] * kape_zPML2[i] * alpe_zPML2[i]) * (be_z2[i] - 1);
        }
    }

    // 计算磁场CPML系数
    for (int i = 0; i < nxPML1 - 1; i++) {
        sigh_xPML1[i] = sig_x_max * pow((nxPML1 - 1.5 - i) / (nxPML1 - 1.0), m);
        kaph_xPML1[i] = 1 + (kap_x_max - 1) * pow((nxPML1 - 1.5 - i) / (nxPML1 - 1.0), m);
        alph_xPML1[i] = alp_x_max * ((nxPML1 - 1.5 - i) / (nxPML1 - 1.0));
        bh_x1[i] = exp(-(sigh_xPML1[i] / (eps0 * kaph_xPML1[i]) + alph_xPML1[i] / eps0) * dt);
        ch_x1[i] = sigh_xPML1[i] / (sigh_xPML1[i] * kaph_xPML1[i] + kaph_xPML1[i] * kaph_xPML1[i] * alph_xPML1[i]) * (bh_x1[i] - 1);
    }

    for (int i = 0; i < nxPML2 - 1; i++) {
        sigh_xPML2[i] = sig_x_max * pow((nxPML2 - 1.5 - i) / (nxPML2 - 1.0), m);
        kaph_xPML2[i] = 1 + (kap_x_max - 1) * pow((nxPML2 - 1.5 - i) / (nxPML2 - 1.0), m);
        alph_xPML2[i] = alp_x_max * ((nxPML2 - 1.5 - i) / (nxPML2 - 1.0));
        bh_x2[i] = exp(-(sigh_xPML2[i] / (eps0 * kaph_xPML2[i]) + alph_xPML2[i] / eps0) * dt);
        ch_x2[i] = sigh_xPML2[i] / (sigh_xPML2[i] * kaph_xPML2[i] + kaph_xPML2[i] * kaph_xPML2[i] * alph_xPML2[i]) * (bh_x2[i] - 1);
    }

    for (int i = 0; i < nyPML1 - 1; i++) {
        sigh_yPML1[i] = sig_y_max * pow((nyPML1 - 1.5 - i) / (nyPML1 - 1.0), m);
        kaph_yPML1[i] = 1 + (kap_y_max - 1) * pow((nyPML1 - 1.5 - i) / (nyPML1 - 1.0), m);
        alph_yPML1[i] = alp_y_max * ((nyPML1 - 1.5 - i) / (nyPML1 - 1.0));
        bh_y1[i] = exp(-(sigh_yPML1[i] / (eps0 * kaph_yPML1[i]) + alph_yPML1[i] / eps0) * dt);
        ch_y1[i] = sigh_yPML1[i] / (sigh_yPML1[i] * kaph_yPML1[i] + kaph_yPML1[i] * kaph_yPML1[i] * alph_yPML1[i]) * (bh_y1[i] - 1);
    }

    for (int i = 0; i < nyPML2 - 1; i++) {
        sigh_yPML2[i] = sig_y_max * pow((nyPML2 - 1.5 - i) / (nyPML2 - 1.0), m);
        kaph_yPML2[i] = 1 + (kap_y_max - 1) * pow((nyPML2 - 1.5 - i) / (nyPML2 - 1.0), m);
        alph_yPML2[i] = alp_y_max * ((nyPML2 - 1.5 - i) / (nyPML2 - 1.0));
        bh_y2[i] = exp(-(sigh_yPML2[i] / (eps0 * kaph_yPML2[i]) + alph_yPML2[i] / eps0) * dt);
        ch_y2[i] = sigh_yPML2[i] / (sigh_yPML2[i] * kaph_yPML2[i] + kaph_yPML2[i] * kaph_yPML2[i] * alph_yPML2[i]) * (bh_y2[i] - 1);
    }

    for (int i = 0; i < nzPML1 - 1; i++) {
        sigh_zPML1[i] = sig_z_max * pow((nzPML1 - 1.5 - i) / (nzPML1 - 1.0), m);
        kaph_zPML1[i] = 1 + (kap_z_max - 1) * pow((nzPML1 - 1.5 - i) / (nzPML1 - 1.0), m);
        alph_zPML1[i] = alp_z_max * ((nzPML1 - 1.5 - i) / (nzPML1 - 1.0));
        bh_z1[i] = exp(-(sigh_zPML1[i] / (eps0 * kaph_zPML1[i]) + alph_zPML1[i] / eps0) * dt);
        ch_z1[i] = sigh_zPML1[i] / (sigh_zPML1[i] * kaph_zPML1[i] + kaph_zPML1[i] * kaph_zPML1[i] * alph_zPML1[i]) * (bh_z1[i] - 1);
    }

    for (int i = 0; i < nzPML2 - 1; i++) {
        sigh_zPML2[i] = sig_z_max * pow((nzPML2 - 1.5 - i) / (nzPML2 - 1.0), m);
        kaph_zPML2[i] = 1 + (kap_z_max - 1) * pow((nzPML2 - 1.5 - i) / (nzPML2 - 1.0), m);
        alph_zPML2[i] = alp_z_max * ((nzPML2 - 1.5 - i) / (nzPML2 - 1.0));
        bh_z2[i] = exp(-(sigh_zPML2[i] / (eps0 * kaph_zPML2[i]) + alph_zPML2[i] / eps0) * dt);
        ch_z2[i] = sigh_zPML2[i] / (sigh_zPML2[i] * kaph_zPML2[i] + kaph_zPML2[i] * kaph_zPML2[i] * alph_zPML2[i]) * (bh_z2[i] - 1);
    }

    // 填充电场系数
    for (int i = 0; i < nxPML1; i++) {
        be_x_tmp[i] = be_x1[i];
        ce_x_tmp[i] = ce_x1[i];
    }
    for (int i = 0; i < nxPML2; i++) {
        be_x_tmp[ihT - nxPML2 + i] = be_x2[nxPML2 - 1 - i];
        ce_x_tmp[ihT - nxPML2 + i] = ce_x2[nxPML2 - 1 - i];
    }

    for (int i = 0; i < nyPML1; i++) {
        be_y_tmp[i] = be_y1[i];
        ce_y_tmp[i] = ce_y1[i];
    }
    for (int i = 0; i < nyPML2; i++) {
        be_y_tmp[jhT - nyPML2 + i] = be_y2[nyPML2 - 1 - i];
        ce_y_tmp[jhT - nyPML2 + i] = ce_y2[nyPML2 - 1 - i];
    }

    for (int i = 0; i < nzPML1; i++) {
        be_z_tmp[i] = be_z1[i];
        ce_z_tmp[i] = ce_z1[i];
    }
    for (int i = 0; i < nzPML2; i++) {
        be_z_tmp[khT - nzPML2 + i] = be_z2[nzPML2 - 1 - i];
        ce_z_tmp[khT - nzPML2 + i] = ce_z2[nzPML2 - 1 - i];
    }

    // 填充磁场系数
    for (int i = 0; i < nxPML1 - 1; i++) {
        bh_x_tmp[i] = bh_x1[i];
        ch_x_tmp[i] = ch_x1[i];
    }
    for (int i = 0; i < nxPML2 - 1; i++) {
        bh_x_tmp[ieT - nxPML2 + i] = bh_x2[nxPML2 - 2 - i];
        ch_x_tmp[ieT - nxPML2 + i] = ch_x2[nxPML2 - 2 - i];
    }

    for (int i = 0; i < nyPML1 - 1; i++) {
        bh_y_tmp[i] = bh_y1[i];
        ch_y_tmp[i] = ch_y1[i];
    }
    for (int i = 0; i < nyPML2 - 1; i++) {
        bh_y_tmp[jeT - nyPML2 + i] = bh_y2[nyPML2 - 2 - i];
        ch_y_tmp[jeT - nyPML2 + i] = ch_y2[nyPML2 - 2 - i];
    }

    for (int i = 0; i < nzPML1 - 1; i++) {
        bh_z_tmp[i] = bh_z1[i];
        ch_z_tmp[i] = ch_z1[i];
    }
    for (int i = 0; i < nzPML2 - 1; i++) {
        bh_z_tmp[keT - nzPML2 + i] = bh_z2[nzPML2 - 2 - i];
        ch_z_tmp[keT - nzPML2 + i] = ch_z2[nzPML2 - 2 - i];
    }
}

// 释放CPML变量内存
void free_cpml_variables() {
    // 释放CPML辅助变量
    if (psi_Exy != NULL) {
        free3DArray(psi_Exy, ieT, jeT);
        free3DArray(psi_Exz, ieT, jeT);
        free3DArray(psi_Eyx, ieT, jeT);
        free3DArray(psi_Eyz, ieT, jeT);
        free3DArray(psi_Ezx, ieT, jeT);
        free3DArray(psi_Ezy, ieT, jeT);
    }

    if (psi_Hxy != NULL) {
        free3DArray(psi_Hxy, ihT, jhT);
        free3DArray(psi_Hxz, ihT, jhT);
        free3DArray(psi_Hyx, ihT, jhT);
        free3DArray(psi_Hyz, ihT, jhT);
        free3DArray(psi_Hzx, ihT, jhT);
        free3DArray(psi_Hzy, ihT, jhT);
    }

    // 释放系数数组
    free1DArray(be_x1);
    free1DArray(ce_x1);
    free1DArray(sige_xPML1);
    free1DArray(kape_xPML1);
    free1DArray(alpe_xPML1);

    free1DArray(be_x2);
    free1DArray(ce_x2);
    free1DArray(sige_xPML2);
    free1DArray(kape_xPML2);
    free1DArray(alpe_xPML2);

    free1DArray(be_y1);
    free1DArray(ce_y1);
    free1DArray(sige_yPML1);
    free1DArray(kape_yPML1);
    free1DArray(alpe_yPML1);

    free1DArray(be_y2);
    free1DArray(ce_y2);
    free1DArray(sige_yPML2);
    free1DArray(kape_yPML2);
    free1DArray(alpe_yPML2);

    free1DArray(be_z1);
    free1DArray(ce_z1);
    free1DArray(sige_zPML1);
    free1DArray(kape_zPML1);
    free1DArray(alpe_zPML1);

    free1DArray(be_z2);
    free1DArray(ce_z2);
    free1DArray(sige_zPML2);
    free1DArray(kape_zPML2);
    free1DArray(alpe_zPML2);

    free1DArray(bh_x1);
    free1DArray(ch_x1);
    free1DArray(sigh_xPML1);
    free1DArray(kaph_xPML1);
    free1DArray(alph_xPML1);

    free1DArray(bh_x2);
    free1DArray(ch_x2);
    free1DArray(sigh_xPML2);
    free1DArray(kaph_xPML2);
    free1DArray(alph_xPML2);

    free1DArray(bh_y1);
    free1DArray(ch_y1);
    free1DArray(sigh_yPML1);
    free1DArray(kaph_yPML1);
    free1DArray(alph_yPML1);

    free1DArray(bh_y2);
    free1DArray(ch_y2);
    free1DArray(sigh_yPML2);
    free1DArray(kaph_yPML2);
    free1DArray(alph_yPML2);

    free1DArray(bh_z1);
    free1DArray(ch_z1);
    free1DArray(sigh_zPML1);
    free1DArray(kaph_zPML1);
    free1DArray(alph_zPML1);

    free1DArray(bh_z2);
    free1DArray(ch_z2);
    free1DArray(sigh_zPML2);
    free1DArray(kaph_zPML2);
    free1DArray(alph_zPML2);

    // 释放完整的系数数组
    free1DArray(be_x_tmp);
    free1DArray(ce_x_tmp);
    free1DArray(be_y_tmp);
    free1DArray(ce_y_tmp);
    free1DArray(be_z_tmp);
    free1DArray(ce_z_tmp);

    free1DArray(bh_x_tmp);
    free1DArray(ch_x_tmp);
    free1DArray(bh_y_tmp);
    free1DArray(ch_y_tmp);
    free1DArray(bh_z_tmp);
    free1DArray(ch_z_tmp);
}