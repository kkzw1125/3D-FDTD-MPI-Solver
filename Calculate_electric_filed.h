#pragma once
#ifndef Calculate_electric_filed_H
#define Calculate_electric_filed_H
#include "Basic_parameters.h"
#include "CPML.h"

void Calculate_electric_filed_Ex() {

	for (int i = 0; i < MAX_X; i++) {
		for (int j = 1; j < MAX_Y; j++) {
			for (int k = 1; k < MAX_Z; k++) {
				/*Ex[i][j][k] = Ex[i][j][k] + CBx[i][j][k] * (Hz[i][j][k] - Hz[i][j - 1][k]) / dy\
					- CBx[i][j][k] * (Hy[i][j][k] - Hy[i][j][k - 1]) / dz;*/
				psi_Exy[i][j][k] = be_y_tmp[j] * psi_Exy[i][j][k] + ce_y_tmp[j] * (Hz[i][j][k] - Hz[i][j - 1][k]) / dy;
				psi_Exz[i][j][k] = be_z_tmp[k] * psi_Exz[i][j][k] + ce_z_tmp[k] * (Hy[i][j][k] - Hy[i][j][k - 1]) / dz;
				Ex[i][j][k] = CA * Ex[i][j][k] + CBx[i][j][k] * (Hz[i][j][k] - Hz[i][j - 1][k]) / dy\
					- CBx[i][j][k] * (Hy[i][j][k] - Hy[i][j][k - 1]) / dz\
					+ CBx[i][j][k] * (psi_Exy[i][j][k] - psi_Exz[i][j][k]);
			}
		}
	}
}

void Calculate_electric_filed_Ey() {

	for (int i = 1; i < MAX_X; i++) {
		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 1; k < MAX_Z; k++) {
				/*Ey[i][j][k] = Ey[i][j][k] + CBy[i][j][k] * (Hx[i][j][k] - Hx[i][j][k - 1]) / dz\
					- CBy[i][j][k] * (Hz[i][j][k] - Hz[i - 1][j][k]) / dx;*/
				psi_Eyx[i][j][k] = be_x_tmp[i] * psi_Eyx[i][j][k] + ce_x_tmp[i] * (Hz[i][j][k] - Hz[i - 1][j][k]) / dx;
				psi_Eyz[i][j][k] = be_z_tmp[k] * psi_Eyz[i][j][k] + ce_z_tmp[k] * (Hx[i][j][k] - Hx[i][j][k - 1]) / dz;
				Ey[i][j][k] = CA * Ey[i][j][k] + CBy[i][j][k] * (Hx[i][j][k] - Hx[i][j][k - 1]) / dz\
					- CBy[i][j][k] * (Hz[i][j][k] - Hz[i - 1][j][k]) / dx\
					+ CBy[i][j][k] * (psi_Eyz[i][j][k] - psi_Eyx[i][j][k]);
			}
		}
	}


	/*// ФґПо
	for (int i = i0 - 1; i < ia; i++) {
		for (int j = j_0 - 1; j < jb - 1; j++) {
			Ey[i][j][k0 - 1] = Ey[i][j][k0 - 1] - CBy[1][j][1] * hx_inc[j];
			Ey[i][j][kc - 1] = Ey[i][j][kc - 1] + CBy[1][j][1] * hx_inc[j];
		}
	}*/

}


void Calculate_electric_filed_Ez() {

	for (int i = 1; i < MAX_X; i++) {
		for (int j = 1; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z; k++) {

				/*Ez[i][j][k] = Ez[i][j][k] + CBz[i][j][k] * (Hy[i][j][k] - Hy[i - 1][j][k]) / dx\
					- CBz[i][j][k] * (Hx[i][j][k] - Hx[i][j - 1][k]) / dy;*/
				psi_Ezx[i][j][k] = be_x_tmp[i] * psi_Ezx[i][j][k] + ce_x_tmp[i] * (Hy[i][j][k] - Hy[i - 1][j][k]) / dx;
				psi_Ezy[i][j][k] = be_y_tmp[j] * psi_Ezy[i][j][k] + ce_y_tmp[j] * (Hx[i][j][k] - Hx[i][j - 1][k]) / dy;
				Ez[i][j][k] = CA * Ez[i][j][k] + CBz[i][j][k] * (Hy[i][j][k] - Hy[i - 1][j][k]) / dx\
					- CBz[i][j][k] * (Hx[i][j][k] - Hx[i][j - 1][k]) / dy\
					+ CBz[i][j][k] * (psi_Ezx[i][j][k] - psi_Ezy[i][j][k]);
			}
		}
	}

	/*for (int i = i0 - 1; i < ia; i++) {
		for (int k = k0 - 1; k < kc - 1; k++) {
			Ez[i][j_0 - 1][k] = Ez[i][j_0 - 1][k] + CBy[1][1][1] * hx_inc[j_0 - 2];
			Ez[i][jb - 1][k] = Ez[i][jb - 1][k] - CBy[1][1][1] * hx_inc[jb - 1];
		}
	}*/

}

void Calculate_electric_filed() {
	Calculate_electric_filed_Ex();
	Calculate_electric_filed_Ey();
	Calculate_electric_filed_Ez();

}






#endif
