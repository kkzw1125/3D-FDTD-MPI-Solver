#pragma once
#ifndef Calculate_magnetic_filed_H
#define Calculate_magnetic_filed_H


void Calculate_magnetic_filed_Hx() {

	for (int i = 0; i < MAX_X + 1; i++) {
		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z; k++) {

				/*Hx[i][j][k] = Hx[i][j][k] + CQx[i][j][k] * (Ey[i][j][k + 1] - Ey[i][j][k]) / dz\
					- CQx[i][j][k] * (Ez[i][j + 1][k] - Ez[i][j][k]) / dy;*/
				psi_Hxy[i][j][k] = bh_y_tmp[j] * psi_Hxy[i][j][k] + ch_y_tmp[j] * (Ez[i][j + 1][k] - Ez[i][j][k]) / dy;
				psi_Hxz[i][j][k] = bh_z_tmp[k] * psi_Hxz[i][j][k] + ch_z_tmp[k] * (Ey[i][j][k + 1] - Ey[i][j][k]) / dz;
				Hx[i][j][k] = CP * Hx[i][j][k] - CQx[i][j][k] * (Ez[i][j + 1][k] - Ez[i][j][k]) / dy\
					+ CQx[i][j][k] * (Ey[i][j][k + 1] - Ey[i][j][k]) / dz\
					- CQx[i][j][k] * (psi_Hxy[i][j][k] - psi_Hxz[i][j][k]);
			}
		}
	}

	/*for (int i = i0 - 1; i < ia; i++) {
		for (int k = k0 - 1; k < kc - 1; k++) {
			Hx[i][j_0 - 2][k] = Hx[i][j_0 - 2][k] + CQy[1][1][1] * ez_inc[j_0 - 1];
			Hx[i][jb - 1][k] = Hx[i][jb - 1][k] - CQy[1][1][1] * ez_inc[jb - 1];
		}
	}*/

}

void Calculate_magnetic_filed_Hy() {

	for (int i = 0; i < MAX_X; i++) {
		for (int j = 0; j < MAX_Y + 1; j++) {
			for (int k = 0; k < MAX_Z; k++) {

				/*Hy[i][j][k] = Hy[i][j][k] + CQy[i][j][k] * (Ez[i + 1][j][k] - Ez[i][j][k]) / dx\
					- CQy[i][j][k] * (Ex[i][j][k + 1] - Ex[i][j][k]) / dz;*/
				psi_Hyx[i][j][k] = bh_x_tmp[i] * psi_Hyx[i][j][k] + ch_x_tmp[i] * (Ez[i + 1][j][k] - Ez[i][j][k]) / dx;
				psi_Hyz[i][j][k] = bh_z_tmp[k] * psi_Hyz[i][j][k] + ch_z_tmp[k] * (Ex[i][j][k + 1] - Ex[i][j][k]) / dz;
				Hy[i][j][k] = CP * Hy[i][j][k] - CQy[i][j][k] * (Ex[i][j][k + 1] - Ex[i][j][k]) / dz\
					+ CQy[i][j][k] * (Ez[i + 1][j][k] - Ez[i][j][k]) / dx\
					- CQy[i][j][k] * (psi_Hyz[i][j][k] - psi_Hyx[i][j][k]);
			}
		}
	}

	/*for (int j = j_0 - 1; j < jb; j++) {
		for (int k = k0 - 1; k < kc - 1; k++) {
			Hy[i0 - 2][j][k] = Hy[i0 - 2][j][k] - CQy[1][1][1] * ez_inc[j];
			Hy[ia - 1][j][k] = Hy[ia - 1][j][k] + CQy[1][1][1] * ez_inc[j];
		}
	}*/

}

void Calculate_magnetic_filed_Hz() {

	for (int i = 0; i < MAX_X; i++) {
		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z + 1; k++) {

				/*Hz[i][j][k] = Hz[i][j][k] + CQz[i][j][k] * (Ex[i][j + 1][k] - Ex[i][j][k]) / dy\
					- CQz[i][j][k] * (Ey[i + 1][j][k] - Ey[i][j][k]) / dx;*/
				psi_Hzx[i][j][k] = bh_x_tmp[i] * psi_Hzx[i][j][k] + ch_x_tmp[i] * (Ey[i + 1][j][k] - Ey[i][j][k]) / dx;
				psi_Hzy[i][j][k] = bh_y_tmp[j] * psi_Hzy[i][j][k] + ch_y_tmp[j] * (Ex[i][j + 1][k] - Ex[i][j][k]) / dy;
				Hz[i][j][k] = CP * Hz[i][j][k] - CQz[i][j][k] * (Ey[i + 1][j][k] - Ey[i][j][k]) / dx\
					+ CQz[i][j][k] * (Ex[i][j + 1][k] - Ex[i][j][k]) / dy\
					- CQz[i][j][k] * (psi_Hzx[i][j][k] - psi_Hzy[i][j][k]);
			}
		}
	}

}

void Calculate_magnetic_filed() {
	Calculate_magnetic_filed_Hx();
	Calculate_magnetic_filed_Hy();
	Calculate_magnetic_filed_Hz();
}


#endif