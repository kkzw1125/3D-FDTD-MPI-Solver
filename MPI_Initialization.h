#pragma once
#ifndef MPI_INITIAIZATION_H
#define MPI_INITIAIZATION_H
//初始化MPI
#include"MPI.h"
#include"Basic_parameters.h"
#include"Basic_functions.h"

//线程序号
//所有的网格
int total_MAX_X, total_MAX_Y, total_MAX_Z;

//各个方向分配的线程
const int MPI_x_direction = 1;
const int MPI_y_direction = 1;
const int MPI_z_direction = 1;

// 每个线程在三维的排布
//int rank_to_id[MPI_x_direction][MPI_y_direction][MPI_z_direction] = { {{0}} };
//int process_index[3] = { 0 };//每个线程的三位坐标

// 结构体定义，记录每个线程实际位置，方便修改参数
typedef struct {
    int x_start;
    int x_end;
    int y_start;
    int y_end;
    int z_start;
    int z_end;
} AxisParameters;

typedef struct {
    int rank_id; // 第几个线程
    AxisParameters Ex;
    AxisParameters Ey;
    AxisParameters Ez;
    AxisParameters Hx;
    AxisParameters Hy;
    AxisParameters Hz;
    int rank_spatial_distribution[MPI_x_direction][MPI_y_direction][MPI_z_direction];//线程空间分布
    int process_index[3];//每个线程的三位坐标
} MPI_Struct;



void MPI_INITIAIZATION( int rank, MPI_Struct* MPI_parameter) {

    int process_index[3] = { 0 };
    MPI_parameter->rank_id = rank;


    for (int i = 0; i < MPI_x_direction; i++) {
        for (int j = 0; j < MPI_y_direction; j++) {
            for (int k = 0; k < MPI_z_direction; k++) {

                MPI_parameter->rank_spatial_distribution[i][j][k] = (i * MPI_y_direction + j) * MPI_z_direction + k;
                //if (rank == 0) {
                //    printf("rank_to_id[%d][%d][%d]=%d\n", i, j, k, rank_to_id[i][j][k]);
                //}

                if (rank == MPI_parameter->rank_spatial_distribution[i][j][k]) {
                    process_index[0] = i;
                    process_index[1] = j;
                    process_index[2] = k;
                }
            }

        }

    }

    MPI_parameter->process_index[0] = process_index[0];
    MPI_parameter->process_index[1] = process_index[1];
    MPI_parameter->process_index[2] = process_index[2];



	total_MAX_X = MAX_X;
	total_MAX_Y = MAX_Y;
	total_MAX_Z = MAX_Z;

    MAX_X = int(ceil((double(total_MAX_X) + MPI_x_direction - 1) / double(MPI_x_direction)));
    MAX_Y = int(ceil((double(total_MAX_Y) + MPI_y_direction - 1) / double(MPI_y_direction)));
    MAX_Z = int(ceil((double(total_MAX_Z) + MPI_z_direction - 1) / double(MPI_z_direction)));


    MPI_parameter->Ex.x_start = (MAX_X - 1) * process_index[0];
    MPI_parameter->Ex.x_end = MPI_parameter->Ex.x_start + MAX_X - 1;
    MPI_parameter->Ex.y_start = (MAX_Y - 1) * process_index[1];
    MPI_parameter->Ex.y_end = MPI_parameter->Ex.y_start + MAX_Y;
    MPI_parameter->Ex.z_start = (MAX_Z - 1) * process_index[2];
    MPI_parameter->Ex.z_end = MPI_parameter->Ex.z_start + MAX_Z;

    //Ey
    MPI_parameter->Ey.x_start = (MAX_X - 1) * process_index[0];
    MPI_parameter->Ey.x_end = MPI_parameter->Ey.x_start + MAX_X;
    MPI_parameter->Ey.y_start = (MAX_Y - 1) * process_index[1];
    MPI_parameter->Ey.y_end = MPI_parameter->Ey.y_start + MAX_Y - 1;
    MPI_parameter->Ey.z_start = (MAX_Z - 1) * process_index[2];
    MPI_parameter->Ey.z_end = MPI_parameter->Ey.z_start + MAX_Z;

    // Ez
    MPI_parameter->Ez.x_start = (MAX_X - 1) * process_index[0];
    MPI_parameter->Ez.x_end = MPI_parameter->Ez.x_start + MAX_X;
    MPI_parameter->Ez.y_start = (MAX_Y - 1) * process_index[1];
    MPI_parameter->Ez.y_end = MPI_parameter->Ez.y_start + MAX_Y;
    MPI_parameter->Ez.z_start = (MAX_Z - 1) * process_index[2];
    MPI_parameter->Ez.z_end = MPI_parameter->Ez.z_start + MAX_Z - 1;

    // Hx
    MPI_parameter->Hx.x_start = (MAX_X - 1) * process_index[0];
    MPI_parameter->Hx.x_end = MPI_parameter->Hx.x_start + MAX_X;
    MPI_parameter->Hx.y_start = (MAX_Y - 1) * process_index[1];
    MPI_parameter->Hx.y_end = MPI_parameter->Hx.y_start + MAX_Y - 1;
    MPI_parameter->Hx.z_start = (MAX_Z - 1) * process_index[2];
    MPI_parameter->Hx.z_end = MPI_parameter->Hx.z_start + MAX_Z - 1;

    // Hy
    MPI_parameter->Hy.x_start = (MAX_X - 1) * process_index[0];
    MPI_parameter->Hy.x_end = MPI_parameter->Hy.x_start + MAX_X - 1;
    MPI_parameter->Hy.y_start = (MAX_Y - 1) * process_index[1];
    MPI_parameter->Hy.y_end = MPI_parameter->Hy.y_start + MAX_Y;
    MPI_parameter->Hy.z_start = (MAX_Z - 1) * process_index[2];
    MPI_parameter->Hy.z_end = MPI_parameter->Hy.z_start + MAX_Z - 1;

    // Hz
    MPI_parameter->Hz.x_start = (MAX_X - 1) * process_index[0];
    MPI_parameter->Hz.x_end = MPI_parameter->Hz.x_start + MAX_X - 1;
    MPI_parameter->Hz.y_start = (MAX_Y - 1) * process_index[1];
    MPI_parameter->Hz.y_end = MPI_parameter->Hz.y_start + MAX_Y - 1;
    MPI_parameter->Hz.z_start = (MAX_Z - 1) * process_index[2];
    MPI_parameter->Hz.z_end = MPI_parameter->Hz.z_start + MAX_Z;

    if (process_index[0] == MPI_x_direction - 1) {
        MPI_parameter->Ex.x_end = total_MAX_X - 1;
        MPI_parameter->Ey.x_end = total_MAX_X;
        MPI_parameter->Ez.x_end = total_MAX_X;
        MPI_parameter->Hx.x_end = total_MAX_X;
        MPI_parameter->Hy.x_end = total_MAX_X - 1;
        MPI_parameter->Hz.x_end = total_MAX_X - 1;
        MAX_X = MPI_parameter->Ex.x_end - MPI_parameter->Ex.x_start + 1;
    }
    if (process_index[1] == MPI_y_direction - 1) {
        MPI_parameter->Ex.y_end = total_MAX_Y;
        MPI_parameter->Ey.y_end = total_MAX_Y - 1;
        MPI_parameter->Ez.y_end = total_MAX_Y;
        MPI_parameter->Hx.y_end = total_MAX_Y - 1;
        MPI_parameter->Hy.y_end = total_MAX_Y;
        MPI_parameter->Hz.y_end = total_MAX_Y - 1;
        MAX_Y = MPI_parameter->Ey.y_end - MPI_parameter->Ey.y_start + 1;
    }

    if (process_index[2] == MPI_z_direction - 1) {
        MPI_parameter->Ex.z_end = total_MAX_Z;
        MPI_parameter->Ey.z_end = total_MAX_Z;
        MPI_parameter->Ez.z_end = total_MAX_Z - 1;
        MPI_parameter->Hx.z_end = total_MAX_Z - 1;
        MPI_parameter->Hy.z_end = total_MAX_Z - 1;
        MPI_parameter->Hz.z_end = total_MAX_Z;
        MAX_Z = MPI_parameter->Ez.z_end - MPI_parameter->Ez.z_start + 1;
    }

}



#endif