#pragma once
#ifndef FDTD_MPI_FUNCTIONS_H
#define FDTD_MPI_FUNCTIONS_H

#include"MPI_Initialization.h"
#include"Basic_parameters.h"
#include"Basic_functions.h"


void MPI_Ex_y_direction_Change(MPI_Struct* MPI_parameter) {

	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { 0 };
	
	double** Ex_send_y_last = createContinuous2DArray(MAX_X, MAX_Z + 1, 0.0);
	double** Ex_recv_y_last = createContinuous2DArray(MAX_X, MAX_Z + 1, 0.0);
	double** Ex_send_y_next = createContinuous2DArray(MAX_X, MAX_Z + 1, 0.0);
	double** Ex_recv_y_next = createContinuous2DArray(MAX_X, MAX_Z + 1, 0.0);

	int number_transmit = MAX_X * (MAX_Z + 1);

	process_index[0] = MPI_parameter->process_index[0];
	process_index[1] = MPI_parameter->process_index[1];
	process_index[2] = MPI_parameter->process_index[2];

	// 发送请求以及确定
	MPI_Request requests[4]{ MPI_REQUEST_NULL };
	MPI_Status statuses[4];
	int count = 0;

	if (process_index[1] >= 0 && process_index[1] < MPI_y_direction - 1) {
		int y_direction_next_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] + 1][process_index[2]];

		for (int i = 0; i < MAX_X; i++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ex_send_y_next[i][k] = Ex[i][MAX_Y - 1][k];
			}
		}
					
		MPI_Isend(Ex_send_y_next[0], number_transmit, MPI_DOUBLE, y_direction_next_rank, y_direction_next_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ex_recv_y_next[0], number_transmit, MPI_DOUBLE, y_direction_next_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}


	if (process_index[1] > 0 && process_index[1] <= MPI_y_direction - 1) {

		int y_direction_last_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]];

		//printf("y-directiong last rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]]);

		for (int i = 0; i < MAX_X; i++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ex_send_y_last[i][k] = Ex[i][1][k];
			}
		}
		

		MPI_Isend(Ex_send_y_last[0], number_transmit, MPI_DOUBLE, y_direction_last_rank, y_direction_last_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ex_recv_y_last[0], number_transmit, MPI_DOUBLE, y_direction_last_rank, rankid, MPI_COMM_WORLD, &requests[count++]);
		
	}


	//printf("count=%d\n", count);

	if (count > 0) {				
		MPI_Waitall(count, requests, statuses);
	}
	
	if (process_index[1] >= 0 && process_index[1] < MPI_y_direction - 1) {

		for (int i = 0; i < MAX_X; i++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ex[i][MAX_Y][k] = Ex_recv_y_next[i][k];
			}
		}

	}
	if (process_index[1] > 0 && process_index[1] <= MPI_y_direction - 1) {

		for (int i = 0; i < MAX_X; i++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ex[i][0][k] = Ex_recv_y_last[i][k];
			}
		}

	}



	free(Ex_recv_y_next[0]);  // 释放存储数据的连续内存块
	free(Ex_recv_y_next);
	free(Ex_send_y_next[0]);  // 释放存储数据的连续内存块
	free(Ex_send_y_next);
	free(Ex_recv_y_last[0]);  // 释放存储数据的连续内存块
	free(Ex_recv_y_last);
	free(Ex_send_y_last[0]);  // 释放存储数据的连续内存块
	free(Ex_send_y_last);


	//printf("---------------------\n");
	
	
}

void MPI_Ex_z_direction_Change(MPI_Struct* MPI_parameter) {

	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { 0 };

	double** Ex_send_z_last = createContinuous2DArray(MAX_X, MAX_Y + 1, 0.0);
	double** Ex_recv_z_last = createContinuous2DArray(MAX_X, MAX_Y + 1, 0.0);
	double** Ex_send_z_next = createContinuous2DArray(MAX_X, MAX_Y + 1, 0.0);
	double** Ex_recv_z_next = createContinuous2DArray(MAX_X, MAX_Y + 1, 0.0);

	int number_transmit = MAX_X * (MAX_Y + 1);

	process_index[0] = MPI_parameter->process_index[0];
	process_index[1] = MPI_parameter->process_index[1];
	process_index[2] = MPI_parameter->process_index[2];

	// 发送请求以及确定
	MPI_Request requests[4]{ MPI_REQUEST_NULL };
	MPI_Status statuses[4];
	int count = 0;

	if (process_index[2] >= 0 && process_index[2] < MPI_z_direction - 1) {
		int z_direction_next_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1]][process_index[2] + 1];

		for (int i = 0; i < MAX_X; i++) {
			for (int j = 0; j < MAX_Y + 1; j++) {
				Ex_send_z_next[i][j] = Ex[i][j][MAX_Z - 1];
			}
		}

		MPI_Isend(Ex_send_z_next[0], number_transmit, MPI_DOUBLE, z_direction_next_rank, z_direction_next_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ex_recv_z_next[0], number_transmit, MPI_DOUBLE, z_direction_next_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}


	if (process_index[2] > 0 && process_index[2] <= MPI_z_direction - 1) {

		int z_direction_last_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1]][process_index[2] - 1];

		//printf("y-directiong last rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]]);

		for (int i = 0; i < MAX_X; i++) {
			for (int j = 0; j < MAX_Y + 1; j++) {
				Ex_send_z_last[i][j] = Ex[i][j][1];
			}
		}


		MPI_Isend(Ex_send_z_last[0], number_transmit, MPI_DOUBLE, z_direction_last_rank, z_direction_last_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ex_recv_z_last[0], number_transmit, MPI_DOUBLE, z_direction_last_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}

	if (count > 0) {
		MPI_Waitall(count, requests, statuses);
	}

	if (process_index[2] >= 0 && process_index[2] < MPI_z_direction - 1) {

		for (int i = 0; i < MAX_X; i++) {
			for (int j = 0; j < MAX_Y + 1; j++) {
				Ex[i][j][MAX_Z] = Ex_recv_z_next[i][j];
			}
		}

	}
	if (process_index[2] > 0 && process_index[2] <= MPI_z_direction - 1) {

		for (int i = 0; i < MAX_X; i++) {
			for (int j = 0; j < MAX_Y + 1; j++) {
				Ex[i][j][0] = Ex_recv_z_last[i][j];
			}
		}

	}

	free(Ex_recv_z_next[0]);  // 释放存储数据的连续内存块
	free(Ex_recv_z_next);
	free(Ex_send_z_next[0]);  // 释放存储数据的连续内存块
	free(Ex_send_z_next);
	free(Ex_recv_z_last[0]);  // 释放存储数据的连续内存块
	free(Ex_recv_z_last);
	free(Ex_send_z_last[0]);  // 释放存储数据的连续内存块
	free(Ex_send_z_last);


}

void MPI_Ey_x_direction_Change(MPI_Struct* MPI_parameter) {

	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { 0 };

	double** Ey_send_x_last = createContinuous2DArray(MAX_Y, MAX_Z + 1, 0.0);
	double** Ey_recv_x_last = createContinuous2DArray(MAX_Y, MAX_Z + 1, 0.0);
	double** Ey_send_x_next = createContinuous2DArray(MAX_Y, MAX_Z + 1, 0.0);
	double** Ey_recv_x_next = createContinuous2DArray(MAX_Y, MAX_Z + 1, 0.0);

	int number_transmit = MAX_Y * (MAX_Z + 1);

	process_index[0] = MPI_parameter->process_index[0];
	process_index[1] = MPI_parameter->process_index[1];
	process_index[2] = MPI_parameter->process_index[2];

	// 发送请求以及确定
	MPI_Request requests[4]{ MPI_REQUEST_NULL };
	MPI_Status statuses[4];
	int count = 0;

	if (process_index[0] >= 0 && process_index[0] < MPI_x_direction - 1) {
		int x_direction_next_rank = MPI_parameter->rank_spatial_distribution[process_index[0] + 1][process_index[1]][process_index[2]];

		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ey_send_x_next[j][k] = Ey[MAX_X - 1 ][j][k];
			}
		}

		MPI_Isend(Ey_send_x_next[0], number_transmit, MPI_DOUBLE, x_direction_next_rank, x_direction_next_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ey_recv_x_next[0], number_transmit, MPI_DOUBLE, x_direction_next_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}


	if (process_index[0] > 0 && process_index[0] <= MPI_x_direction - 1) {

		int x_direction_last_rank = MPI_parameter->rank_spatial_distribution[process_index[0] - 1][process_index[1]][process_index[2]];

		//printf("y-directiong last rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]]);

		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ey_send_x_last[j][k] = Ey[1][j][k];
			}
		}


		MPI_Isend(Ey_send_x_last[0], number_transmit, MPI_DOUBLE, x_direction_last_rank, x_direction_last_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ey_recv_x_last[0], number_transmit, MPI_DOUBLE, x_direction_last_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}

	if (count > 0) {
		MPI_Waitall(count, requests, statuses);
	}

	if (process_index[0] >= 0 && process_index[0] < MPI_x_direction - 1) {

		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ey[MAX_X][j][k] = Ey_recv_x_next[j][k];
			}
		}

	}
	if (process_index[0] > 0 && process_index[0] <= MPI_x_direction - 1) {

		for (int j = 0; j < MAX_Y; j++) {
			for (int k = 0; k < MAX_Z + 1; k++) {
				Ey[0][j][k] = Ey_recv_x_last[j][k];
			}
		}

	}

	free(Ey_recv_x_next[0]);  // 释放存储数据的连续内存块
	free(Ey_recv_x_next);
	free(Ey_send_x_next[0]);  // 释放存储数据的连续内存块
	free(Ey_send_x_next);
	free(Ey_recv_x_last[0]);  // 释放存储数据的连续内存块
	free(Ey_recv_x_last);
	free(Ey_send_x_last[0]);  // 释放存储数据的连续内存块
	free(Ey_send_x_last);


}

void MPI_Ey_z_direction_Change(MPI_Struct* MPI_parameter) {

	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { 0 };

	double** Ey_send_z_last = createContinuous2DArray(MAX_X + 1, MAX_Y, 0.0);
	double** Ey_recv_z_last = createContinuous2DArray(MAX_X + 1, MAX_Y, 0.0);
	double** Ey_send_z_next = createContinuous2DArray(MAX_X + 1, MAX_Y, 0.0);
	double** Ey_recv_z_next = createContinuous2DArray(MAX_X + 1, MAX_Y, 0.0);

	int number_transmit = (MAX_X + 1) * MAX_Y;

	process_index[0] = MPI_parameter->process_index[0];
	process_index[1] = MPI_parameter->process_index[1];
	process_index[2] = MPI_parameter->process_index[2];

	// 发送请求以及确定
	MPI_Request requests[4]{ MPI_REQUEST_NULL };
	MPI_Status statuses[4];
	int count = 0;

	if (process_index[2] >= 0 && process_index[2] < MPI_z_direction - 1) {
		int z_direction_next_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1]][process_index[2] + 1];

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int j = 0; j < MAX_Y; j++) {
				Ey_send_z_next[i][j] = Ey[i][j][MAX_Z - 1];
			}
		}

		MPI_Isend(Ey_send_z_next[0], number_transmit, MPI_DOUBLE, z_direction_next_rank, z_direction_next_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ey_recv_z_next[0], number_transmit, MPI_DOUBLE, z_direction_next_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}


	if (process_index[2] > 0 && process_index[2] <= MPI_z_direction - 1) {

		int z_direction_last_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1]][process_index[2] - 1];

		//printf("y-directiong last rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]]);

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int j = 0; j < MAX_Y; j++) {
				Ey_send_z_last[i][j] = Ey[i][j][1];
			}
		}


		MPI_Isend(Ey_send_z_last[0], number_transmit, MPI_DOUBLE, z_direction_last_rank, z_direction_last_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ey_recv_z_last[0], number_transmit, MPI_DOUBLE, z_direction_last_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}

	if (count > 0) {
		MPI_Waitall(count, requests, statuses);
	}

	if (process_index[2] >= 0 && process_index[2] < MPI_z_direction - 1) {

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int j = 0; j < MAX_Y; j++) {
				Ey[i][j][MAX_Z] = Ey_recv_z_next[i][j];
			}
		}

	}
	if (process_index[2] > 0 && process_index[2] <= MPI_z_direction - 1) {

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int j = 0; j < MAX_Y; j++) {
				Ey[i][j][0] = Ey_recv_z_last[i][j];
			}
		}

	}

	free(Ey_recv_z_next[0]);  // 释放存储数据的连续内存块
	free(Ey_recv_z_next);
	free(Ey_send_z_next[0]);  // 释放存储数据的连续内存块
	free(Ey_send_z_next);
	free(Ey_recv_z_last[0]);  // 释放存储数据的连续内存块
	free(Ey_recv_z_last);
	free(Ey_send_z_last[0]);  // 释放存储数据的连续内存块
	free(Ey_send_z_last);


}

void MPI_Ez_x_direction_Change(MPI_Struct* MPI_parameter) {

	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { 0 };

	double** Ez_send_x_last = createContinuous2DArray(MAX_Y + 1, MAX_Z, 0.0);
	double** Ez_recv_x_last = createContinuous2DArray(MAX_Y + 1, MAX_Z, 0.0);
	double** Ez_send_x_next = createContinuous2DArray(MAX_Y + 1, MAX_Z, 0.0);
	double** Ez_recv_x_next = createContinuous2DArray(MAX_Y + 1, MAX_Z, 0.0);

	int number_transmit = (MAX_Y + 1) * (MAX_Z);

	process_index[0] = MPI_parameter->process_index[0];
	process_index[1] = MPI_parameter->process_index[1];
	process_index[2] = MPI_parameter->process_index[2];

	// 发送请求以及确定
	MPI_Request requests[4]{ MPI_REQUEST_NULL };
	MPI_Status statuses[4];
	int count = 0;

	if (process_index[0] >= 0 && process_index[0] < MPI_x_direction - 1) {
		int x_direction_next_rank = MPI_parameter->rank_spatial_distribution[process_index[0] + 1][process_index[1]][process_index[2]];

		for (int j = 0; j < MAX_Y + 1; j++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez_send_x_next[j][k] = Ez[MAX_X - 1 ][j][k];
			}
		}

		MPI_Isend(Ez_send_x_next[0], number_transmit, MPI_DOUBLE, x_direction_next_rank, x_direction_next_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ez_recv_x_next[0], number_transmit, MPI_DOUBLE, x_direction_next_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}


	if (process_index[0] > 0 && process_index[0] <= MPI_x_direction - 1) {

		int x_direction_last_rank = MPI_parameter->rank_spatial_distribution[process_index[0] - 1][process_index[1]][process_index[2]];

		//printf("y-directiong last rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]]);

		for (int j = 0; j < MAX_Y + 1; j++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez_send_x_last[j][k] = Ez[1][j][k];
			}
		}


		MPI_Isend(Ez_send_x_last[0], number_transmit, MPI_DOUBLE, x_direction_last_rank, x_direction_last_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ez_recv_x_last[0], number_transmit, MPI_DOUBLE, x_direction_last_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}

	if (count > 0) {
		MPI_Waitall(count, requests, statuses);
	}

	if (process_index[0] >= 0 && process_index[0] < MPI_x_direction - 1) {

		for (int j = 0; j < MAX_Y + 1; j++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez[MAX_X][j][k] = Ez_recv_x_next[j][k];
			}
		}

	}
	if (process_index[0] > 0 && process_index[0] <= MPI_x_direction - 1) {

		for (int j = 0; j < MAX_Y + 1; j++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez[0][j][k] = Ez_recv_x_last[j][k];
			}
		}

	}

	free(Ez_recv_x_next[0]);  // 释放存储数据的连续内存块
	free(Ez_recv_x_next);
	free(Ez_send_x_next[0]);  // 释放存储数据的连续内存块
	free(Ez_send_x_next);
	free(Ez_recv_x_last[0]);  // 释放存储数据的连续内存块
	free(Ez_recv_x_last);
	free(Ez_send_x_last[0]);  // 释放存储数据的连续内存块
	free(Ez_send_x_last);


}

void MPI_Ez_y_direction_Change(MPI_Struct* MPI_parameter) {

	//printf("---------Ez----------\n");
	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { 0 };

	int a = 3; int b = 3;

	double** Ez_send_y_last = createContinuous2DArray(MAX_X + 1, MAX_Z, 0.0);
	double** Ez_recv_y_last = createContinuous2DArray(MAX_X + 1, MAX_Z, 0.0);
	double** Ez_send_y_next = createContinuous2DArray(MAX_X + 1, MAX_Z, 0.0);
	double** Ez_recv_y_next = createContinuous2DArray(MAX_X + 1, MAX_Z, 0.0);



	int number_transmit = (MAX_X + 1) * (MAX_Z);




	process_index[0] = MPI_parameter->process_index[0];
	process_index[1] = MPI_parameter->process_index[1];
	process_index[2] = MPI_parameter->process_index[2];

	// 发送请求以及确定
	MPI_Request requests[4]{ MPI_REQUEST_NULL };
	MPI_Status statuses[4];
	int count = 0;

	if (process_index[1] >= 0 && process_index[1] < MPI_y_direction - 1) {
		int y_direction_next_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] + 1][process_index[2]];
		//printf("y-directiong next rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] + 1][process_index[2]]);

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez_send_y_next[i][k] = Ez[i][MAX_Y - 1][k];
			}
		}


		MPI_Isend(Ez_send_y_next[0], number_transmit, MPI_DOUBLE, y_direction_next_rank, y_direction_next_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ez_recv_y_next[0], number_transmit, MPI_DOUBLE, y_direction_next_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}
	


	if (process_index[1] > 0 && process_index[1] <= MPI_y_direction - 1) {

		int y_direction_last_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]];

		//printf("y-directiong last rank:%d\n", MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]]);

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez_send_y_last[i][k] = Ez[i][1][k];
			}
		}


		MPI_Isend(Ez_send_y_last[0], number_transmit, MPI_DOUBLE, y_direction_last_rank, y_direction_last_rank, MPI_COMM_WORLD, &requests[count++]);
		MPI_Irecv(Ez_recv_y_last[0], number_transmit, MPI_DOUBLE, y_direction_last_rank, rankid, MPI_COMM_WORLD, &requests[count++]);

	}
	


	//printf("count=%d\n", count);

	if (count > 0) {

		MPI_Waitall(count, requests, statuses);
		//printf("rank:%d finsh transmit,count = %d\n", rankid, count);
	}
	//MPI_Barrier(MPI_COMM_WORLD);

	if (process_index[1] >= 0 && process_index[1] < MPI_y_direction - 1) {

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez[i][MAX_Y][k] = Ez_recv_y_next[i][k];
			}
		}

	}
	if (process_index[1] > 0 && process_index[1] <= MPI_y_direction - 1) {

		for (int i = 0; i < MAX_X + 1; i++) {
			for (int k = 0; k < MAX_Z; k++) {
				Ez[i][0][k] = Ez_recv_y_last[i][k];
			}
		}

	}



	free(Ez_recv_y_next[0]);  // 释放存储数据的连续内存块
	free(Ez_recv_y_next);
	free(Ez_send_y_next[0]);  // 释放存储数据的连续内存块
	free(Ez_send_y_next);
	free(Ez_recv_y_last[0]);  // 释放存储数据的连续内存块
	free(Ez_recv_y_last);
	free(Ez_send_y_last[0]);  // 释放存储数据的连续内存块
	free(Ez_send_y_last);


	//printf("---------------------\n");
	

}


void MPI_E_Change(MPI_Struct* MPI_parameter) {

	MPI_Ex_y_direction_Change(MPI_parameter);
	MPI_Ex_z_direction_Change(MPI_parameter);
	MPI_Ey_x_direction_Change(MPI_parameter);
	MPI_Ey_z_direction_Change(MPI_parameter);
	MPI_Ez_x_direction_Change(MPI_parameter);
	MPI_Ez_y_direction_Change(MPI_parameter);

}


// =========================================================================
//                          Magnetic Field Exchange
// =========================================================================

// 1. Hx 在 Y 方向交换 (需要 Hx 在 Y 维有 Ghost Cell)
void MPI_Hx_y_direction_Change(MPI_Struct* MPI_parameter) {
	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { MPI_parameter->process_index[0], MPI_parameter->process_index[1], MPI_parameter->process_index[2] };

	// Hx 的 X 维是 MAX_X+1，Z 维是 MAX_Z (建议扩充为 MAX_Z+1 以防万一)
	int rows = MAX_X + 1;
	int cols = MAX_Z;

	double** send_next = createContinuous2DArray(rows, cols, 0.0);
	double** recv_next = createContinuous2DArray(rows, cols, 0.0);
	double** send_last = createContinuous2DArray(rows, cols, 0.0);
	double** recv_last = createContinuous2DArray(rows, cols, 0.0);
	int num = rows * cols;

	MPI_Request reqs[4]{ MPI_REQUEST_NULL }; MPI_Status stats[4]; int c = 0;

	// 向 Y+1 发送 (取 MAX_Y-1 的数据)
	if (process_index[1] < MPI_y_direction - 1) {
		int next_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] + 1][process_index[2]];
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++)
			send_next[i][k] = Hx[i][MAX_Y - 1][k];

		MPI_Isend(send_next[0], num, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_next[0], num, MPI_DOUBLE, next_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}

	// 向 Y-1 发送 (取索引 1 的数据)
	if (process_index[1] > 0) {
		int last_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1] - 1][process_index[2]];
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++)
			send_last[i][k] = Hx[i][1][k];

		MPI_Isend(send_last[0], num, MPI_DOUBLE, last_rank, last_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_last[0], num, MPI_DOUBLE, last_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}

	if (c > 0) MPI_Waitall(c, reqs, stats);

	// 填回 Ghost Cell
	if (process_index[1] < MPI_y_direction - 1)
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) Hx[i][MAX_Y][k] = recv_next[i][k];

	if (process_index[1] > 0)
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) Hx[i][0][k] = recv_last[i][k];

	free(recv_next[0]); free(recv_next); free(send_next[0]); free(send_next);
	free(recv_last[0]); free(recv_last); free(send_last[0]); free(send_last);
}

// 2. Hx 在 Z 方向交换
void MPI_Hx_z_direction_Change(MPI_Struct* MPI_parameter) {
	int rankid = MPI_parameter->rank_id;
	int process_index[3] = { MPI_parameter->process_index[0], MPI_parameter->process_index[1], MPI_parameter->process_index[2] };

	int rows = MAX_X + 1;
	int cols = MAX_Y; // Hx 的 Y 维度 (注意检查 Basic_parameters.h)

	double** send_next = createContinuous2DArray(rows, cols, 0.0);
	double** recv_next = createContinuous2DArray(rows, cols, 0.0);
	double** send_last = createContinuous2DArray(rows, cols, 0.0);
	double** recv_last = createContinuous2DArray(rows, cols, 0.0);
	int num = rows * cols;

	MPI_Request reqs[4]{ MPI_REQUEST_NULL }; MPI_Status stats[4]; int c = 0;

	if (process_index[2] < MPI_z_direction - 1) {
		int next_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1]][process_index[2] + 1];
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) send_next[i][j] = Hx[i][j][MAX_Z - 1];
		MPI_Isend(send_next[0], num, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_next[0], num, MPI_DOUBLE, next_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (process_index[2] > 0) {
		int last_rank = MPI_parameter->rank_spatial_distribution[process_index[0]][process_index[1]][process_index[2] - 1];
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) send_last[i][j] = Hx[i][j][1];
		MPI_Isend(send_last[0], num, MPI_DOUBLE, last_rank, last_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_last[0], num, MPI_DOUBLE, last_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (c > 0) MPI_Waitall(c, reqs, stats);

	if (process_index[2] < MPI_z_direction - 1)
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) Hx[i][j][MAX_Z] = recv_next[i][j];
	if (process_index[2] > 0)
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) Hx[i][j][0] = recv_last[i][j];

	free(recv_next[0]); free(recv_next); free(send_next[0]); free(send_next);
	free(recv_last[0]); free(recv_last); free(send_last[0]); free(send_last);
}

// 3. Hy 在 X 方向交换
void MPI_Hy_x_direction_Change(MPI_Struct* MPI_parameter) {
	int rankid = MPI_parameter->rank_id;
	int p_idx[3] = { MPI_parameter->process_index[0], MPI_parameter->process_index[1], MPI_parameter->process_index[2] };

	int rows = MAX_Y + 1; // Hy 的 Y 维
	int cols = MAX_Z;

	double** send_next = createContinuous2DArray(rows, cols, 0.0);
	double** recv_next = createContinuous2DArray(rows, cols, 0.0);
	double** send_last = createContinuous2DArray(rows, cols, 0.0);
	double** recv_last = createContinuous2DArray(rows, cols, 0.0);
	int num = rows * cols;

	MPI_Request reqs[4]{ MPI_REQUEST_NULL }; MPI_Status stats[4]; int c = 0;

	if (p_idx[0] < MPI_x_direction - 1) {
		int next_rank = MPI_parameter->rank_spatial_distribution[p_idx[0] + 1][p_idx[1]][p_idx[2]];
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) send_next[j][k] = Hy[MAX_X - 1][j][k];
		MPI_Isend(send_next[0], num, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_next[0], num, MPI_DOUBLE, next_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (p_idx[0] > 0) {
		int last_rank = MPI_parameter->rank_spatial_distribution[p_idx[0] - 1][p_idx[1]][p_idx[2]];
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) send_last[j][k] = Hy[1][j][k];
		MPI_Isend(send_last[0], num, MPI_DOUBLE, last_rank, last_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_last[0], num, MPI_DOUBLE, last_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (c > 0) MPI_Waitall(c, reqs, stats);

	if (p_idx[0] < MPI_x_direction - 1)
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) Hy[MAX_X][j][k] = recv_next[j][k];
	if (p_idx[0] > 0)
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) Hy[0][j][k] = recv_last[j][k];

	free(recv_next[0]); free(recv_next); free(send_next[0]); free(send_next);
	free(recv_last[0]); free(recv_last); free(send_last[0]); free(send_last);
}

// 4. Hy 在 Z 方向交换
void MPI_Hy_z_direction_Change(MPI_Struct* MPI_parameter) {
	int rankid = MPI_parameter->rank_id;
	int p_idx[3] = { MPI_parameter->process_index[0], MPI_parameter->process_index[1], MPI_parameter->process_index[2] };

	int rows = MAX_X; // Hy 的 X 维 (需确认 Basic_parameters.h 是否扩充)
	int cols = MAX_Y + 1;

	double** send_next = createContinuous2DArray(rows, cols, 0.0);
	double** recv_next = createContinuous2DArray(rows, cols, 0.0);
	double** send_last = createContinuous2DArray(rows, cols, 0.0);
	double** recv_last = createContinuous2DArray(rows, cols, 0.0);
	int num = rows * cols;

	MPI_Request reqs[4]{ MPI_REQUEST_NULL }; MPI_Status stats[4]; int c = 0;

	if (p_idx[2] < MPI_z_direction - 1) {
		int next_rank = MPI_parameter->rank_spatial_distribution[p_idx[0]][p_idx[1]][p_idx[2] + 1];
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) send_next[i][j] = Hy[i][j][MAX_Z - 1];
		MPI_Isend(send_next[0], num, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_next[0], num, MPI_DOUBLE, next_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (p_idx[2] > 0) {
		int last_rank = MPI_parameter->rank_spatial_distribution[p_idx[0]][p_idx[1]][p_idx[2] - 1];
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) send_last[i][j] = Hy[i][j][1];
		MPI_Isend(send_last[0], num, MPI_DOUBLE, last_rank, last_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_last[0], num, MPI_DOUBLE, last_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (c > 0) MPI_Waitall(c, reqs, stats);

	if (p_idx[2] < MPI_z_direction - 1)
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) Hy[i][j][MAX_Z] = recv_next[i][j];
	if (p_idx[2] > 0)
		for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) Hy[i][j][0] = recv_last[i][j];

	free(recv_next[0]); free(recv_next); free(send_next[0]); free(send_next);
	free(recv_last[0]); free(recv_last); free(send_last[0]); free(send_last);
}

// 5. Hz 在 X 方向交换
void MPI_Hz_x_direction_Change(MPI_Struct* MPI_parameter) {
	int rankid = MPI_parameter->rank_id;
	int p_idx[3] = { MPI_parameter->process_index[0], MPI_parameter->process_index[1], MPI_parameter->process_index[2] };

	int rows = MAX_Y; // Hz 的 Y 维
	int cols = MAX_Z + 1; // Hz 的 Z 维

	double** send_next = createContinuous2DArray(rows, cols, 0.0);
	double** recv_next = createContinuous2DArray(rows, cols, 0.0);
	double** send_last = createContinuous2DArray(rows, cols, 0.0);
	double** recv_last = createContinuous2DArray(rows, cols, 0.0);
	int num = rows * cols;

	MPI_Request reqs[4]{ MPI_REQUEST_NULL }; MPI_Status stats[4]; int c = 0;

	if (p_idx[0] < MPI_x_direction - 1) {
		int next_rank = MPI_parameter->rank_spatial_distribution[p_idx[0] + 1][p_idx[1]][p_idx[2]];
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) send_next[j][k] = Hz[MAX_X - 1][j][k];
		MPI_Isend(send_next[0], num, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_next[0], num, MPI_DOUBLE, next_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (p_idx[0] > 0) {
		int last_rank = MPI_parameter->rank_spatial_distribution[p_idx[0] - 1][p_idx[1]][p_idx[2]];
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) send_last[j][k] = Hz[1][j][k];
		MPI_Isend(send_last[0], num, MPI_DOUBLE, last_rank, last_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_last[0], num, MPI_DOUBLE, last_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (c > 0) MPI_Waitall(c, reqs, stats);

	if (p_idx[0] < MPI_x_direction - 1)
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) Hz[MAX_X][j][k] = recv_next[j][k];
	if (p_idx[0] > 0)
		for (int j = 0; j < rows; j++) for (int k = 0; k < cols; k++) Hz[0][j][k] = recv_last[j][k];

	free(recv_next[0]); free(recv_next); free(send_next[0]); free(send_next);
	free(recv_last[0]); free(recv_last); free(send_last[0]); free(send_last);
}

// 6. Hz 在 Y 方向交换
void MPI_Hz_y_direction_Change(MPI_Struct* MPI_parameter) {
	int rankid = MPI_parameter->rank_id;
	int p_idx[3] = { MPI_parameter->process_index[0], MPI_parameter->process_index[1], MPI_parameter->process_index[2] };

	int rows = MAX_X; // Hz 的 X 维
	int cols = MAX_Z + 1;

	double** send_next = createContinuous2DArray(rows, cols, 0.0);
	double** recv_next = createContinuous2DArray(rows, cols, 0.0);
	double** send_last = createContinuous2DArray(rows, cols, 0.0);
	double** recv_last = createContinuous2DArray(rows, cols, 0.0);
	int num = rows * cols;

	MPI_Request reqs[4]{ MPI_REQUEST_NULL }; MPI_Status stats[4]; int c = 0;

	if (p_idx[1] < MPI_y_direction - 1) {
		int next_rank = MPI_parameter->rank_spatial_distribution[p_idx[0]][p_idx[1] + 1][p_idx[2]];
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) send_next[i][k] = Hz[i][MAX_Y - 1][k];
		MPI_Isend(send_next[0], num, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_next[0], num, MPI_DOUBLE, next_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (p_idx[1] > 0) {
		int last_rank = MPI_parameter->rank_spatial_distribution[p_idx[0]][p_idx[1] - 1][p_idx[2]];
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) send_last[i][k] = Hz[i][1][k];
		MPI_Isend(send_last[0], num, MPI_DOUBLE, last_rank, last_rank, MPI_COMM_WORLD, &reqs[c++]);
		MPI_Irecv(recv_last[0], num, MPI_DOUBLE, last_rank, rankid, MPI_COMM_WORLD, &reqs[c++]);
	}
	if (c > 0) MPI_Waitall(c, reqs, stats);

	if (p_idx[1] < MPI_y_direction - 1)
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) Hz[i][MAX_Y][k] = recv_next[i][k];
	if (p_idx[1] > 0)
		for (int i = 0; i < rows; i++) for (int k = 0; k < cols; k++) Hz[i][0][k] = recv_last[i][k];

	free(recv_next[0]); free(recv_next); free(send_next[0]); free(send_next);
	free(recv_last[0]); free(recv_last); free(send_last[0]); free(send_last);
}

// ==========================================
//          Total Control Function
// ==========================================

void MPI_H_Change(MPI_Struct* MPI_parameter) {
	// 依次调用 6 个方向的交换
	MPI_Hx_y_direction_Change(MPI_parameter);
	MPI_Hx_z_direction_Change(MPI_parameter);

	MPI_Hy_x_direction_Change(MPI_parameter);
	MPI_Hy_z_direction_Change(MPI_parameter);

	MPI_Hz_x_direction_Change(MPI_parameter);
	MPI_Hz_y_direction_Change(MPI_parameter);
}



#endif