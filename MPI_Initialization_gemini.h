#pragma once
#ifndef MPI_INITIAIZATION_H
#define MPI_INITIAIZATION_H

#include "MPI.h"
#include "Basic_parameters.h"
#include "Basic_functions.h"
#include <limits.h>

// 所有的全局网格数
int total_MAX_X, total_MAX_Y, total_MAX_Z;

// 各个方向分配的进程数 (改为变量，不再是 const)
extern int MPI_x_direction = 2;
extern int MPI_y_direction = 1;
extern int MPI_z_direction = 1;

// 结构体定义
typedef struct {
    int x_start;
    int x_end;
    int y_start;
    int y_end;
    int z_start;
    int z_end;
} AxisParameters;

typedef struct {
    int rank_id; // 当前进程ID
    AxisParameters Ex;
    AxisParameters Ey;
    AxisParameters Ez;
    AxisParameters Hx;
    AxisParameters Hy;
    AxisParameters Hz;
    // 修改为三级指针以支持动态分配，因为维度在运行时确定
    int*** rank_spatial_distribution;
    int process_index[3]; // 当前进程的三维坐标 [px, py, pz]
} MPI_Struct;

// --- 论文 5.1.1 自适应区域划分算法 ---
// 寻找最优的 Px, Py, Pz，使得子区域表面积最小
void get_optimal_partition(int num_procs, int Nx, int Ny, int Nz, int* px, int* py, int* pz) {
    double min_surface_area = 1e30; // 初始化为极大值
    int best_px = 1, best_py = 1, best_pz = 1;

    // 遍历所有可能的划分组合
    for (int i = 1; i <= num_procs; i++) {
        if (num_procs % i == 0) { // i 是 Px 的因子
            int remainder_xy = num_procs / i;
            for (int j = 1; j <= remainder_xy; j++) {
                if (remainder_xy % j == 0) { // j 是 Py 的因子
                    int k = remainder_xy / j; // k 是 Pz

                    // 当前划分下的子区域平均尺寸
                    double lx = (double)Nx / i;
                    double ly = (double)Ny / j;
                    double lz = (double)Nz / k;

                    // 计算子区域表面积 S = 2 * (xy + xz + yz)
                    // 通信代价与表面积成正比
                    double current_area = 2.0 * (lx * ly + lx * lz + ly * lz);

                    if (current_area < min_surface_area) {
                        min_surface_area = current_area;
                        best_px = i;
                        best_py = j;
                        best_pz = k;
                    }
                }
            }
        }
    }

    *px = best_px;
    *py = best_py;
    *pz = best_pz;

    // 如果是Rank 0，打印划分结果
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf("[自适应划分] 总进程数: %d\n", num_procs);
        printf("[自适应划分] 最优拓扑 (Px, Py, Pz): %d x %d x %d\n", best_px, best_py, best_pz);
    }
}

// --- 论文 5.1.1 均匀划分策略 ---
// 计算局部网格范围，处理余数分布
// 参数: global_n(全局网格数), num_p(该方向分割数), coord(当前进程在该方向的坐标)
// 输出: start(起始索引), end(结束索引, 不包含), local_n(局部大小)
void calculate_local_range(int global_n, int num_p, int coord, int* start, int* end, int* local_n) {
    int base = global_n / num_p;      // 基准大小
    int remainder = global_n % num_p; // 余数

    if (coord < remainder) {
        *local_n = base + 1;
        *start = coord * (base + 1);
    }
    else {
        *local_n = base;
        *start = remainder * (base + 1) + (coord - remainder) * base;
    }
    *end = *start + *local_n;
}

void MPI_INITIAIZATION(int rank, MPI_Struct* MPI_parameter) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_parameter->rank_id = rank;

    // 1. 保存全局网格大小 (从 Basic_parameters.h 获取)
    total_MAX_X = MAX_X;
    total_MAX_Y = MAX_Y;
    total_MAX_Z = MAX_Z;

    // 2. 执行自适应区域划分
    get_optimal_partition(size, total_MAX_X, total_MAX_Y, total_MAX_Z,&MPI_x_direction, &MPI_y_direction, &MPI_z_direction);

    // 3. 动态分配 rank_spatial_distribution 数组
    MPI_parameter->rank_spatial_distribution = (int***)malloc(MPI_x_direction * sizeof(int**));
    for (int i = 0; i < MPI_x_direction; i++) {
        MPI_parameter->rank_spatial_distribution[i] = (int**)malloc(MPI_y_direction * sizeof(int*));
        for (int j = 0; j < MPI_y_direction; j++) {
            MPI_parameter->rank_spatial_distribution[i][j] = (int*)malloc(MPI_z_direction * sizeof(int));
        }
    }

    // 4. 建立笛卡尔虚拟拓扑映射
    int process_index[3] = { 0 };
    for (int i = 0; i < MPI_x_direction; i++) {
        for (int j = 0; j < MPI_y_direction; j++) {
            for (int k = 0; k < MPI_z_direction; k++) {
                // 简单的线性映射，也可以使用 MPI_Cart_create
                int current_rank = (i * MPI_y_direction + j) * MPI_z_direction + k;
                MPI_parameter->rank_spatial_distribution[i][j][k] = current_rank;

                if (rank == current_rank) {
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

    // 5. 计算当前进程的局部网格范围（应用余数均匀分配策略）
    // 注意：Basic_parameters.h 定义了 MAX_X 是 grid cells 数量。
    // Field 数组通常比 grid cells 多 1 (例如 Ex 在 x 方向是 MAX_X, 但 Ey 在 x 方向是 MAX_X + 1)
    // 这里我们主要基于网格单元(Cells)进行划分，具体的场分量索引需微调

    int local_nx, local_ny, local_nz;
    int x_s, x_e, y_s, y_e, z_s, z_e;

    calculate_local_range(total_MAX_X, MPI_x_direction, process_index[0], &x_s, &x_e, &local_nx);
    calculate_local_range(total_MAX_Y, MPI_y_direction, process_index[1], &y_s, &y_e, &local_ny);
    calculate_local_range(total_MAX_Z, MPI_z_direction, process_index[2], &z_s, &z_e, &local_nz);

    // 更新全局变量 MAX_X 为当前进程的局部大小 (为了兼容后续计算代码)
    MAX_X = local_nx;
    MAX_Y = local_ny;
    MAX_Z = local_nz;

    // 6. 设置具体的场分量范围 
    // 注意：Ex 位于 (i+0.5, j, k)，在 X 方向上如果不处于最后边界，通常索引数量等于网格数
    // 下面的逻辑基于原代码逻辑进行适配，但在边界处理上更精确

    // Ex: 定义在 X 轴边上
    MPI_parameter->Ex.x_start = x_s;
    MPI_parameter->Ex.x_end = x_e; // 注意原代码逻辑似乎包含边界
    MPI_parameter->Ex.y_start = y_s;
    MPI_parameter->Ex.y_end = y_e;
    MPI_parameter->Ex.z_start = z_s;
    MPI_parameter->Ex.z_end = z_e;
    // 边界修正: Ex 在 Y 和 Z 方向是节点，如果是最后一个区域，可能需要多存一个点
    // 但原代码 MAX_Y 是指网格数，Ex[MAX_X][MAX_Y+1][MAX_Z+1]。
    // 此处保持 start/end 为全局索引，用于源的定位等

    // Ey
    MPI_parameter->Ey.x_start = x_s;
    MPI_parameter->Ey.x_end = x_e;
    MPI_parameter->Ey.y_start = y_s;
    MPI_parameter->Ey.y_end = y_e;
    MPI_parameter->Ey.z_start = z_s;
    MPI_parameter->Ey.z_end = z_e;

    // Ez
    MPI_parameter->Ez.x_start = x_s;
    MPI_parameter->Ez.x_end = x_e;
    MPI_parameter->Ez.y_start = y_s;
    MPI_parameter->Ez.y_end = y_e;
    MPI_parameter->Ez.z_start = z_s;
    MPI_parameter->Ez.z_end = z_e;

    // Hx, Hy, Hz 范围同理，通常磁场位于元胞中心或对偶网格
    MPI_parameter->Hx.x_start = x_s; MPI_parameter->Hx.x_end = x_e;
    MPI_parameter->Hx.y_start = y_s; MPI_parameter->Hx.y_end = y_e;
    MPI_parameter->Hx.z_start = z_s; MPI_parameter->Hx.z_end = z_e;

    MPI_parameter->Hy.x_start = x_s; MPI_parameter->Hy.x_end = x_e;
    MPI_parameter->Hy.y_start = y_s; MPI_parameter->Hy.y_end = y_e;
    MPI_parameter->Hy.z_start = z_s; MPI_parameter->Hy.z_end = z_e;

    MPI_parameter->Hz.x_start = x_s; MPI_parameter->Hz.x_end = x_e;
    MPI_parameter->Hz.y_start = y_s; MPI_parameter->Hz.y_end = y_e;
    MPI_parameter->Hz.z_start = z_s; MPI_parameter->Hz.z_end = z_e;

    // 边界条件修正 (针对原有代码逻辑的兼容)
    // 如果是该方向的最后一个进程，确保覆盖到物理边界
    if (process_index[0] == MPI_x_direction - 1) {
        MPI_parameter->Ex.x_end = total_MAX_X - 1; // Ex 在 x 方向通常只有 Nx 个
        MPI_parameter->Ey.x_end = total_MAX_X;     // Ey 在 x 方向边缘有值
        MPI_parameter->Ez.x_end = total_MAX_X;
    }
    if (process_index[1] == MPI_y_direction - 1) {
        MPI_parameter->Ex.y_end = total_MAX_Y;
        MPI_parameter->Ey.y_end = total_MAX_Y - 1;
        MPI_parameter->Ez.y_end = total_MAX_Y;
    }
    if (process_index[2] == MPI_z_direction - 1) {
        MPI_parameter->Ex.z_end = total_MAX_Z;
        MPI_parameter->Ey.z_end = total_MAX_Z;
        MPI_parameter->Ez.z_end = total_MAX_Z - 1;
    }
}

// 辅助函数：释放 MPI_Struct 中动态分配的内存 (建议在 main 函数结束前调用)
void FREE_MPI_STRUCT(MPI_Struct* MPI_parameter) {
    if (MPI_parameter->rank_spatial_distribution != NULL) {
        for (int i = 0; i < MPI_x_direction; i++) {
            for (int j = 0; j < MPI_y_direction; j++) {
                free(MPI_parameter->rank_spatial_distribution[i][j]);
            }
            free(MPI_parameter->rank_spatial_distribution[i]);
        }
        free(MPI_parameter->rank_spatial_distribution);
    }
}

#endif#pragma once
