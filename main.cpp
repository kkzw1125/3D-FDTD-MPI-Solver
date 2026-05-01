#include <stdio.h>
#include"FDTD_headfiles.h"
//#include"MPI_Initialization.h"
#include"MPI_Initialization_gemini.h"
//#include"mpi.h"
#include "FDTD_MPI_Functions.h"
#include <time.h>

#include <iostream>
//#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include "MatlabVisualizer2.h" 
#include "CPML.h" 

using namespace std;




int main(int argc, char** argv) {

    // 初始化MATLAB可视化
   
    MatlabVisualizer matlabVis;
    if (!matlabVis.isOpen()) {
        cout << "警告: MATLAB可视化不可用，将继续运行无可视化" << endl;
    }
    else {
        cout << "MATLAB可视化已启用" << endl;
    }
    
    
    int rank, size;
    MPI_Struct MPI_parameter;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_INITIAIZATION(rank, &MPI_parameter);
    //MPI_INITIAIZATION(rank,size,&MPI_parameter);//doubao
    printf("now rank:%d\n", rank);
    //printf("MPI_parameter_rank: %d\t process[%d][%d][%d]\n",\
    //    MPI_parameter.rank_id, MPI_parameter.process_index[0], MPI_parameter.process_index[1], MPI_parameter.process_index[2]);
    //printf("totalMAX:%d,%d,%d\n", total_MAX_X, total_MAX_Y, total_MAX_Z);
    printf("rank_MAX:%d,%d,%d\n", MAX_X, MAX_Y, MAX_Z);

    //printf("x_start:%d,%d,%d\n", MPI_parameter.Ex.x_start, MPI_parameter.Ex.y_start, MPI_parameter.Ex.z_start);
    //printf("x_end:%d,%d,%d\n", MPI_parameter.Ex.x_end, MPI_parameter.Ex.y_end, MPI_parameter.Ex.z_end);

    Initialization_grid();  
    calculate_source();
    MPI_Barrier(MPI_COMM_WORLD);//等待函数
    
    // ==========================================
    // [修复代码]：根据 MPI 拓扑关闭内部边界的 PML
    // ==========================================
    int px = MPI_parameter.process_index[0];
    int py = MPI_parameter.process_index[1];
    int pz = MPI_parameter.process_index[2];

    // X 方向处理
    if (MPI_x_direction > 1) {
        // 如果左边有邻居，关闭左侧 PML
        if (px > 0) nxPML1 = 0;
        // 如果右边有邻居，关闭右侧 PML
        if (px < MPI_x_direction - 1) nxPML2 = 0;
    }
    // Y 方向处理
    if (MPI_y_direction > 1) {
        if (py > 0) nyPML1 = 0;
        if (py < MPI_y_direction - 1) nyPML2 = 0;
    }
    // Z 方向处理
    if (MPI_z_direction > 1) {
        if (pz > 0) nzPML1 = 0;
        if (pz < MPI_z_direction - 1) nzPML2 = 0;
    }
 
    printf("Rank %d PML Settings: X[%d, %d], Y[%d, %d], Z[%d, %d]\n",
        rank, nxPML1, nxPML2, nyPML1, nyPML2, nzPML1, nzPML2);

    // ==========================================
    // [修复结束]
    // ==========================================
   
    

    // 初始化并计算CPML
    init_cpml_variables();
    calculate_cpml();
   /*
    double ez_low_m1 = 0.0, ez_low_m2 = 0.0;
    double ez_high_m1 = 0.0, ez_high_m2 = 0.0;
   */


    clock_t start_time = clock();
    
    for (int t = 0; t < timesteps; t++) {
        /*
        // 更新入射波
        for (int j = 1; j < MAX_Y; j++) {
            ez_inc[j] = ez_inc[j] - CBy[1][j][1] * (hx_inc[j] - hx_inc[j - 1]) / dz;
        }

        ez_inc[0] = ez_low_m2;
        ez_low_m2 = ez_low_m1;
        ez_low_m1 = ez_inc[1];

        ez_inc[MAX_Y - 1] = ez_high_m2;
        ez_high_m2 = ez_high_m1;
        ez_high_m1 = ez_inc[MAX_Y - 2];

        double Tw = 2e-10, t0 = 4 * Tw;
        double pulse = exp(-pow((t * dt - t0) / Tw, 2));
        ez_inc[7] = ez_inc[7] + pulse;
*/



        Calculate_electric_filed();
        
        
         // =======================================================
         // [修复代码]：自动判断哪个进程应该添加源，并转换坐标
         // =======================================================

         // 1. 定义源的全局坐标
        int source_global_x = 70;
        int source_global_y = 35;
        int source_global_z = 20;

        // 2. 检查当前进程是否包含该坐标
        // 必须同时检查 X, Y, Z (虽然目前主要是 X 分割，但这样写最稳健)
        bool in_x_range = (source_global_x >= MPI_parameter.Ez.x_start && source_global_x < MPI_parameter.Ez.x_end);
        bool in_y_range = (source_global_y >= MPI_parameter.Ez.y_start && source_global_y < MPI_parameter.Ez.y_end);
        bool in_z_range = (source_global_z >= MPI_parameter.Ez.z_start && source_global_z < MPI_parameter.Ez.z_end);

        if (in_x_range && in_y_range && in_z_range) {
            // 3. 将全局坐标转换为本地坐标
            // 本地坐标 = 全局坐标 - 当前进程的起始偏移量
            int local_i = source_global_x - MPI_parameter.Ez.x_start;
            int local_j = source_global_y - MPI_parameter.Ez.y_start;
            int local_k = source_global_z - MPI_parameter.Ez.z_start;

            // 4. 添加源
            Ez[local_i][local_j][local_k] = source[t];

            // 调试打印 
            if (t == 0) printf("Rank %d added source at Local[%d][%d][%d]\n", rank, local_i, local_j, local_k);
        }
        // =======================================================
       

        MPI_E_Change(&MPI_parameter);
       
       /* // 更新入射磁场
        for (int j = 0; j < MAX_Y - 1; j++) {
            hx_inc[j] = hx_inc[j] + CQy[1][1][1] * (ez_inc[j] - ez_inc[j + 1]) / dx;
        }
*/
        Calculate_magnetic_filed();

        MPI_H_Change(&MPI_parameter);

        MPI_Barrier(MPI_COMM_WORLD);
        

        // --- 可视化部分开始 ---
        int vis_step = 25;
        int target_global_z = 20;// 全局 Z 切面这一层

        if ((t + 1) % vis_step == 0) {

            // 1. 准备全局切片数据的容器 (仅 Rank 0 需要分配内存)
            // 使用 vector 方便传给你的 Visualization 类
            // 注意：total_MAX_X 需要在 main 开头从 MPI_Initialization 获取或定义为全局
            std::vector<std::vector<double>> GlobalSlice;

            if (rank == 0) {
                cout << "时间步 " << t + 1 << " - 正在收集数据并绘图..." << endl;
                GlobalSlice.resize(total_MAX_X); // 全局大小
                for (int i = 0; i < total_MAX_X; ++i) GlobalSlice[i].resize(total_MAX_Y);
            }

            // 2. 临时存储当前进程的切片数据
            // 判断当前进程是否包含目标 Z 层
            bool i_have_z_slice = (target_global_z >= MPI_parameter.Ez.z_start && target_global_z < MPI_parameter.Ez.z_end);
            int local_z_idx = target_global_z - MPI_parameter.Ez.z_start;

            // 3. 数据收集逻辑
            if (rank == 0) {
                // A. 先填入 Rank 0 自己的数据 (如果有的话)
                if (i_have_z_slice) {
                    for (int i = 0; i < MAX_X; i++) {
                        for (int j = 0; j < MAX_Y; j++) {
                            // 算出全局坐标，填入 GlobalSlice
                            int g_i = i + MPI_parameter.Ez.x_start;
                            int g_j = j + MPI_parameter.Ez.y_start;
                            if (g_i < total_MAX_X && g_j < total_MAX_Y) // 安全检查
                                GlobalSlice[g_i][g_j] = Ez[i][j][local_z_idx];
                        }
                    }
                }

                // B. 接收其他 Rank 的数据
                for (int r = 1; r < size; r++) {
                    // 这里为了简单，我们需要知道发送方的数据大小和偏移量。
                    // 严谨的做法应该先 Recv 对方的坐标范围，再 Recv 数据。
                    // 简易做法：假设所有从属进程依次发送它们的那一块切片

                    // 接收头部信息：[x_start, x_len, y_start, y_len, has_data]
                    int header[5];
                    MPI_Status status;
                    MPI_Recv(header, 5, MPI_INT, r, 999, MPI_COMM_WORLD, &status);

                    if (header[4] == 1) { // 对方有数据
                        int r_xs = header[0]; int r_xl = header[1];
                        int r_ys = header[2]; int r_yl = header[3];

                        // 接收实际数据 buffer
                        std::vector<double> recv_buf(r_xl * r_yl);
                        MPI_Recv(recv_buf.data(), r_xl * r_yl, MPI_DOUBLE, r, 888, MPI_COMM_WORLD, &status);

                        // 填入 GlobalSlice
                        int count = 0;
                        for (int i = 0; i < r_xl; i++) {
                            for (int j = 0; j < r_yl; j++) {
                                GlobalSlice[r_xs + i][r_ys + j] = recv_buf[count++];
                            }
                        }
                    }
                }

                // 4. Rank 0 调用 Matlab 显示
                matlabVis.showField2D(GlobalSlice, "Ez Field (Global) at T=" + to_string(t + 1));

            }
            else {
                // --- 从属进程 (Rank > 0) ---
                int header[5] = { 0 };
                std::vector<double> send_buf;

                if (i_have_z_slice) {
                    header[0] = MPI_parameter.Ez.x_start;
                    header[1] = MAX_X; // 本地 X 长度
                    header[2] = MPI_parameter.Ez.y_start;
                    header[3] = MAX_Y; // 本地 Y 长度
                    header[4] = 1;     // 标记有数据

                    // 压平数据准备发送
                    send_buf.resize(MAX_X * MAX_Y);
                    int count = 0;
                    for (int i = 0; i < MAX_X; i++) {
                        for (int j = 0; j < MAX_Y; j++) {
                            send_buf[count++] = Ez[i][j][local_z_idx];
                        }
                    }
                }
                else {
                    header[4] = 0; // 标记无数据 
                }

                // 发送头部
                MPI_Send(header, 5, MPI_INT, 0, 999, MPI_COMM_WORLD);
                // 发送数据 (如果有)
                if (header[4] == 1) {
                    MPI_Send(send_buf.data(), send_buf.size(), MPI_DOUBLE, 0, 888, MPI_COMM_WORLD);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // 绘图时同步一下，防止跑太快
        // --- 可视化部分结束 ---
           
       

    }

    clock_t end_time = clock();
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // 输出仿真时间
    printf("程序执行时间: %f 秒\n", time_taken);


    free_array();    
    free_cpml_variables();
  
    MPI_Finalize();

    FREE_MPI_STRUCT(&MPI_parameter);

    
    return 0;
}


