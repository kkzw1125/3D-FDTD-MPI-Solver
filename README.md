# 3D-FDTD-MPI-Solver

基于 C++ 和 MPI (Message Passing Interface) 实现的高性能三维时域有限差分 (FDTD) 电磁波仿真求解器。本项目采用区域分解法进行并行计算，集成了 CPML (卷积完美匹配层) 吸收边界条件，并通过 MATLAB Engine 实现仿真过程中的电磁场实时动态可视化。

## ✨ 主要特性 (Key Features)

* **三维全波电磁仿真**: 基于 Yee 网格的完整三维 Maxwell 方程组求解 (Ex, Ey, Ez, Hx, Hy, Hz)。
* **MPI 分布式并行加速**: 
  * 支持三维笛卡尔虚拟拓扑的自适应区域划分。
  * 自动处理多进程间 6 个方向边界的场分量交换 (`FDTD_MPI_Functions.h`)。
* **CPML 吸收边界条件**: 
  * 实现了高阶 CPML 截断边界。
  * **拓扑感知优化**: 能够根据当前 MPI 进程的相对空间位置，自动识别并关闭 MPI 子区域内部通信面上的虚拟 PML 边界，避免多余计算。
* **实时动态可视化**: 集成 MATLAB Engine API，在 C++ 迭代计算期间实时采集并汇总全局场切片数据 (如 Z 面切片)，调用 MATLAB 绘制二维场分布动态图 (`MatlabVisualizer2.h`)。
* **高斯脉冲激励源**: 内置软源激励注入机制 (`sorce.h`)，并支持全局坐标到子进程局部坐标的自动映射。

## 🛠️ 环境依赖 (Dependencies)

本项目需要链接外部并行库和可视化引擎，推荐使用 **Visual Studio 2022** 进行项目管理。

* **C++ 编译器**: 支持 C++11 及以上标准。
* **MPI 环境**: MS-MPI (Microsoft MPI) v10.0 或更高版本。
* **MATLAB**: 需安装 MATLAB 并确保其支持 C/C++ Engine API (`engine.h`)。

## 📁 核心项目结构 (Project Structure)

* `main.cpp`: 主程序入口，包含 MPI 初始化、时间步迭代、CPML 计算调用以及 MATLAB 绘图的数据汇总逻辑。
* `Basic_parameters.h / .cpp`: 定义全局物理常数、网格尺寸 (`MAX_X`, `MAX_Y`, `MAX_Z`)、空间步长 (`dx`, `dy`, `dz`) 及时间步长 (`dt`)。
* `MPI_Initialization_gemini.h`: 包含 MPI 自适应区域划分算法与局部网格映射逻辑。
* `FDTD_MPI_Functions.h`: 处理分布在不同进程间的电场 (E) 和磁场 (H) 切面数据的非阻塞发送与接收 (MPI_Isend/MPI_Irecv)。
* `CPML.h / .cpp`: CPML 参数初始化、辅助系数数组计算及内存释放。
* `Calculate_electric/magnetic_filed.h`: FDTD 电磁场核心更新方程。
* `MatlabVisualizer2.h`: 封装的 MATLAB 引擎类，用于建立连接并下发 `figure`, `mesh`, `drawnow` 等绘图指令。

## 🚀 编译与运行 (Build and Run)

### 1. 配置 Visual Studio 项目属性
在编译本项目前，请确保在 IDE 中正确配置了以下路径：
* **包含目录 (Include Directories)**:
  * MS-MPI Include 路径 (例如: `C:\Program Files (x86)\Microsoft SDKs\MPI\Include`)
  * MATLAB Engine Include 路径 (例如: `C:\Program Files\MATLAB\R202x\extern\include`)
* **库目录 (Library Directories)**:
  * MS-MPI Lib 路径 (例如: `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64`)
  * MATLAB Engine Lib 路径 (例如: `C:\Program Files\MATLAB\R202x\extern\lib\win64\microsoft`)
* **附加依赖项 (Additional Dependencies)**:
  * `msmpi.lib`, `libeng.lib`, `libmx.lib`, `libmat.lib`

### 2. 参数修改
若需调整仿真模型，请修改以下文件：
* **仿真空间与步数**：在 `Basic_parameters.cpp` 中修改 `MAX_X`, `MAX_Y`, `MAX_Z` 以及 `timesteps`。
* **并行拓扑设置**：程序默认开启自适应划分，如需强制指定进程拓扑结构，可在 `MPI_Initialization_gemini.h` 中修改 `MPI_x_direction`, `MPI_y_direction`, `MPI_z_direction`。
* **激励源位置**：在 `main.cpp` 的时间步循环中，修改 `source_global_x`, `y`, `z` 坐标。

### 3. 运行程序
编译生成可执行文件 (如 `FDTD_MPI.exe`) 后，在命令行使用 `mpiexec` 启动并行计算。例如，使用 4 个进程运行：

```bash
mpiexec -n 4 FDTD_MPI.exe
