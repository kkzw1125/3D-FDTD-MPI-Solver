#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "engine.h"

class MatlabVisualizer {
private:
    Engine* ep;
    bool isEngineOpen;

public:
    MatlabVisualizer() : ep(nullptr), isEngineOpen(false) {
        // 启动MATLAB引擎
        if (!(ep = engOpen(nullptr))) {
            std::cerr << "错误: 无法启动MATLAB引擎!" << std::endl;
            std::cerr << "请确保MATLAB已安装且路径正确配置" << std::endl;
            return;
        }
        isEngineOpen = true;
        std::cout << "MATLAB引擎启动成功" << std::endl;

        // 设置MATLAB工作目录为当前目录
        engEvalString(ep, "cd(pwd);");

    }

    ~MatlabVisualizer() {
        close();
    }

    bool isOpen() const {
        return isEngineOpen;
    }

    void close() {
        if (ep && isEngineOpen) {
            engClose(ep);
            ep = nullptr;
            isEngineOpen = false;
            std::cout << "MATLAB引擎已关闭" << std::endl;
        }
    }

    // 显示二维场数据
    void showField2D(const std::vector<std::vector<double>>& fieldData,
        const std::string& title = "Field Visualization",
        const std::string& colormap = "jet") {
        if (!isEngineOpen) {
            std::cerr << "MATLAB引擎未启动" << std::endl;
            return;
        }

        int rows = fieldData.size();
        if (rows == 0) return;
        int cols = fieldData[0].size();

        // 创建MATLAB矩阵
        mxArray* mat = mxCreateDoubleMatrix(rows, cols, mxREAL);
        double* data = mxGetPr(mat);

        // 填充数据 (注意MATLAB是列优先)
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                data[j * rows + i] = fieldData[i][j];
            }
        }

        // 发送到MATLAB工作区
        engPutVariable(ep, "field_data", mat);

        // 构建MATLAB命令
        std::string command;
        command += "figure('Name','" + title + "','NumberTitle','off'); ";
        command += "mesh(field_data); ";
        //command += "colormap(" + colormap + "); ";
        //command += "colorbar; ";
        command += "title('" + title + "'); ";
        //command += "axis equal; ";
        command += "axis([1 80 1 100 -0.001 0.001]); ";
        //command += "xlabel('X'); ylabel('Y'); ";
        command += "drawnow; ";  // 强制立即更新图形

        // 执行命令
        engEvalString(ep, command.c_str());

        // 清理内存
        mxDestroyArray(mat);
    }

    // 显示3D场数据的切片
    void showFieldSlice(const std::vector<std::vector<std::vector<double>>>& field3D,
        int slice_dim, int slice_index,
        const std::string& title = "3D Field Slice") {
        if (!isEngineOpen) return;

        std::vector<std::vector<double>> sliceData;

        if (slice_dim == 0) { // X切片
            int x_size = field3D.size();
            if (x_size == 0 || slice_index >= x_size) return;
            int y_size = field3D[0].size();
            int z_size = field3D[0][0].size();

            sliceData.resize(y_size, std::vector<double>(z_size));
            for (int j = 0; j < y_size; j++) {
                for (int k = 0; k < z_size; k++) {
                    sliceData[j][k] = field3D[slice_index][j][k];
                }
            }
        }
        else if (slice_dim == 1) { // Y切片
            int x_size = field3D.size();
            if (x_size == 0) return;
            int y_size = field3D[0].size();
            if (slice_index >= y_size) return;
            int z_size = field3D[0][0].size();

            sliceData.resize(x_size, std::vector<double>(z_size));
            for (int i = 0; i < x_size; i++) {
                for (int k = 0; k < z_size; k++) {
                    sliceData[i][k] = field3D[i][slice_index][k];
                }
            }
        }
        else if (slice_dim == 2) { // Z切片 (通常使用这个)
            int x_size = field3D.size();
            if (x_size == 0) return;
            int y_size = field3D[0].size();
            if (y_size == 0) return;
            int z_size = field3D[0][0].size();
            if (slice_index >= z_size) return;

            sliceData.resize(x_size, std::vector<double>(y_size));
            for (int i = 0; i < x_size; i++) {
                for (int j = 0; j < y_size; j++) {
                    sliceData[i][j] = field3D[i][j][slice_index];
                }
            }
        }

        showField2D(sliceData, title);
    }


    // 执行任意MATLAB命令
    void executeCommand(const std::string& command) {
        if (isEngineOpen) {
            engEvalString(ep, command.c_str());
        }
    }

    // 阻塞程序，直到用户手动关闭当前 Figure 窗口
    void waitForFigureClose() {
        if (!isEngineOpen) return;
        std::cout << "请在 MATLAB 窗口中手动关闭图形以结束程序..." << std::endl;

        // uiwait(gcf) 会阻塞 MATLAB 引擎，直到当前图形窗口被关闭
        // 由于 engEvalString 是同步的，C++ 也会在这里等待
        engEvalString(ep, "uiwait(gcf);");

        std::cout << "图形窗口已关闭，程序继续退出。" << std::endl;
    }

};
