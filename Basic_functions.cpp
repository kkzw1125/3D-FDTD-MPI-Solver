#include "Basic_functions.h"

double*** create3DArray(int x, int y, int z, double num) {
    double*** array = (double***)malloc(x * sizeof(double**));
    if (array == NULL) {
        fprintf(stderr, "create3DArray:Memory allocation failed for array.\n");
        exit(EXIT_FAILURE);
    }
    else {
        for (int i = 0; i < x; i++) {
            array[i] = (double**)malloc(y * sizeof(double*));
            for (int j = 0; j < y; j++) {
                array[i][j] = (double*)malloc(z * sizeof(double));
                // 初始化为指定的值
                for (int k = 0; k < z; k++) {
                    array[i][j][k] = num;
                }
            }
        }
        return array;
    }
}

double** create2DArray(int rows, int cols, double num) {
    double** array = (double**)malloc(rows * sizeof(double*));
    if (array == NULL) {
        fprintf(stderr, "create2DArray:Memory allocation failed for array.\n");
        exit(EXIT_FAILURE);
    }
    else {
        for (int i = 0; i < rows; i++) {
            array[i] = (double*)malloc(cols * sizeof(double));
            // 初始化为指定的值
            for (int j = 0; j < cols; j++) {
                array[i][j] = num;
            }
        }
        return array;
    }
}


double* create1DArray(int size, double num) {
    // 【关键修复】
    // 必须改为 <= 0。
    // 因为当 nxPML=0 时，CPML.cpp 会请求 (nxPML - 1) 即 -1 的大小。
    // 此时必须返回 NULL，否则 malloc 会崩溃。
    if (size <= 0) {
        return NULL;
    }

    double* array = (double*)malloc(size * sizeof(double));

    if (array == NULL) {
        fprintf(stderr, "create1DArray:Memory allocation failed for array.\n");
        // 打印出请求的大小，方便调试（可选）
        fprintf(stderr, "Requested size: %d\n", size);
        exit(EXIT_FAILURE);
    }
    else {
        for (int i = 0; i < size; i++) {
            array[i] = num;
        }
        return array;
    }
}
/*double* create1DArray(int size, double num) {
    double* array = (double*)malloc(size * sizeof(double));
    if (array == NULL) {
        fprintf(stderr, "create1DArray:Memory allocation failed for array.\n");
        exit(EXIT_FAILURE);
    }
    else {
        for (int i = 0; i < size; i++) {
            array[i] = num;
        }
        return array;
    }
}
*/
double** createContinuous2DArray(int rows, int cols, double num) {
    // 分配存储行指针的空间
    double** array = (double**)malloc(rows * sizeof(double*));
    if (array == NULL) {
        fprintf(stderr, "createContinuous2DArray: Memory allocation failed for array.\n");
        exit(EXIT_FAILURE);
    }

    // 分配实际连续的内存块
    array[0] = (double*)malloc(rows * cols * sizeof(double));
    if (array[0] == NULL) {
        fprintf(stderr, "createContinuous2DArray: Memory allocation failed for array data.\n");
        free(array); // 释放之前分配的行指针
        exit(EXIT_FAILURE);
    }

    // 设置每一行的指针
    for (int i = 1; i < rows; i++) {
        array[i] = array[0] + i * cols;
    }

    // 初始化数组为指定值
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] = num;
        }
    }

    return array;
}

void free3DArray(double*** array, int x, int y) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

void free2DArray(double** array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

void free1DArray(double* array) {
    free(array);
}