#pragma once
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H
#include <stdlib.h>
#include <stdio.h>

// º¯ÊýÉùÃ÷
double*** create3DArray(int x, int y, int z, double num);
double** create2DArray(int rows, int cols, double num);
double* create1DArray(int size, double num);
double** createContinuous2DArray(int rows, int cols, double num);
void free3DArray(double*** array, int x, int y);
void free2DArray(double** array, int rows);
void free1DArray(double* array);

#endif