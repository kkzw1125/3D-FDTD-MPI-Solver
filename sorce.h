#pragma once
#ifndef SORCE_H
#define SORCE_H

#include"Basic_functions.h"
#include "Basic_parameters.h"
#include"math.h"

double* source;
double tau = 1.0618e-9;
double t0 = 1 * tau;

void calculate_source() {
	source = create1DArray(timesteps, 0);
	for (int t = 0; t < timesteps; t++) {
		source[t] = exp(-4 * PI * pow(t * dt - t0, 2) / pow(tau, 2));
	}
}


#endif