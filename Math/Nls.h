#pragma once

#include "Identifier.h"

double derive(
	double (*func)(double, const Eigen::ArrayXd&),
	double input,
	Eigen::ArrayXd params,
	int n
);

void nls_LevenbergMarquardt(
	double (*func)(double, const Eigen::ArrayXd&),
	const Eigen::ArrayXd& inputs,
	const Eigen::ArrayXd& outputs,
	Eigen::ArrayXd& params
);

