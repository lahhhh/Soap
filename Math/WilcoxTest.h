#pragma once

#include "Identifier.h"

/*
Modified from R wilcox.test
*/

/*
alternative - 0 : two-sides ; 1 : less ; 2 : greater
*/
double wilcox_test(
	const Eigen::ArrayXd& x,
	const Eigen::ArrayXd& y,
	bool paired = false,
	int alternative = 0,
	double mu = 0.0,
	bool correct = true
);