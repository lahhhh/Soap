#pragma once

#include "Identifier.h"

struct LM_RES {
	Eigen::ArrayXd coefficients;
};

LM_RES lm(
	const Eigen::ArrayXd& y,
	const Eigen::ArrayXd& x,
	bool fit_intercept = true
);