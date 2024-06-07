#pragma once

#include "Identifier.h"


struct GLM_RES {
	Eigen::ArrayXd coefficients;
	Eigen::ArrayXd residuals;
	Eigen::ArrayXd fitted_values;
	Eigen::ArrayXd effects;
	int rank{ 0 };
	struct {
		Eigen::MatrixXd qr;
		int rank;
		Eigen::ArrayXd qraux;
		Eigen::ArrayXi pivot;
		double tol;
	}qr;
	Eigen::ArrayXd linear_predictors;
	double deviance{ 0.0 };
	double null_deviance{ 0.0 };
	int iter{ 0 };
	int df_residual{ 0 };
	int df_null{ 0 };
	Eigen::ArrayXd y;
	bool converged{ false };
	bool boundary{ false };
};

struct GLM_SUMMARY_RES {

	bool hasnan{ false };
	Eigen::ArrayXd estimate;
	Eigen::ArrayXd std_err;
	Eigen::ArrayXd t_value;
	Eigen::ArrayXd p;
};

/*
	note : x now does not have intercept column, we will add it at 1st column later
*/
GLM_RES glm_gaussian(
	const Eigen::ArrayXd& y,
	const Eigen::MatrixXd& x,
	int maxit = 25
);

GLM_SUMMARY_RES glm_summary(
	const GLM_RES& object
);


