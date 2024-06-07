#pragma once

#include "Identifier.h"


struct DQRLS_RES {
	Eigen::MatrixXd qr;
	Eigen::ArrayXd coefficients;
	Eigen::ArrayXd residuals;
	Eigen::ArrayXd effects;
	int rank;
	Eigen::ArrayXi pivot;
	Eigen::ArrayXd qraux;
	double tol;
	bool pivoted;
};

/*
	specialized for single response
*/
DQRLS_RES cdqrls(
	const Eigen::MatrixXd& x,
	const Eigen::ArrayXd& y,
	double tol
);

/*
	modified from r package dqrdc2.f
*/

/*
	nrow of x >= ncol of x
*/

void dqrls(
	Eigen::MatrixXd& x,
	Eigen::ArrayXd& y,
	double tol,
	Eigen::ArrayXd& b,
	Eigen::ArrayXd& rsd,
	Eigen::ArrayXd& qty,
	int& k,
	Eigen::ArrayXi& jpvt,
	Eigen::ArrayXd& qraux,
	Eigen::MatrixXd& work
);

void dqrdc2(
	Eigen::MatrixXd& x,
	double tol,
	int& k,
	Eigen::ArrayXd& qraux,
	Eigen::ArrayXi& jpvt,
	Eigen::MatrixXd& work
);

void dqrsl(
	Eigen::MatrixXd& x,
	int n,
	int& k,
	Eigen::ArrayXd& qraux,
	Eigen::ArrayXd& y,
	Eigen::ArrayXd& qy,
	Eigen::ArrayXd& qty,
	Eigen::ArrayXd& b,
	Eigen::ArrayXd& rsd,
	Eigen::ArrayXd& xb,
	int job
);
