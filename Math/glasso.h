#pragma once

#include "Identifier.h"
/*
	start type : cold
	no trace
*/
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> glasso(
	const Eigen::MatrixXd& sss,
	const Eigen::MatrixXd& rho,
	double threshold = 1e-4,
	int maxit = 1e4,
	bool approx = false,
	bool penalize_diagonal = true
);

void lasinv1(
	const Eigen::MatrixXd& ss,
	const Eigen::MatrixXd& rho,
	bool approx,
	int start_type,
	double threshold,
	int maxit,
	bool penalize_diagonal,
	Eigen::MatrixXd& ww,
	Eigen::MatrixXd& wwi,
	int& niter,
	double& del
);

void connect(
	const Eigen::MatrixXd& ss,
	const Eigen::MatrixXd& rho,
	int& nc,
	Eigen::MatrixXi& ic,
	Eigen::ArrayXi& ir,
	Eigen::ArrayXi& ie
);

void setup(
	int m,
	int n,
	const Eigen::MatrixXd& ss,
	const Eigen::MatrixXd& rho,
	const Eigen::MatrixXd& ww,
	Eigen::MatrixXd& vv,
	Eigen::ArrayXd& s,
	Eigen::ArrayXd& r
);

void lasso(
	const Eigen::ArrayXd& rho,
	int n,
	const Eigen::MatrixXd& vv,
	Eigen::ArrayXd& s,
	double threshold,
	Eigen::ArrayXd& x,
	Eigen::ArrayXd& z,
	Eigen::ArrayXi& mm
);

void fatmul(
	int it,
	int n,
	const Eigen::MatrixXd& vv,
	const Eigen::ArrayXd& x,
	Eigen::ArrayXd& s,
	Eigen::ArrayXd& z,
	Eigen::ArrayXi& m
);

void inv(
	const Eigen::MatrixXd& ww,
	Eigen::MatrixXd& xs,
	Eigen::MatrixXd& wwi
);