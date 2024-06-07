#pragma once

#include "Identifier.h"

/*
	functions here are not the same as lapack.
*/

Eigen::MatrixXd chol2inv(const Eigen::MatrixXd& mat);

void dpotri_u(int& n,
	Eigen::MatrixXd& a,
	int& lda
);

void dtrtri_u_nu(
	int& n,
	Eigen::MatrixXd& a,
	int& lda
);

void dtrti2_u_nu(
	int& n,
	Eigen::MatrixXd& a,
	int& lda
);

void dtrmv_u_not_nu(
	int& n,
	Eigen::MatrixXd& a,
	int& lda,
	int jj
);

void dlauu2_u(
	int& n,
	Eigen::MatrixXd& a,
	int& lda
);

void dgemv_not(
	int m,
	int n,
	Eigen::MatrixXd& a,
	double beta,
	int ii
);