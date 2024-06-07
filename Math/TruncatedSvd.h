#pragma once
#include "Identifier.h"

// modified from R package irlba
// input [m * n]
// matrix u : [m * nu], s [nu] , v [ n * nu]
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>
tsvd(
	const Eigen::MatrixXd* mat,
	int nu = 50,
	unsigned int random_state = 0,
	double tol = 0.0001,
	int maxit = 100
);

// modified from R package irlba
// input [m * n]
// matrix u : [m * nu], s [nu] , v [ n * nu]
std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd>
tsvd(
	const Eigen::SparseMatrix<double>* mat,
	int nu = 50,
	unsigned int random_state = 0,
	double tol = 0.0001,
	int maxit = 100
);

