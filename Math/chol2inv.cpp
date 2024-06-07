#include "chol2inv.h"

/*
	mat should be positive-definite matrix
*/
Eigen::MatrixXd chol2inv(const Eigen::MatrixXd& mat) {
	
	Eigen::MatrixXd ans = mat;

	int size = mat.rows();

	for (int j = 0; j < size; ++j) {
		for (int i = j + 1; i < size; ++i) {
			ans(i, j) = 0.0;
		}
	}

	dpotri_u(size, ans, size);

	for (int j = 0; j < size; ++j) {
		for (int i = j + 1; i < size; ++i) {
			ans(i, j) = ans(j, i);
		}
	}

	return ans;
}

void dpotri_u(int& n,
	Eigen::MatrixXd& a,
	int& lda) {
	dtrtri_u_nu(
		n, a, lda
	);

	dlauu2_u(n, a, lda);
}

void dtrtri_u_nu(
	int& n,
	Eigen::MatrixXd& a,
	int& lda
) {
	//skip ilaenv()

	dtrti2_u_nu(
		n, a, lda
	);
}

void dtrti2_u_nu(
	int& n,
	Eigen::MatrixXd& a,
	int& lda
) {
	for (int j = 0; j < n; ++j) {

		a(j, j) = 1.0 / a(j, j);
		double ajj = -a(j, j);

		dtrmv_u_not_nu(
			j,
			a,
			lda,
			j
		);
		a.col(j).segment(0, j) *= ajj;
	}
}

void dtrmv_u_not_nu(
	int& n,
	Eigen::MatrixXd& a,
	int& lda,
	int jj
) {
	for (int j = 0; j < n; ++j) {
		if (a(j, jj) != 0.0) {
			double temp = a(j, jj);
			for (int i = 0; i < j; ++i) {
				a(i, jj) += temp * a(i, j);
			}
			a(j, jj) *= a(j, j);
		}
	}
}

void dlauu2_u(
	int& n,
	Eigen::MatrixXd& a,
	int& lda
) {
	for (int i = 0; i < n; ++i) {
		double aii = a(i, i);
		if (i < n - 1) {
			a(i, i) = (a.row(i).segment(i, n - i)).squaredNorm();
			dgemv_not(
				i, n - i - 1, a, aii, i
			);
		}
		else {
			a.col(i).segment(0, i + 1) *= aii;
		}
	}

}


/*
	incy == 1, alpha = 1.0
*/

void dgemv_not(
	int m,
	int n,
	Eigen::MatrixXd& a,
	double beta,
	int ii
) {
	if (beta != 1.0) {
		if (beta == 0.0) {
			for (int i = 0; i < m; ++i) {
				a(i, ii) = 0.0;
			}
		}
		else {
			for (int i = 0; i < m; ++i) {
				a(i, ii) *= beta;
			}
		}
	}

	for (int j = 0; j < n; ++j) {
		double temp = a(ii, ii + j + 1);
		for (int i = 0; i < m; ++i) {
			a(i, ii) += temp * a(i, ii + j + 1);
		}
	}
}