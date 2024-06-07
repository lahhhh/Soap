#include "TruncatedSvd.h"
#include "Custom.h"

inline double invcheck(double s) {
	constexpr double eps2 = std::numeric_limits<double>::epsilon() * 2;
	if (s > eps2) {
		return 1.0 / s;
	}
	else {
		return 0;
	}
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> tsvd(
	const Eigen::MatrixXd* mat,
	int n,
	unsigned int random_state,
	double tol,
	int maxit
) {
	const int nu = n;
	const int m = mat->rows();
	n = mat->cols();
	const int m_b = _Cs min(nu + 20, 3 * nu, n);
	int it = 0;
	int j = 0;
	int k = nu;
	double smax = 1;

	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, m_b);
	Eigen::MatrixXd W = Eigen::MatrixXd::Zero(m, m_b);
	Eigen::VectorXd F = Eigen::VectorXd::Zero(n);
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_b, m_b);

	srand(random_state);
	V.col(0) = Eigen::VectorXd::Random(n);

	V.col(0).array() /= V.norm();

	Eigen::JacobiSVD<Eigen::MatrixXd> svd;

	while (it < maxit) {

		if (it > 0) {
			j = k;
		}

		W.col(j) = (*mat) * V.col(j);

		if (it > 0) {
			W.col(j) = W.col(j) - W.block(0, 0, m, j) * (W.block(0, 0, m, j).transpose() * W.col(j));
		}

		double s = W.col(j).norm();
		double sinv = invcheck(s);

		W.col(j).array() *= sinv;

		double fn, fninv;
		while (j < m_b) {
			F = (W.col(j).transpose() * (*mat)).transpose() - s * V.col(j);
			F = F - V.block(0, 0, n, j + 1) * (V.block(0, 0, n, j + 1).transpose() * F);

			fn = F.norm();
			fninv = invcheck(fn);
			F.array() *= fninv;

			if (j < m_b - 1) {
				V.col(j + 1) = F;
				B(j, j) = s;
				B(j, j + 1) = fn;
				W.col(j + 1) = (*mat) * V.col(j + 1) - fn * W.col(j);
				W.col(j + 1) = W.col(j + 1) - W.block(0, 0, m, j + 1) * (W.block(0, 0, m, j + 1).transpose() * W.col(j + 1));

				s = W.col(j + 1).norm();
				sinv = invcheck(s);
				W.col(j + 1).array() *= sinv;
			}
			else {
				B(j, j) = s;
			}
			++j;
		}

		svd.compute(B, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::VectorXd R = fn * svd.matrixU().row(m_b - 1);

		if (it < 1) {
			smax = svd.singularValues()[0];
		}
		else {
			smax = std::max(smax, svd.singularValues()[0]);
		}

		int conv{ 0 };
		for (int i = 0; i < nu; ++i) {
			if (std::abs(R[i]) < tol * smax) {
				++conv;
			}
		}

		if (conv < nu) {
			k = std::max(conv + nu, k);
			k = std::min(k, m_b - 3);
		}
		else {
			break;
		}

		V.block(0, 0, n, k) = V * (svd.matrixV().transpose().block(0, 0, m_b, k));
		V.col(k) = F;

		B.setZero();

		for (int i = 0; i < k; ++i) {
			B(i, i) = svd.singularValues()[i];
			B(i, k) = R[i];
		}

		W.block(0, 0, m, k) = W * svd.matrixU().block(0, 0, m_b, k);

		++it;
	}

	Eigen::MatrixXd mat_u = W * svd.matrixU().block(0, 0, m_b, nu);
	Eigen::MatrixXd mat_v = V * svd.matrixV().transpose().block(0, 0, m_b, nu);
	Eigen::VectorXd eigen_singular_values = svd.singularValues().segment(0, nu);
	return std::make_tuple(mat_u, eigen_singular_values, mat_v);
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> tsvd(
	const Eigen::SparseMatrix<double>* mat,
	int n,
	unsigned int random_state,
	double tol,
	int maxit
) {
	const int nu = n;
	const int m = mat->rows();
	n = mat->cols();
	const int m_b = _Cs min(nu + 20, 3 * nu, n);
	int it = 0;
	int j = 0;
	int k = nu;
	double smax = 1;

	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, m_b);
	Eigen::MatrixXd W = Eigen::MatrixXd::Zero(m, m_b);
	Eigen::VectorXd F = Eigen::VectorXd::Zero(n);
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_b, m_b);

	srand(random_state);
	V.col(0) = Eigen::VectorXd::Random(n);

	V.col(0).array() /= V.norm();

	Eigen::JacobiSVD<Eigen::MatrixXd> svd;

	while (it < maxit) {

		if (it > 0) {
			j = k;
		}

		W.col(j) = (*mat) * V.col(j);

		if (it > 0) {
			W.col(j) = W.col(j) - W.block(0, 0, m, j) * (W.block(0, 0, m, j).transpose() * W.col(j));
		}

		double s = W.col(j).norm();
		double sinv = invcheck(s);

		W.col(j).array() *= sinv;

		double fn, fninv;
		while (j < m_b) {
			F = (W.col(j).transpose() * (*mat)).transpose() - s * V.col(j);
			F = F - V.block(0, 0, n, j + 1) * (V.block(0, 0, n, j + 1).transpose() * F);

			fn = F.norm();
			fninv = invcheck(fn);
			F.array() *= fninv;

			if (j < m_b - 1) {
				V.col(j + 1) = F;
				B(j, j) = s;
				B(j, j + 1) = fn;
				W.col(j + 1) = (*mat) * V.col(j + 1) - fn * W.col(j);
				W.col(j + 1) = W.col(j + 1) - W.block(0, 0, m, j + 1) * (W.block(0, 0, m, j + 1).transpose() * W.col(j + 1));

				s = W.col(j + 1).norm();
				sinv = invcheck(s);
				W.col(j + 1).array() *= sinv;
			}
			else {
				B(j, j) = s;
			}
			++j;
		}

		svd.compute(B, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::VectorXd R = fn * svd.matrixU().row(m_b - 1);

		if (it < 1) {
			smax = svd.singularValues()[0];
		}
		else {
			smax = std::max(smax, svd.singularValues()[0]);
		}

		int conv{ 0 };
		for (int i = 0; i < nu; ++i) {
			if (std::abs(R[i]) < tol * smax) {
				++conv;
			}
		}

		if (conv < nu) {
			k = std::max(conv + nu, k);
			k = std::min(k, m_b - 3);
		}
		else {
			break;
		}

		V.block(0, 0, n, k) = V * (svd.matrixV().transpose().block(0, 0, m_b, k));
		V.col(k) = F;

		B.setZero();

		for (int i = 0; i < k; ++i) {
			B(i, i) = svd.singularValues()[i];
			B(i, k) = R[i];
		}

		W.block(0, 0, m, k) = W * svd.matrixU().block(0, 0, m_b, k);

		++it;
	}

	Eigen::MatrixXd mat_u = W * svd.matrixU().block(0, 0, m_b, nu);
	Eigen::MatrixXd mat_v = V * svd.matrixV().transpose().block(0, 0, m_b, nu);
	Eigen::VectorXd eigen_singular_values = svd.singularValues().segment(0, nu);
	return std::make_tuple(mat_u, eigen_singular_values, mat_v);
}