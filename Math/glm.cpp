#include "glm.h"

#include "dqrdc2.h"
#include "custom.h"
#include "chol2inv.h"
#include "pval.h"

GLM_RES glm_gaussian(
	const Eigen::ArrayXd& y,
	const Eigen::MatrixXd& x,
	int maxit
) {
	bool conv{ false };
	bool boundary{ false };

	Eigen::MatrixXd _x = Eigen::MatrixXd::Ones(x.rows(), x.cols() + 1);
	_x.block(0, 1, x.rows(), x.cols()) = x;

	int nobs = y.size();
	int nvars = _x.cols();

	Eigen::ArrayXd mustart = y;
	Eigen::ArrayXd eta = y;
	Eigen::ArrayXd mu = eta;
	Eigen::ArrayXd coef, coefold;

	double devold{ 0.0 };

	DQRLS_RES fit;

	double dev{ 0.0 };

	int iter{ 0 };
	for (; iter < maxit; ++iter) {
		Eigen::ArrayXd z = eta + y - mu;

		fit = cdqrls(_x, z, 1e-11);

		if (fit.coefficients.isInf().any()) {			
			break;
		}

		if (nobs < fit.rank) {
			break;
		}

		const int n_coef = fit.coefficients.size();
		Eigen::ArrayXd start(n_coef);
		for (int i = 0; i < n_coef; ++i) {
			start[fit.pivot[i]] = fit.coefficients[i];
		}
		eta = (_x * start.matrix()).col(0);
		mu = eta;
		dev = (y - mu).matrix().squaredNorm();

		boundary = false;

		if (std::isinf(dev)) {
			// impossible dev -> mu -> eta ->start ->fit.coefficients
			break;
		}

		if (std::abs(dev - devold) / (0.1 + std::abs(dev)) < 1e-8) {
			conv = true;
			coef = start;
			break;
		}
		else {
			devold = dev;
			coef = start;
		}
	}

	if (fit.rank < nvars) {
		for (int i = fit.rank; i < nvars; ++i) {
			coef[fit.pivot[i]] = std::nan("0");
		}
	}

	Eigen::ArrayXd residuals = y - mu;
	//Eigen::MatrixXd Rmat = fit.qr.block(0, 0, nvars, nvars);
	//for (int i = 0; i < nvars; ++i) {
	//	for (int j = 0; j < nvars; ++j) {
	//		if (i > j) {
	//			Rmat(i, j) = 0.0;
	//		}
	//	}
	//}
	double wtdmu = y.sum() / y.size();
	double nulldev = (y - wtdmu).matrix().squaredNorm();
	int nok = nobs;
	int nulldf = nok - 1;
	int resdf = nok - fit.rank;

	GLM_RES res;
	res.coefficients = coef;
	res.residuals = residuals;
	res.fitted_values = mu;
	res.effects = fit.effects;
	res.rank = fit.rank;
	res.qr.qr = fit.qr;
	res.qr.rank = fit.rank;
	res.qr.qraux = fit.qraux;
	res.qr.pivot = fit.pivot;
	res.qr.tol = fit.tol;
	res.linear_predictors = eta;
	res.deviance = dev;
	res.null_deviance = nulldev;
	res.iter = iter + 1;
	res.df_residual = resdf;
	res.df_null = nulldf;
	res.y = y;
	res.converged = conv;
	res.boundary = boundary; // always false

	return res;
}

GLM_SUMMARY_RES glm_summary(
	const GLM_RES& object
) {
	int df_r = object.df_residual;

	GLM_SUMMARY_RES res;

	if (df_r <= 0) {
		res.hasnan = true;
		return res;
	}

	double dispersion = object.residuals.matrix().squaredNorm() / df_r;

	int p = object.rank;

	// simplified for porting
	if (p > 0) {
		auto Qr = object.qr;
		Eigen::ArrayXd coef_p = _Cs reordered(object.coefficients, _Cs reordered(Qr.pivot, _Cs seq_n(0, p)));
		Eigen::MatrixXd covmat_unscaled = chol2inv(Qr.qr.block(0, 0, p, p));
		Eigen::MatrixXd covmat = covmat_unscaled * dispersion;
		Eigen::ArrayXd var_cf = covmat.diagonal();
		Eigen::ArrayXd s_err = var_cf.sqrt();
		Eigen::ArrayXd t_value = coef_p / s_err;

		Eigen::ArrayXd p_value = 2.0 * p_students_t<Eigen::ArrayXd>(-t_value.abs(), df_r);

		res.estimate = coef_p;
		res.std_err = s_err;
		res.t_value = t_value;
		res.p = p_value;

		return res;
	}
	else {
		res.hasnan = true;
		return res;
	}
}