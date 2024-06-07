#include "lm.h"

#include "dqrdc2.h"

LM_RES lm(
	const Eigen::ArrayXd& y,
	const Eigen::ArrayXd& x,
	bool fit_intercept
) {

	Eigen::MatrixXd _x = Eigen::MatrixXd::Ones(x.rows(), 2);

	if (fit_intercept) {
		_x = Eigen::MatrixXd::Ones(x.rows(), 2);
		_x.col(1) = x;
	}
	else {
		_x = x.matrix();
	}

	auto fit = cdqrls(
		_x,
		y,
		1e-7
	);

	const int n_coef = fit.coefficients.size();
	Eigen::ArrayXd coef(n_coef);
	for (int i = 0; i < n_coef; ++i) {
		coef[fit.pivot[i]] = fit.coefficients[i];
	}

	LM_RES res;
	res.coefficients = coef;

	return res;
};