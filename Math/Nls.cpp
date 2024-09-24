#include "Nls.h"

constexpr double DERIV_STEP = 1e-5;
constexpr int MAX_ITER = 100;

double derive(double (*func)(double, const Eigen::ArrayXd&), double input, Eigen::ArrayXd params, int n) {
	Eigen::ArrayXd params1 = params;
	params[n] -= DERIV_STEP;
	double r1 = func(input, params);
	params[n] += 2 * DERIV_STEP;
	double r2 = func(input, params);
	return (r2 - r1) / (2 * DERIV_STEP);
};

void nls_LevenbergMarquardt(
	double (*func)(double, const Eigen::ArrayXd&),
	const Eigen::ArrayXd& inputs, 
	const Eigen::ArrayXd& outputs,
	Eigen::ArrayXd& params) 
{
	int size_of_inputs = inputs.size(), size_of_parameters = params.size();
	Eigen::ArrayXd r(size_of_inputs), r_tmp(size_of_inputs), params_tmp(size_of_parameters);
	Eigen::MatrixXd jf(size_of_inputs, size_of_parameters);
	double input;

	double last_mse = 0;
	double u = 1., v = 2.;
	Eigen::MatrixXd I = Eigen::MatrixXd::Ones(size_of_parameters, size_of_parameters);

	for (int i = 0; i < MAX_ITER; ++i) {
		double mse = 0., mse_temp = 0.;

		for (int j = 0; j < size_of_inputs; ++j) {
			input = inputs[j];

			r[j] = outputs[j] - func(input, params);
			mse += r[j] * r[j];

			for (int k = 0; k < size_of_parameters; ++k) {
				jf(j, k) = derive(func, input, params, k);
			}
		}

		mse /= size_of_inputs;
		params_tmp = params;

		Eigen::MatrixXd hlm = (jf.transpose() * jf + u * I).inverse() * jf.transpose() * r.matrix();
		params_tmp += hlm.array();

		for (int j = 0; j < size_of_inputs; ++j) {
			r_tmp[j] = outputs[j] - func(input, params_tmp);
			mse_temp += r_tmp[j] * r_tmp[j];
		}
		mse_temp /= size_of_inputs;

		Eigen::MatrixXd q = 0.5 * hlm.transpose() * (u * hlm - jf.transpose() * r.matrix());
		double q_val = (mse - mse_temp) / q(0, 0);

		if (q_val > 0) {
			double s = 1. / 3.;
			v = 2;
			mse = mse_temp;
			params = params_tmp;
			double temp = 1 - std::pow(2 * q_val - 1, 3);
			if (s > temp) {
				u *= s;
			}
			else {
				u *= temp;
			}
		}
		else {
			u *= v;
			v *= 2;
			params = params_tmp;
		}

		if (fabs(mse - last_mse) < 1e-8) {
			break;
		}
		last_mse = mse;
	}
};
