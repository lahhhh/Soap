#pragma once

double __log_space_add(double logX, double logY);

double __poisson_cdf(int n, double lambda);

double __poisson_cdf_q(int n, double lambda);

double __poisson_cdf_large_lambda(int n, double lambda);

double __poisson_cdf_q_large_lambda(int n, double lambda);

double __log10_poisson_cdf_p_large_lambda(int n, double lambda);

double __log10_poisson_cdf_q_large_lambda(int n, double lambda);

// modified from macs3 source code
template <bool lower = false, bool log10 = false>
double poisson_cdf(int n, double lambda) {
	if (lambda <= 0) {
		return 0;//error
	}
	if (log10) {
		if (lower) {
			return __log10_poisson_cdf_p_large_lambda(n, lambda);
		}
		else {
			return __log10_poisson_cdf_q_large_lambda(n, lambda);
		}
	}
	if (lower) {
		if (lambda > 700) {
			return __poisson_cdf_large_lambda(n, lambda);
		}
		else {
			return __poisson_cdf(n, lambda);
		}
	}
	else {
		if (lambda > 700) {
			return __poisson_cdf_q_large_lambda(n, lambda);
		}
		else {
			return __poisson_cdf_q(n, lambda);
		}
	}
};