#include "Poisson.h"
#include <math.h>
#include <cmath>

static const int LSTEP = 200;
static const double EXPTHRES = exp(LSTEP);
static const double EXPSTEP = exp(-LSTEP);
static const int bigx = 20;

constexpr double M_LN10 = 2.30258509299404568402;

double __log_space_add(double log_x, double log_y) {
	if (log_x > log_y) {
		return log_x + log1p(exp(log_y - log_x));
	}
	else {
		return log_y + log1p(exp(log_x - log_y));
	}
};

double __poisson_cdf(int n, double lambda) {
	if (n < 0)return 0;

	double next_cdf = exp(-lambda);
	double cdf = next_cdf;

	for (int i = 1; i < n + 1; ++i) {
		double last_cdf = next_cdf;
		next_cdf = last_cdf * lambda / i;
		cdf += next_cdf;
	}
	if (cdf > 1) {
		return 1;
	}
	else {
		return cdf;
	}
};

double __poisson_cdf_q(int n, double lambda) {
	if (n < 0)return 1;
	double next_cdf = exp(-lambda);

	for (int i = 1; i < n + 1; ++i) {
		double last_cdf = next_cdf;
		next_cdf = last_cdf * lambda / i;
	}
	double cdf = 0;
	int i = n + 1;
	while (next_cdf > 0) {
		double last_cdf = next_cdf;
		next_cdf = last_cdf * lambda / i;
		cdf += next_cdf;
		++i;
	}
	return cdf;
};

double __poisson_cdf_large_lambda(int n, double lambda) {
	if (n < 0)return 0;
	int num_parts = (int)(lambda / LSTEP);
	double last_exp = exp(-fmod(lambda, LSTEP));
	double next_cdf = EXPSTEP;
	double cdf = 0;
	--num_parts;
	for (int i = 1; i < n + 1; ++i) {
		double last_cdf = next_cdf;
		next_cdf = last_cdf * lambda / i;
		cdf += next_cdf;
		if (next_cdf > EXPTHRES || cdf > EXPTHRES) {
			if (num_parts >= 1) {
				cdf *= EXPSTEP;
				next_cdf *= EXPSTEP;
				--num_parts;
			}
			else {
				cdf *= last_exp;
				last_exp = 1;
			}
		}

	}
	for (int i = 0; i < num_parts; ++i) {
		cdf *= EXPSTEP;
	}
	cdf *= last_exp;
	return cdf;
}

double __poisson_cdf_q_large_lambda(int n, double lambda) {
	if (n < 0)return 1; //???

	int num_parts = (int)(lambda / LSTEP);
	double last_exp = exp(-fmod(lambda, LSTEP));
	double next_cdf = EXPSTEP;
	double cdf = 0;
	--num_parts;
	for (int i = 1; i < n + 1; ++i) {
		double last_cdf = next_cdf;
		next_cdf = last_cdf * lambda / i;
		if (next_cdf > EXPTHRES) {
			if (num_parts > 1) {
				next_cdf *= EXPSTEP;
				--num_parts;
			}
			else {
				return 0;//error
			}
		}
	}
	int i = n + 1;
	while (next_cdf > 0) {
		double last_cdf = next_cdf;
		next_cdf = last_cdf * lambda / i;
		cdf += next_cdf;
		++i;
		if (next_cdf > EXPTHRES || cdf > EXPTHRES) {
			if (num_parts >= 1) {
				cdf *= EXPSTEP;
				next_cdf *= EXPSTEP;
				--num_parts;
			}
			else {
				cdf *= last_exp;
				last_exp = 1;
			}
		}
	}
	for (int i = 0; i < num_parts; ++i) {
		cdf *= EXPSTEP;
	}
	cdf *= last_exp;
	return cdf;
};

double __log10_poisson_cdf_p_large_lambda(int n, double lambda) {
	double residue = 0, log_x = 0, ln_lambda = log(lambda);
	int m = n, i = 0;
	double sum_ln_m = 0;
	for (int j = 1; j < m + 1; ++j) {
		sum_ln_m += log(j);
	}
	log_x = m * ln_lambda - sum_ln_m;
	residue = log_x;
	while (m > 1) {
		--m;
		double log_y = log_x - ln_lambda + log(m);
		double pre_residue = residue;
		residue = __log_space_add(pre_residue, log_y);
		if (fabs(pre_residue - residue) < 1e-10) {
			break;
		}
		log_x = log_y;
	}
	return (residue - lambda) / M_LN10;
};

double __log10_poisson_cdf_q_large_lambda(int n, double lambda) {
	double residue = 0, log_x = 0, ln_lambda = log(lambda);
	int m = n + 1;
	double sum_ln_m = 0;
	for (int i = 1; i < m + 1; ++i) {
		sum_ln_m += log(i);
	}
	log_x = m * ln_lambda - sum_ln_m;
	residue = log_x;
	while (true) {
		++m;
		double log_y = log_x + ln_lambda - log(m);
		double pre_residue = residue;
		residue = __log_space_add(pre_residue, log_y);
		if (fabs(pre_residue - residue) < 1e-5) {
			break;
		}
		log_x = log_y;
	}
	return (residue - lambda) / M_LN10;
};