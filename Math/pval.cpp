#include "pval.h"

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#define WILCOX_MAX 50

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

static void
w_free(int m, int n);

static void
w_init_maybe(int m, int n);

static double*** w; /* to store  cwilcox(i,j,k) -> w[i][j][k] */
static int allocated_m, allocated_n;

static void
w_free(int m, int n)
{
	int i, j;

	for (i = m; i >= 0; i--) {
		for (j = n; j >= 0; j--) {
			if (w[i][j] != 0)
				free((void*)w[i][j]);
		}
		free((void*)w[i]);
	}
	free((void*)w);
	w = 0; allocated_m = allocated_n = 0;
}

static void
w_init_maybe(int m, int n)
{
	int i;

	if (m > n) {
		i = n;
		n = m;
		m = i;
	}
	if (w && (m > allocated_m || n > allocated_n))
		w_free(allocated_m, allocated_n); /* zeroes w */

	if (!w) { /* initialize w[][] */
		m = std::max(m, WILCOX_MAX);
		n = std::max(n, WILCOX_MAX);
		w = (double***)calloc((size_t)m + 1, sizeof(double**));
		for (i = 0; i <= m; ++i) {
			w[i] = (double**)calloc((size_t)n + 1, sizeof(double*));
		}
		allocated_m = m; allocated_n = n;
	}
}

static
double cwilcox(int k, int m, int n) {
	int c, u, i, j, l;
	u = m * n;
	if (k < 0 || k > u) {
		return 0;
	}
	c = (int)(u / 2);
	if (k > c) {
		k = u - k; /* hence  k <= floor(u / 2) */
	}
	if (m < n) {
		i = m;
		j = n;
	}
	else {
		i = n;
		j = m;
	} /* hence  i <= j */

	if (j == 0) { /* and hence i == 0 */
		return (k == 0);
	}

	if (j > 0 && k < j) {
		return cwilcox(k, i, k);
	}

	if (w[i][j] == 0) {
		w[i][j] = (double*)calloc((size_t)c + 1, sizeof(double));

		for (l = 0; l <= c; l++)
			w[i][j][l] = -1;
	}
	if (w[i][j][k] < 0) {
		if (j == 0) /* and hence i == 0 */
			w[i][j][k] = (k == 0);
		else
			w[i][j][k] = cwilcox(k - j, i - 1, j) + cwilcox(k, i, j - 1);

	}
	return(w[i][j][k]);
};

static double 
chebyshev_eval(double x, const double* a, const int n) // 0 < x < 10/31
{
	double b0 = 0, b1 = 0, b2 = 0, twox = x * 2;
	for (int i = 1; i <= n; ++i) {
		b2 = b1;
		b1 = b0;
		b0 = twox * b1 - b2 + a[n - i];
	}
	return (b0 - b2) * 0.5;
}

static double 
lgammacor(double x) // x should be more than 31 and less than 1e6
{
	constexpr static double algmcs[15] = {
	+.1666389480451863247205729650822e+0,
	-.1384948176067563840732986059135e-4,
	+.9810825646924729426157171547487e-8,
	-.1809129475572494194263306266719e-10,
	+.6221098041892605227126015543416e-13,
	-.3399615005417721944303330599666e-15,
	+.2683181998482698748957538846666e-17,
	-.2868042435334643284144622399999e-19,
	+.3962837061046434803679306666666e-21,
	-.6831888753985766870111999999999e-23,
	+.1429227355942498147573333333333e-24,
	-.3547598158101070547199999999999e-26,
	+.1025680058010470912000000000000e-27,
	-.3401102254316748799999999999999e-29,
	+.1276642195630062933333333333333e-30
	};

	double tmp = 10 / x;
	return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, 5) / x;
}

static double 
lbeta(double a, double b) { // a >= 31, b >= 31
	double p = a, q = a;
	if (b < p) p = b;
	if (b > q) q = b;
	double corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
	return std::log(q) * -0.5 + M_LN_SQRT_2PI + corr
		+ (p - 0.5) * std::log(p / (p + q)) + q * std::log1p(-p / (p + q));
}


static double 
lfastchoose(double n, double k) {
	return -std::log(n + 1.) - lbeta(n - k + 1., k + 1.);
};

static double 
choose(double n, double k) { // k -> integer --- this function is manually checked to simplify --- input : n > k , n > 0
	double r;
	if (k < 30) {
		if (n - k < k && n >= 0) {
			k = n - k;
		}
		if (k < 0) {
			return 0;
		}
		if (k == 0) {
			return 1;
		}
		r = n;
		for (int j = 2; j <= k; ++j) {
			r *= (n - j + 1) / j;
		}
		return r;
	}

	if (n < k) {
		return 0;
	}
	if (n - k < 30) {
		return choose(n, n - k);
	}
	return std::nearbyint(std::exp(lfastchoose(n, k)));

}

double p_wilcox(double q, double m, double n, bool lower_tail) {
	int i;
	double c, p;
	q = floor(q + 1e-7);
	if (q < 0.0) {
		return lower_tail ? 0 : 1;
	}
	if (q >= m * n) {
		return lower_tail ? 1 : 0;
	}
	int mm = (int)m, nn = (int)n;
	w_init_maybe(mm, nn);
	c = choose(m + n, n);
	p = 0;
	if (q <= (m * n / 2)) {
		for (i = 0; i <= q; ++i) {
			p += cwilcox(i, mm, nn) / c;
		}
		return lower_tail ? p : (1 - p);
	}
	else {
		q = m * n - q;
		for (i = 0; i < q; ++i)
			p += cwilcox(i, mm, nn) / c;
		return lower_tail ? (1 - p) : p;
	}
};


double p_normal(double x, double mu, double sigma, bool lower_tail) {

	boost::math::normal distribution(mu, sigma);

	double p = boost::math::cdf(distribution, x);
	
	return lower_tail ? p : (1.0 - p);
}

double dnorm(double x, double mu, double sigma) {

	boost::math::normal distribution(mu, sigma);

	double d = boost::math::pdf(distribution, x);

	return d;
};

Eigen::MatrixXd dnorm(const Eigen::MatrixXd& mat, double mu, double sigma) {

	int nrow = mat.rows(), ncol = mat.cols();

	Eigen::MatrixXd ret(nrow, ncol);

	boost::math::normal distribution(mu, sigma);

	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j) {
			ret(i, j) = boost::math::pdf(distribution, mat(i, j));
		}
	}

	return ret;
}

double p_students_t(double x, double n, bool lower_tail) {

	boost::math::students_t distribution(n);

	double p = boost::math::cdf(distribution, x);

	return lower_tail ? p : (1.0 - p);
};

double p_fisher_f(double f, double m, double n, bool lower_tail) {

	boost::math::fisher_f distribution(m, n);

	double p = boost::math::cdf(distribution, f);

	return lower_tail ? p : (1.0 - p);
};

double p_logistic(double q, double location, double scale, bool lower_tail) {

	boost::math::logistic distribution(location, scale);

	double p = boost::math::cdf(distribution, q);

	return lower_tail ? p : (1.0 - p);
};

double p_chisq(double q, int df, bool lower_tail) {

	boost::math::chi_squared chi_sq_dist(df);

	if (lower_tail) {
		return boost::math::cdf(chi_sq_dist, q);  
	}
	else {
		return boost::math::cdf(complement(chi_sq_dist, q));  
	}
}