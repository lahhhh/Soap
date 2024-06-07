#include "Phyper.h"

#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#include <math.h>
#include <limits>

double stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

	/*
	  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
	*/
	const static double sferr_halves[31] = {
		0.0, /* n=0 - wrong, place holder only */
		0.1534264097200273452913848,  /* 0.5 */
		0.0810614667953272582196702,  /* 1.0 */
		0.0548141210519176538961390,  /* 1.5 */
		0.0413406959554092940938221,  /* 2.0 */
		0.03316287351993628748511048, /* 2.5 */
		0.02767792568499833914878929, /* 3.0 */
		0.02374616365629749597132920, /* 3.5 */
		0.02079067210376509311152277, /* 4.0 */
		0.01848845053267318523077934, /* 4.5 */
		0.01664469118982119216319487, /* 5.0 */
		0.01513497322191737887351255, /* 5.5 */
		0.01387612882307074799874573, /* 6.0 */
		0.01281046524292022692424986, /* 6.5 */
		0.01189670994589177009505572, /* 7.0 */
		0.01110455975820691732662991, /* 7.5 */
		0.010411265261972096497478567, /* 8.0 */
		0.009799416126158803298389475, /* 8.5 */
		0.009255462182712732917728637, /* 9.0 */
		0.008768700134139385462952823, /* 9.5 */
		0.008330563433362871256469318, /* 10.0 */
		0.007934114564314020547248100, /* 10.5 */
		0.007573675487951840794972024, /* 11.0 */
		0.007244554301320383179543912, /* 11.5 */
		0.006942840107209529865664152, /* 12.0 */
		0.006665247032707682442354394, /* 12.5 */
		0.006408994188004207068439631, /* 13.0 */
		0.006171712263039457647532867, /* 13.5 */
		0.005951370112758847735624416, /* 14.0 */
		0.005746216513010115682023589, /* 14.5 */
		0.005554733551962801371038690  /* 15.0 */
	};
	double nn;

	if (n <= 15.0) {
		nn = n + n;
		return(sferr_halves[(int)nn]);
	}

	nn = n * n;
	if (n > 500) return((S0 - S1 / nn) / n);
	if (n > 80) return((S0 - (S1 - S2 / nn) / nn) / n);
	if (n > 35) return((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
	/* 15 < n <= 35 : */
	return((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
}

double bd0(double x, double np)
{

	if (fabs(x - np) < 0.1 * (x + np)) {
		double
			v = (x - np) / (x + np),  // might underflow to 0
			s = (x - np) * v;
		double ej = 2 * x * v;
		v *= v; // "v = v^2"
		for (int j = 1; j < 1000; ++j) { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
			ej *= v;// = 2 x v^(2j+1)
			double s_ = s;
			s += ej / ((j << 1) + 1);
			if (s == s_) { /* last term was effectively 0 */
				return s;
			}
		}

	}
	/* else:  | x - np |  is not too small */
	return(x * log(x / np) + np - x);
}

double dbinom_raw(double x, double n, double p, double q)
{
	if (p == 0) return((x == 0) ? 1 : 0);
	if (q == 0) return((x == n) ? 1 : 0);

	double lc;
	if (x == 0) {
		if (n == 0) return 1;
		lc = (p < 0.1) ? -bd0(n, n * q) - n * p : n * log(q);
		return(exp(lc));
	}
	if (x == n) {
		lc = (q < 0.1) ? -bd0(n, n * p) - n * q : n * log(p);
		return(exp(lc));
	}
	if (x < 0 || x > n) return(0);

	/* n*p or n*q can underflow to zero if n and p or q are small.  This
	   used to occur in dbeta, and gives NaN as from R 2.3.0.  */
	lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);

	/* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
	/* Upto R 2.7.1:
	 * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
	 * -- following is much better for  x << n : */
	double lf = M_LN_2PI + log(x) + log1p(-x / n);

	return exp(lc - 0.5 * lf);
}

double dhyper(double x, double r, double b, double n)
{
	double p, q, p1, p2, p3;

	if (x < 0) return(0);

	if (n < x || r < x || n - x > b) return(0);
	if (n == 0) return((x == 0) ? 1 : 0);

	p = ((double)n) / ((double)(r + b));
	q = ((double)(r + b - n)) / ((double)(r + b));

	p1 = dbinom_raw(x, r, p, q);
	p2 = dbinom_raw(n - x, b, p, q);
	p3 = dbinom_raw(n, r + b, p, q);

	return(p1 * p2 / p3);
}

static double pdhyper(double x, double NR, double NB, double n)
{
	/*
	 * Calculate
	 *
	 *	    phyper (x, NR, NB, n, TRUE, FALSE)
	 *   [log]  ----------------------------------
	 *	       dhyper (x, NR, NB, n, FALSE)
	 *
	 * without actually calling phyper.  This assumes that
	 *
	 *     x * (NR + NB) <= n * NR
	 *
	 */
	double sum = 0;
	double term = 1;

	while (x > 0 && term >= std::numeric_limits<double>::epsilon() * sum) {
		term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
		sum += term;
		x--;
	}

	double ss = (double)sum;
	return 1 + ss;
}

double phyper(double x, double NR, double NB, double n,
	bool lower_tail) // -> int
{
	/* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */

	double d, pd;

	if (x * (NR + NB) > n * NR) {
		/* Swap tails.	*/
		double oldNB = NB;
		NB = NR;
		NR = oldNB;
		x = n - x - 1;
		lower_tail = !lower_tail;
	}

	/* support of dhyper() as a function of its parameters
	 * R:  .suppHyper <- function(m,n,k) max(0, k-n) : min(k, m)
	 * --  where R's (m,n, k) == (NR,NB, n)  here */
	if (x < 0 || x < n - NB)
		return lower_tail ? 0 : 1;
	if (x >= NR || x >= n)
		return lower_tail ? 1 : 0;
	d = dhyper(x, NR, NB, n);
	// dhyper(.., log_p=FALSE) > 0 mathematically, but not always numerically :
	if (d == 0.)
		return lower_tail ? 0 : 1;
	pd = pdhyper(x, NR, NB, n);

	return lower_tail ? (d * pd) : (0.5 - (d * pd) + 0.5);
}
