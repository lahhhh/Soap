#include "SignRank.h"

#include <algorithm>
#include <cmath>

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

static double* w;
static int allocated_n;

static void
w_free(void)
{
    if (!w) return;

    free((void*)w);
    w = 0;
    allocated_n = 0;
}

void signrank_free(void)
{
    w_free();
}

static void
w_init_maybe(int n)
{
    int u, c;

    u = n * (n + 1) / 2;
    c = (u / 2);

    if (w) {
        if (n != allocated_n) {
            w_free();
        }
        else return;
    }


    if (!w) {
        w = (double*)calloc((size_t)c + 1, sizeof(double));
        allocated_n = n;
    }
}


static double
csignrank(int k, int n)
{
    int u = n * (n + 1) / 2;
    int c = (u / 2);

    if (k < 0 || k > u) {
        return 0;
    }
    if (k > c) {
        k = u - k;
    }

    if (n == 1) {
        return 1.0;
    }
    if (w[0] == 1.) {
        return w[k];
    }

    w[0] = w[1] = 1.0;

    for (int j = 2; j < n + 1; ++j) {

        int end = std::min(j * (j + 1) / 2, c);

        for (int i = end; i >= j; --i) {
            w[i] += w[i - j];
        }
    }

    return w[k];
}


double psignrank(double x, double n, bool lower_tail)
{
    int i;
    double f, p;

    n = std::nearbyint(n);
    if (n <= 0) return 1;

    x = std::nearbyint(x + 1e-7);
    if (x < 0.0) {
        return(lower_tail ? 0.0 : 1.0);
    }
    if (x >= n * (n + 1) / 2) {
        return(lower_tail ? 1.0 : 0.0);
    }

    int nn = (int)n;
    w_init_maybe(nn);
    f = std::exp(-n * M_LN2);
    p = 0;
    if (x <= (n * (n + 1) / 4)) {
        for (i = 0; i <= x; i++) {
            p += csignrank(i, nn) * f;
        }
    }
    else {
        x = n * (n + 1) / 2 - x;
        for (i = 0; i < x; i++)
            p += csignrank(i, nn) * f;
        lower_tail = !lower_tail; /* p = 1 - p; */
    }

    return lower_tail ? p : (1.0 - p);
} /* psignrank() */
