#include "dqrdc2.h"

DQRLS_RES cdqrls(
	const Eigen::MatrixXd& x,
	const Eigen::ArrayXd& y,
	double tol
) {

	const int n = x.rows();
	const int p = x.cols();

	Eigen::MatrixXd _x = x;
	Eigen::ArrayXd _y = y;
	Eigen::ArrayXd coefficients(p);
	Eigen::ArrayXd residuals = _y;
	Eigen::ArrayXd effects = _y;
	int rank;
	Eigen::ArrayXi pivot = Eigen::ArrayXi::LinSpaced(p, 0, p - 1);
	Eigen::ArrayXd qraux(p);

	Eigen::MatrixXd work(2, p);
	dqrls(
		_x,
		_y,
		tol,
		coefficients,
		residuals,
		effects,
		rank,
		pivot,
		qraux,
		work
	);

	bool pivoted{ false };
	for (int i = 0; i < p; ++i) {
		if (pivot[i] != i) {
			pivot = true;
			break;
		}
	}

	DQRLS_RES ans;
	ans.qr = _x;
	ans.coefficients = coefficients;
	ans.residuals = residuals;
	ans.effects = effects;
	ans.rank = rank;
	ans.pivot = pivot;
	ans.qraux = qraux;
	ans.tol = tol;
	ans.pivoted = pivoted;
	return ans;
}

void dqrls(
	Eigen::MatrixXd& x,
	Eigen::ArrayXd& y,
	double tol,
	Eigen::ArrayXd& b,
	Eigen::ArrayXd& rsd,
	Eigen::ArrayXd& qty,
	int& k,
	Eigen::ArrayXi& jpvt,
	Eigen::ArrayXd& qraux,
	Eigen::MatrixXd& work
) {
	dqrdc2(
		x,
		tol,
		k,
		qraux,
		jpvt,
		work
	);

	const int n = x.rows();
	const int p = x.cols();

	if (k > 0) {
		dqrsl(
			x,
			n,
			k,
			qraux,
			y,
			rsd,
			qty,
			b,
			rsd,
			rsd,
			1110
		);
	}
	else {
		for (int i = 0; i < n; ++i) {
			rsd(i) = y(i);
		}
	}

	for (int j = k; j < p; ++j) {
		b(j) = 0.0;
	}
}

void dqrdc2(
	Eigen::MatrixXd& x,
	double tol,
	int& k,
	Eigen::ArrayXd& qraux,
	Eigen::ArrayXi& jpvt,
	Eigen::MatrixXd& work
) {
	const int n = x.rows();
	const int p = x.cols();

	qraux = x.colwise().norm();
	work.row(0) = qraux;
	work.row(1) = qraux;

	for (int j = 0; j < p; ++j) {
		if (work(j, 1) == 0.0) {
			work(j, 1) = 1.0;
		}
	}

	int lup = std::min(n, p);
	k = p + 1;

	for (int l = 0; l < lup; ++l) {

		while (l < k - 1 && qraux[l] < work(l, 1) * tol) {
			for (int i = 0; i < n; ++i) {
				double t = x(i, l);

				for (int j = l + 1; j < p; ++j) {
					x(i, j - 1) = x(i, j);
				}
				x(i, p - 1) = t;
			}

			int i = jpvt[l];
			double t = qraux[l];
			double tt = work(l, 0);
			double ttt = work(l, 1);

			for (int j = l + 1; j < p; ++j) {
				jpvt[j - 1] = jpvt[j];
				qraux[j - 1] = qraux[j];
				work(j - 1, 0) = work(j, 0);
				work(j - 1, 1) = work(j, 1);
			}

			jpvt[p - 1] = i;
			qraux[p - 1] = t;
			work(p - 1, 0) = tt;
			work(p - 1, 1) = ttt;
			--k;
		}

		if (l == n - 1) {
			break;
		}

		double nrmxl = x.col(l).segment(l, n - l).norm();
		if (nrmxl > 0.0) {
			if (x(l, l) != 0.0) {
				nrmxl = x(l, l) < 0.0 ? -std::abs(nrmxl) : std::abs(nrmxl);
			}
			x.col(l).segment(l, n - l) /= nrmxl;
			x(l, l) = x(l, l) + 1.0;

			if (p >= l + 2) {
				for (int j = l + 1; j < p; ++j) {
					double t = -(x.col(l).segment(l, n - l).array() * x.col(j).segment(l, n - l).array()).sum() / x(l, l);

					x.col(j).segment(l, n - l) += t * x.col(l).segment(l, n - l);

					if (qraux[j] != 0.0) {
						double tt = 1.0 - (std::abs(x(l, j)) / qraux[j]) * (std::abs(x(l, j)) / qraux[j]);
						tt = std::max(tt, 0.0);
						t = tt;

						if (std::abs(t) < 1e-6) {
							qraux[j] = x.col(j).segment(l + 1, n - l - 1).norm();
							work(j, 0) = qraux[j];
						}
						else {
							qraux[j] *= std::sqrt(t);
						}
					}
				}
			}

			qraux[l] = x(l, l);
			x(l, l) = -nrmxl;
		}
	}

	k = std::min(k - 1, n);
};

void dqrsl(
	Eigen::MatrixXd& x,
	int n,
	int& k,
	Eigen::ArrayXd& qraux,
	Eigen::ArrayXd& y,
	Eigen::ArrayXd& qy,
	Eigen::ArrayXd& qty,
	Eigen::ArrayXd& b,
	Eigen::ArrayXd& rsd,
	Eigen::ArrayXd& xb,
	int job
) {

	bool cqy = (job / 10000) != 0;
	bool cqty = (job % 10000) != 0;
	bool cb = ((job % 1000) / 100) != 0;
	bool cr = ((job % 100) / 10) != 0;
	bool cxb = (job % 10) != 0;

	int ju = std::min(k, n - 1);

	if (ju == 0) {
		if (cqy) {
			qy[0] = y[0];
		}
		if (cqty) {
			qty[0] = y[0];
		}
		if (cxb) {
			xb[0] = y[0];
		}
		if (cb) {
			if (x(0, 0) != 0.0) {
				b[0] = y[0] / x(0, 0);
			}
		}
		if (cr) {
			rsd[0] = 0.0;
		}

		return;
	}

	if (cqy) {
		qy = y;
	}

	if (cqty) {
		qty = y;
	}

	if (cqy) {
		for (int jj = 0; jj < ju; ++jj) {
			int j = ju - jj - 1;
			if (qraux[j] != 0.0) {
				double temp = x(j, j);
				x(j, j) = qraux[j];
				double t = -(x.col(j).segment(j, n - j).array() * qy.segment(j, n - j).array()).sum() / x(j, j);
				qy.segment(j, n - j).array() += t * x.col(j).segment(j, n - j).array();
				x(j, j) = temp;
			}
		}
	}

	if (cqty) {
		for (int j = 0; j < ju; ++j) {
			if (qraux[j] != 0.0) {
				double temp = x(j, j);
				x(j, j) = qraux[j];
				double t = -(x.col(j).segment(j, n - j).array() * qty.segment(j, n - j).array()).sum() / x(j, j);
				qty.segment(j, n - j).array() += t * x.col(j).segment(j, n - j).array();
				x(j, j) = temp;
			}
		}
	}

	if (cb) {
		b.segment(0, k) = qty.segment(0, k);
	}
	if (cxb) {
		xb.segment(0, k) = qty.segment(0, k);
	}
	if (cr && k < n) {
		rsd.segment(k, n - k) = qty.segment(k, n - k);
	}
	if (cxb && k <= n - 1) {
		for (int i = k; i < n; ++i) {
			xb(i) = 0.0;
		}
	}

	if (cr) {
		for (int i = 0; i < k; ++i) {
			rsd(i) = 0.0;
		}
	}

	if (cb) {
		for (int jj = 0; jj < k; ++jj) {
			int j = k - jj - 1;
			if (x(j, j) == 0.0) {
				break;
			}

			b(j) = b(j) / x(j, j);
			if (j != 0) {
				double t = -b(j);
				b.segment(0, j).array() += t * x.col(j).segment(0, j).array();
			}
		}
	}

	if (cr || cxb) {
		for (int jj = 0; jj < ju; ++jj) {
			int j = ju - jj - 1;
			if (qraux[j] != 0.0) {
				double temp = x(j, j);
				x(j, j) = qraux[j];
				if (cr) {
					double t = -(x.col(j).segment(j, n - j).array() * rsd.segment(j, n - j).array()).sum() / x(j, j);
					rsd.segment(j, n - j).array() += t * x.col(j).segment(j, n - j).array();
				}
				if (cxb) {
					double t = -(x.col(j).segment(j, n - j).array() * xb.segment(j, n - j).array()).sum() / x(j, j);
					xb.segment(j, n - j).array() += t * x.col(j).segment(j, n - j).array();
				}
				x(j, j) = temp;
			}
		}
	}
}