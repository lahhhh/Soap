#include "glasso.h"

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> glasso(
	const Eigen::MatrixXd& sss,
	const Eigen::MatrixXd& rrho,
	double threshold,
	int maxit,
	bool approx,
	bool penalize_diagonal
) {

	int n1 = sss.rows();

	Eigen::MatrixXd www = Eigen::MatrixXd::Zero(n1, n1);
	Eigen::MatrixXd wwwi = Eigen::MatrixXd::Zero(n1, n1);

	double ddel{ 0.0 };
	int nniter{ 0 };

	if (approx) {

		lasinv1(
			sss,
			rrho,
			approx,
			0,
			threshold,
			maxit,
			penalize_diagonal,
			www,
			wwwi,
			nniter,
			ddel
		);

		return std::make_pair(www, wwwi);
	}

	Eigen::MatrixXi ic = Eigen::MatrixXi::Zero(2, n1);
	Eigen::ArrayXi ir = Eigen::ArrayXi::Zero(n1);
	Eigen::ArrayXi ie = Eigen::ArrayXi::Zero(n1);
	int nc{ 0 };
	connect(sss, rrho, nc, ic, ir, ie);

	Eigen::MatrixXd ss, rho, ww, wwi;
	int l{ 0 };

	for (int kc = 0; kc < nc; ++kc) {

		int n = ic(1, kc) - ic(0, kc) + 1;

		if (n <= 1) {

			int k = ir[ic(0, kc) - 1];
			www.col(k).setZero();
			www.row(k).setZero();
			wwwi.col(k).setZero();
			wwwi.row(k).setZero();
			continue;
		}
		ss = Eigen::MatrixXd::Zero(n, n);
		rho = Eigen::MatrixXd::Zero(n, n);
		ww = Eigen::MatrixXd::Zero(n, n);
		wwi = Eigen::MatrixXd::Zero(n, n);

		int kb = ic(0, kc);
		int ke = ic(1, kc);
		l = 0;
		for (int k = kb; k < ke + 1; ++k, ++l) {
			int ik = ir(k - 1);

			for (int j = kb; j < ke + 1; ++j) {
				int ij = ir(j - 1);

				ss(j - kb, l) = sss(ij, ik);
				rho(j - kb, l) = rrho(ij, ik);
				ww(j - kb, l) = www(ij, ik);
				wwi(j - kb, l) = wwwi(ij, ik);
			}
		}

		int niter{ 0 };
		double del{ 0.0 };
		lasinv1(
			ss,
			rho,
			approx,
			0,
			threshold,
			maxit,
			penalize_diagonal,
			ww,
			wwi,
			niter,
			del
		);
		nniter += niter;
		ddel += del;

		for (int j = kb; j < ke + 1; ++j) {
			int k = ir[j - 1];
			www.col(k).setZero();
			www.row(k).setZero();
			wwwi.col(k).setZero();
			wwwi.row(k).setZero();
		}
		l = 0;
		for (int k = kb; k < ke + 1; ++k) {
			int ik = ir(k - 1);
			for (int j = kb; j < ke + 1; ++j) {
				wwwi(ir[j - 1], ik) = wwi(j - kb, l);
			}
			++l;
		}
		if (!approx) {
			l = 0;
			for (int k = kb; k < ke + 1; ++k) {
				int ik = ir(k - 1);
				for (int j = kb; j < ke + 1; ++j) {
					www(ir[j - 1], ik) = ww(j - kb, l);
				}
				++l;
			}
		}
	}
	ddel /= nc;
	if (!approx) {
		for (int j = 0; j < n1; ++j) {
			if (www(j, j) != 0.0) {
				continue;
			}
			if (!penalize_diagonal) {
				www(j, j) = sss(j, j);
			}
			else {
				www(j, j) = sss(j, j) + rrho(j, j);
			}
			wwwi(j, j) = 1.0 / www(j, j);
		}
	}
	return std::make_pair(www, wwwi);
};

void connect(
	const Eigen::MatrixXd& ss,
	const Eigen::MatrixXd& rho,
	int& nc,
	Eigen::MatrixXi& ic,
	Eigen::ArrayXi& ir,
	Eigen::ArrayXi& ie
) {
	ie.setZero();
	nc = 0;
	int is{ 1 };

	int n = ss.rows();
	for (int k = 0; k < n; ++k) {
		if (ie[k] > 0) {
			continue;
		}

		ir[is - 1] = k;
		++nc;
		ie[k] = nc;
		ic(0, nc - 1) = is;
		++is;
		int na{ 0 };
		int nr{ 1 };
		{
			na = 0;
			for (int l = 0; l < nr; ++l) {
				int kk = ir[is + l - 2];
				for (int j = 0; j < n; ++j) {
					if (ie[j] > 0) {
						continue;
					}

					if (j == kk) {
						continue;
					}

					if (abs(ss(j, kk)) < rho(j, kk)) {
						continue;
					}

					++na;
					ir[is + na - 2] = j;
					ie[j] = nc;
				}
			}
		}

		if (na == 0) {
			ic(1, nc - 1) = is - 1;
			continue;
		}

		int il{ 0 };
		while (true) {
			int nas = na;
			int iss = is;
			il = iss + nas - 1;

			if (il >= n) {
				break;
			}

			is += na;
			nr = nas;
			{
				na = 0;
				for (int l = 0; l < nr; ++l) {
					int kk = ir[iss + l - 1];
					for (int j = 0; j < n; ++j) {
						if (ie[j] > 0) {
							continue;
						}

						if (j == kk) {
							continue;
						}

						if (abs(ss(j, kk)) < rho(j, kk)) {
							continue;
						}

						++na;
						ir[is + na - 2] = j;
						ie[j] = nc;
					}
				}

			}

			if (na == 0) {
				break;
			}
		}

		ic(1, nc - 1) = il;
	}
};

void lasinv1(
	const Eigen::MatrixXd& ss,
	const Eigen::MatrixXd& rho,
	bool approx,
	int start_type,
	double threshold,
	int maxit,
	bool penalize_diagonal,
	Eigen::MatrixXd& ww,
	Eigen::MatrixXd& wwi,
	int& niter,
	double& del
) {

	constexpr double eps{ 1e-7 };

	int n = ss.rows();

	int nm1 = n - 1;
	Eigen::MatrixXd vv = Eigen::MatrixXd::Zero(nm1, nm1);
	Eigen::MatrixXd xs;
	if (!approx) {
		xs = Eigen::MatrixXd::Zero(nm1, n);
	}

	Eigen::ArrayXd s = Eigen::ArrayXd::Zero(nm1);
	Eigen::ArrayXd so = Eigen::ArrayXd::Zero(nm1);
	Eigen::ArrayXd x = Eigen::ArrayXd::Zero(nm1);
	Eigen::ArrayXd z = Eigen::ArrayXd::Zero(nm1);
	Eigen::ArrayXd ro = Eigen::ArrayXd::Zero(nm1);
	Eigen::ArrayXi mm = Eigen::ArrayXi::Zero(nm1);

	Eigen::ArrayXd ws;
	if (!approx) {
		ws = Eigen::ArrayXd::Zero(n);
	}

	double shr{ 0.0 };
	for (int j = 0; j < n; ++j) {
		for (int k = 0; k < n; ++k) {
			if (j == k) {
				continue;
			}

			shr += std::abs(ss(j, k));
		}
	}
	if (shr == 0.0) {
		ww.setZero();
		wwi.setZero();
		for (int j = 0; j < n; ++j) {
			if (penalize_diagonal) {
				ww(j, j) = ss(j, j) + rho(j, j);
			}
			else {
				ww(j, j) = ss(j, j);
			}

			wwi(j, j) = 1.0 / std::max(ww(j, j), eps);
		}

		return;
	}

	shr = threshold * shr / nm1;

	if (approx) {
		if (start_type == 0) {
			wwi.setZero();
		}
		for (int m = 0; m < n; ++m) {
			setup(
				m,
				n,
				ss,
				rho,
				ss,
				vv,
				s,
				ro
			);

			int l{ 0 };
			for (int j = 0; j < n; ++j) {
				if (j == m) {
					continue;
				}

				x[l] = wwi(j, m);
				++l;
			}

			lasso(
				ro,
				nm1,
				vv,
				s,
				shr / n,
				x,
				z,
				mm
			);

			l = 0;
			for (int j = 0; j < n; ++j) {
				if (j == m) {
					continue;
				}

				wwi(j, m) = x[l];
				++l;
			}
		}
		niter = 1;
		return;
	}
	else {

		if (start_type == 0) {
			ww = ss;
			xs.setZero();
		}
		else {

			double xjj{ 0.0 };
			int l{ 0 };

			for (int j = 0; j < n; ++j) {
				xjj = -wwi(j, j);
				l = 0;
				for (int k = 0; k < n; ++k) {
					if (k == j) {
						continue;
					}

					xs(l, j) = wwi(k, j) / xjj;
					++l;
				}
			}
		}

		for (int j = 0; j < n; ++j) {
			if (!penalize_diagonal) {
				ww(j, j) = ss(j, j);
			}
			else {
				ww(j, j) = ss(j, j) + rho(j, j);
			}
		}

		niter = 0;
		double dlx{ 0.0 };

		while (true) {

			dlx = 0.0;

			for (int m = 0; m < n; ++m) {

				x = xs.col(m);
				ws = ww.col(m);

				setup(
					m,
					n,
					ss,
					rho,
					ww,
					vv,
					s,
					ro
				);

				so = s;

				double v = vv.cwiseAbs().sum();
				lasso(
					ro,
					nm1,
					vv,
					s,
					shr / v,
					x,
					z,
					mm
				);

				int l{ 0 };
				for (int j = 0; j < n; ++j) {
					if (j == m) {
						continue;
					}

					ww(j, m) = so[l] - s[l];
					ww(m, j) = ww(j, m);
					++l;
				}

				dlx = std::max(dlx, (ww.col(m).array() - ws).abs().sum());
				xs.col(m) = x;
			}

			++niter;
			if (niter >= maxit) {
				break;
			}

			if (dlx < shr) {
				break;
			}
		}

		del = dlx / nm1;
		inv(ww, xs, wwi);
	}
};

void setup(
	int m,
	int n,
	const Eigen::MatrixXd& ss,
	const Eigen::MatrixXd& rho,
	const Eigen::MatrixXd& ww,
	Eigen::MatrixXd& vv,
	Eigen::ArrayXd& s,
	Eigen::ArrayXd& r
) {

	int l{ 0 };
	for (int j = 0; j < n; ++j) {
		if (j == m) {
			continue;
		}

		r[l] = rho(j, m);
		s[l] = ss(j, m);

		int i{ 0 };

		for (int k = 0; k < n; ++k) {
			if (k == m) {
				continue;
			}

			vv(i, l) = ww(k, j);
			++i;
		}

		++l;
	}
};

void lasso(
	const Eigen::ArrayXd& rho,
	int n,
	const Eigen::MatrixXd& vv,
	Eigen::ArrayXd& s,
	double threshold,
	Eigen::ArrayXd& x,
	Eigen::ArrayXd& z,
	Eigen::ArrayXi& mm
) {

	fatmul(
		2, n, vv, x, s, z, mm
	);

	while (true) {

		double dlx{ 0.0 };

		for (int j = 0; j < n; ++j) {
			double xj = x[j];
			x[j] = 0.0;
			double t = s[j] + vv(j, j) * xj;

			double val = std::abs(t) - rho[j];
			if (val > 0.0) {
				x[j] = (t >= 0 ? val / vv(j, j) : -val / vv(j, j));
			}

			if (x[j] != xj) {
				double del = x[j] - xj;
				dlx = std::max(dlx, std::abs(del));
				s -= del * vv.col(j).array();
			}

		}

		if (dlx < threshold) {
			break;
		}
	}
};

void fatmul(
	int it,
	int n,
	const Eigen::MatrixXd& vv,
	const Eigen::ArrayXd& x,
	Eigen::ArrayXd& s,
	Eigen::ArrayXd& z,
	Eigen::ArrayXi& m
) {
	constexpr double fac = 0.2;

	int l{ 0 };
	for (int j = 0; j < n; ++j) {
		if (x[j] == 0.0) {
			continue;
		}

		m[l] = j;
		z[l] = x[j];
		++l;
	}

	if (l > int(fac * n)) {
		if (it == 1) {
			s = vv * x.matrix();
		}
		else {
			s -= (x.transpose().matrix() * vv).array();
		}
	}
	else {
		if (it == 1) {
			for (int j = 0; j < n; ++j) {
				s[j] = (vv.row(j)(m.segment(0, l)).array() * z.segment(0, l)).sum();
			}
		}
		else {
			for (int j = 0; j < n; ++j) {
				s[j] -= (vv.col(j)(m.segment(0, l)).array() * z.segment(0, l)).sum();
			}
		}
	}
};

void inv(
	const Eigen::MatrixXd& ww,
	Eigen::MatrixXd& xs,
	Eigen::MatrixXd& wwi
) {
	int n = ww.rows() - 1;
	xs = -xs;
	wwi(0, 0) = 1.0 / (ww(0, 0) + (xs.col(0).array() * ww.col(0).segment(1, n).array()).sum());
	wwi.col(0).segment(1, n) = wwi(0, 0) * xs.col(0);
	wwi(n, n) = 1.0 / (ww(n, n) + (xs.col(n).array() * ww.col(n).segment(0, n).array()).sum());
	wwi.col(n).segment(0, n) = wwi(n, n) * xs.col(n);

	for (int j = 1; j < n; ++j) {

		wwi(j, j) = 1.0 / (ww(j, j) +
			(xs.col(j).segment(0, j).array() * ww.col(j).segment(0, j).array()).sum() +
			(xs.col(j).segment(j, n - j).array() * ww.col(j).segment(j + 1, n - j).array()).sum()
			);
		wwi.col(j).segment(0, j) = wwi(j, j) * xs.col(j).segment(0, j);
		wwi.col(j).segment(j + 1, n - j) = wwi(j, j) * xs.col(j).segment(j, n - j);
	}
};