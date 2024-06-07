#include "WilcoxTest.h"

#include "Custom.h"
#include "pval.h"
#include "SignRank.h"

double wilcox_test(
	const Eigen::ArrayXd& x,
	const Eigen::ArrayXd& y,
	bool paired,
	int alternative,
	double mu,	
	bool correct
) {
	double correction = 0;

	if (paired) {
		Eigen::ArrayXd query = x - y - mu;

		Eigen::ArrayX<bool> not_zero = (query != 0);
		bool has_zero = std::ranges::any_of(not_zero, [](auto no) {return !no; });

		if (has_zero) {
			query = _Cs sliced(query, not_zero);
		}

		double n = query.size();

		if (n == 0) {
			return 1.0; // should return std::nan 
		}

		bool exact = n < 50;

		auto [r, ties, n_tie] = _Cs rank_average(query.abs()); // ntie = sum(NTIES^3 - NTIES)
		double statistic = r(_Cs which(query > 0)).sum();

		if (exact && !ties && !has_zero) {
			if (alternative == 0) {
				double p = statistic > (n * (n + 1) / 4) ? 
					psignrank(statistic - 1, n, false) : psignrank(statistic, n, true);
				return std::min(2 * p, 1.0);
			}
			else if (alternative == 2) {
				return psignrank(statistic - 1, n, false);
			}
			else /* alternative == 1 */ {
				return psignrank(statistic, n, true);
			}
		}
		else {
			
			double z = statistic - n * (n + 1) / 4;
			double sigma = std::sqrt(n * (n + 1) * (2 * n + 1) / 24.0 - n_tie / 48.0);

			if (sigma == 0.0)return 1.0; // should return std::nan

			if (correct) {
				if (alternative == 0) {
					correction = _Cs sign(z) * 0.5;
				}
				else if (alternative == 2) {
					correction = 0.5;
				}
				else /* alternative == 1 */ {
					correction = -0.5;
				}
			}

			z = (z - correction) / sigma;

			if (std::isnan(z)) {
				return 1.0;
			}

			if (alternative == 0) {
				return 2 * std::min(p_normal(z), p_normal(z, 0.0, 1.0, false));
			}
			else if (alternative == 2) {
				return p_normal(z, 0.0, 1.0, false);
			}
			else /* alternative == 1 */ {
				return p_normal(z);
			}
		}
	}
	else {
		const int n_x = x.size();
		const int n_y = y.size();
		Eigen::ArrayXd r(n_x + n_y);
		r << x - mu, y;
		auto [rank, ties, n_tie] = _Cs rank_average(r);
		bool exact = (n_x < 50) && (n_y < 50);
		double statistic = rank.segment(0, n_x).sum() - n_x * (n_x + 1) / 2.0;

		if (exact && !ties) {
			if (alternative == 0) {
				double p = statistic > ((n_x * n_y) / 2.0) ?
					p_wilcox(statistic - 1, n_x, n_y, false) : p_wilcox(statistic, n_x, n_y);
				return std::min(2 * p, 1.0);
			}
			else if (alternative == 2) {
				return p_wilcox(statistic - 1, n_x, n_y, false);
			}
			else /* alternative == 1 */ {
				return p_wilcox(statistic, n_x, n_y);
			}
		}
		else {
			double z = statistic - n_x * n_y / 2.0;
			double sigma = std::sqrt((n_x * n_y / 12.0) * ((n_x + n_y + 1.0) - n_tie / ((n_x + n_y) * (n_x + n_y - 1.0))));
			
			if (sigma == 0.0)return 1.0; // should return std::nan

			if (correct) {
				if (alternative == 0) {
					correction = _Cs sign(z) * 0.5;
				}
				else if (alternative == 2) {
					correction = 0.5;
				}
				else /* alternative == 1 */ {
					correction = -0.5;
				}
			}

			z = (z - correction) / sigma;

			if (std::isnan(z)) {
				return 1.0;
			}

			if (alternative == 0) {
				return 2 * std::min(p_normal(z), p_normal(z, 0.0, 1.0, false));
			}
			else if (alternative == 2) {
				return p_normal(z, 0.0, 1.0, false);
			}
			else /* alternative == 1 */ {
				return p_normal(z);
			}
		}
	}
};

