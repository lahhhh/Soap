#include "StatisticalTest.h"

#include "pval.h"

bool is_normal_distribution_anderson_darling(const Eigen::ArrayXd& data, double significance_level) {

	Eigen::ArrayXd x = custom::sorted(data);

	int n = x.size();
	if (n < 8) {
		return false;
	}

	double sd = custom::sd(x);
	if (sd == 0.0) {
		return false;
	}

	Eigen::ArrayXd logp1 = log(p_normal<Eigen::ArrayXd>((x - x.mean()) / sd));
	Eigen::ArrayXd logp2 = log(p_normal<Eigen::ArrayXd>(-(x - x.mean()) / sd));

	Eigen::ArrayXd h = (2 * Eigen::ArrayXd::LinSpaced(n, 1, n) - 1) * (logp1 + logp2.reverse());
	double A = -n - h.mean();
	double AA = (1 + 0.75 / n + 2.25 / (n * n)) * A;

	double p_val{ 3.7e-24 };

	if (AA < 0.2) {
		p_val = 1 - std::exp(-13.436 + 101.14 * AA - 223.73 * AA * AA);
	}
	else if (AA < 0.34) {
		p_val = 1 - std::exp(-8.318 + 42.796 * AA - 59.938 * AA * AA);
	}
	else if (AA < 0.6) {
		p_val = std::exp(0.9177 - 4.279 * AA - 1.38 * AA * AA);
	}
	else if (AA < 10) {
		p_val = std::exp(1.2937 - 5.709 * AA + 0.0186 * AA * AA);
	}

	return p_val < significance_level;
}

bool is_variance_homogeneous_bartlett_test(const std::vector<Eigen::ArrayXd>& data, double significance_level) {

	int k = data.size();
	if (k < 2) {
		return false;
	}

	if (custom::any(custom::sapply(data, [](auto&& t) {return t.size() < 2; }), true)) {
		return false;
	}

	auto n = custom::cast<Eigen::ArrayX>(custom::cast<double>(custom::sapply(data, [](auto&& t) {return t.size() - 1; })));
	int n_total = n.sum();

	auto v = custom::cast<Eigen::ArrayX>(custom::sapply(data, [](auto&& t) {return custom::var(t); }));
	if (custom::any(v, 0.0)) {
		return false;
	}
	auto v_total = (n * v).sum() / n_total;
	double STATISTIC = (n_total * std::log(v_total) - (n * log(v)).sum()) / (1 + ((1 / n).sum() - 1.0 / n_total) / (3 * (k - 1)));
	double p = p_chisq(STATISTIC, k - 1, false);

	return p > significance_level;
};

bool anova(const std::vector<Eigen::ArrayXd>& groups, double significance_level) {
	int num_groups = groups.size();
	int total_samples = 0;

	double grand_mean = 0.0;
	std::vector<double> means(num_groups);
	std::vector<int> group_sizes(num_groups);

	for (int i = 0; i < num_groups; ++i) {
		group_sizes[i] = groups[i].size();

		if (group_sizes[i] < 2) {
			return false;
		}

		means[i] = groups[i].mean();
		total_samples += groups[i].size();
		grand_mean += means[i] * groups[i].size();
	}
	grand_mean /= total_samples;

	double ssb = 0.0;
	for (int i = 0; i < num_groups; ++i) {
		ssb += group_sizes[i] * (means[i] - grand_mean) * (means[i] - grand_mean);
	}

	double ssw = 0.0;
	for (int i = 0; i < num_groups; ++i) {
		ssw += custom::var(groups[i]) * (group_sizes[i] - 1);
	}

	double dfb = num_groups - 1;
	double dfw = total_samples - num_groups;

	double msb = ssb / dfb;
	double msw = ssw / dfw;

	double f_statistic = msb / msw;

	boost::math::fisher_f_distribution<double> f_dist(dfb, dfw);
	double p_value = boost::math::cdf(boost::math::complement(f_dist, f_statistic));

	return p_value < significance_level;
}

bool Welch_anova(const std::vector<Eigen::ArrayXd>& groups, double significance_level) {
	int num_groups = groups.size();

	if (num_groups < 2) {
		return false;
	}

	// 计算每组均值、总体均值和总样本数
	Eigen::ArrayXd group_means(num_groups);
	Eigen::ArrayXd group_vars(num_groups);
	Eigen::ArrayXd group_sizes(num_groups);

	for (int i = 0; i < num_groups; ++i) {
		group_sizes[i] = groups[i].size();

		if (group_sizes[i] < 2) {
			return false;
		}

		group_means[i] = groups[i].mean();
		group_vars[i] = custom::var(groups[i]);

		if (group_vars[i] == 0.0) {
			return false; // do not handle this, because it can not be normal distributed
		}
	}

	Eigen::ArrayXd wi = group_sizes / group_vars;

	double sumwi = wi.sum();

	double tmp = ((1 - wi / sumwi) * (1 - wi / sumwi) / (group_sizes - 1)).sum() / (num_groups * num_groups - 1);
	double m = (wi * group_means).sum() / sumwi;
	double f_statistic = (wi * (group_means - m) * (group_means - m)).sum() / ((num_groups - 1) * (1 + 2 * (num_groups - 2) * tmp));

	double df_numerator = num_groups - 1;
	double df_denominator = 1 / (3 * tmp);

	boost::math::fisher_f_distribution<double> f_dist(df_numerator, df_denominator);
	double p_value = boost::math::cdf(boost::math::complement(f_dist, f_statistic));

	return p_value < significance_level;
}

bool Kruskal_Wallis(const std::vector<Eigen::ArrayXd>& groups, double significance_level) {

	int k = groups.size();
	if (k < 2) {
		return false;
	}

	auto x = custom::concatenated(groups);
	int n = x.size();
	if (n < 2) {
		return false;
	}

	auto [r, ties, ntie] = custom::rank_average(x);
	int ind{ 0 };
	double STATISTIC{ 0.0 };
	for (int i = 0; i < k; ++i) {
		int l = groups[i].size();
		double tmp = r.segment(ind, l).sum();
		STATISTIC += tmp * tmp / l;
		ind += l;
	}

	STATISTIC = ((12.0 * STATISTIC / (n * (n + 1)) - 3.0 * (n + 1)) / (1 - ntie / (n * n * n - n)));

	double p_val = p_chisq(STATISTIC, k - 1, false);

	return p_val < significance_level;
};

// with Bonferroni correction
DUNN_RESULT Dunn_test(const std::vector<Eigen::ArrayXd>& groups) {

	DUNN_RESULT res;

	int k = groups.size();
	if (k < 2) {
		return res;
	}

	auto x = custom::concatenated(groups);
	int n = x.size();
	if (n < 2) {
		return res;
	}

	auto sizes = custom::sapply(groups, [](auto&& t) {return t.size(); });

	auto [r, ties, ntie] = custom::rank_average(x);
	int ind{ 0 };
	Eigen::ArrayXd rank_mean(k);
	for (int i = 0; i < k; ++i) {
		int l = groups[i].size();
		rank_mean[i] = r.segment(ind, l).mean();
		ind += l;
	}
	boost::math::normal_distribution<> norm(0, 1);
	for (int i = 0; i < k; ++i) {
		for (int j = i + 1; j < k; ++j) {
			double numerator = rank_mean[i] - rank_mean[j];
			double denominator = std::sqrt((n * (n + 1) / 12.0 - ntie / ((n - 1) * 12)) * (1.0 / sizes[i] + 1.0 / sizes[j]));
			double z = numerator / denominator;
			double p = 1.0 - boost::math::cdf(norm, std::abs(z));

			p = std::max(1.0, p * i * (i - 1) / 2);

			res.r.emplace_back(i, j, p);
		}
	}

	return res;
};

double t_test(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y) {
	double mean1 = x.mean();
	double mean2 = y.mean();
	double sd1 = custom::sd(x);
	double sd2 = custom::sd(y);
	double n1 = x.size();
	double n2 = y.size();

	double pooled_se = std::sqrt((sd1 * sd1 / n1) + (sd2 * sd2 / n2));
	double t_stat = (mean1 - mean2) / pooled_se;

	double numerator = (sd1 * sd1 / n1 + sd2 * sd2 / n2) * (sd1 * sd1 / n1 + sd2 * sd2 / n2);
	double denominator = ((sd1 * sd1 / n1) * (sd1 * sd1 / n1) / (n1 - 1)) +
		((sd2 * sd2 / n2) * (sd2 * sd2 / n2) / (n2 - 1));
	double df = numerator / denominator;

	boost::math::students_t dist(df);
	double p_value = 2 * (1 - boost::math::cdf(dist, std::fabs(t_stat)));

	return p_value;
};