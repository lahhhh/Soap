#pragma once

#include "Custom.h"

bool is_normal_distribution_anderson_darling(const Eigen::ArrayXd& data, double significance_level = 0.05);

bool is_variance_homogeneous_bartlett_test(const std::vector<Eigen::ArrayXd>& data, double significance_level = 0.05);

bool anova(const std::vector<Eigen::ArrayXd>& groups, double significance_level = 0.05);

bool Welch_anova(const std::vector<Eigen::ArrayXd>& groups, double significance_level = 0.05);

bool Kruskal_Wallis(const std::vector<Eigen::ArrayXd>& groups, double significance_level = 0.05);

struct DUNN_RESULT {
	std::vector<std::tuple<int, int, double>> r;
};

DUNN_RESULT Dunn_test(const std::vector<Eigen::ArrayXd>& groups);

double t_test(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y);
