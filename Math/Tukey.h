#pragma once

#include "Identifier.h"
#include <vector>

struct TUKEYHSD_RESULT {
	std::vector<std::tuple<int, int, double>> r;
};

TUKEYHSD_RESULT Tukey_HSD(const std::vector<Eigen::ArrayXd>& data);

struct GAMESHOWELL_RESULT {
	std::vector<std::tuple<int, int, double>> r;
};

GAMESHOWELL_RESULT Games_Howell(const std::vector<Eigen::ArrayXd>& data);

