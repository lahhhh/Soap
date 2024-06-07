#pragma once

#include "Identifier.h"

Eigen::ArrayXd evaluate_KDE(const Eigen::ArrayXd& data, const Eigen::ArrayXd& loc);