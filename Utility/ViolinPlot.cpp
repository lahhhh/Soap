#include "ViolinPlot.h"

#include <cmath>
#include <numbers>
#include <algorithm>
#include <numeric>
#include <vector>

double gaussian_kernel(double x) {
    return 1.0 / sqrt(2 * std::numbers::pi) * exp(-0.5 * x * x);
}

double calculate_bandwidth(const Eigen::ArrayXd& data) {
    double n = data.size();
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / n - mean * mean);
    return pow(4 * pow(stdev, 5) / (3 * n), 0.2);
}

double KDE(const Eigen::ArrayXd& data, double x, double bandwidth) {
    double sum = 0.0;
    for (double xi : data) {
        sum += gaussian_kernel((x - xi) / bandwidth);
    }
    return sum / (data.size() * bandwidth);
}

Eigen::ArrayXd evaluate_KDE(const Eigen::ArrayXd& data, const Eigen::ArrayXd& loc) {

    const int n_point = loc.size();
    Eigen::ArrayXd result(n_point);
    double bandwidth = calculate_bandwidth(data);
    for (int i = 0; i < n_point; ++i) {
        result[i] = KDE(data, loc[i], bandwidth);
    }
    return result;
}