#pragma once

#include "Custom.h"

#include <boost\math\distributions\logistic.hpp>
#include <boost\math\distributions\students_t.hpp>
#include <boost\math\distributions\fisher_f.hpp>
#include <boost\math\distributions\chi_squared.hpp>

double p_wilcox(double q, double m, double n, bool lower_tail = true);

double p_chisq(double q, int df, bool lower_tail = true);

double p_normal(double x, double mu = 0.0, double sigma = 1.0, bool lower_tail = true);

double dnorm(double x, double mu = 0.0, double sigma = 1.0);

Eigen::MatrixXd dnorm(const Eigen::MatrixXd& mat, double mu = 0.0, double sigma = 1.0);

// n > 0, x - finite, n - finite
double p_students_t(double x, double n, bool lower_tail = true);

double p_logistic(double q, double location = 0.0, double scale = 1.0, bool lower_tail = true);

double p_fisher_f(double f, double m, double n, bool lower_tail = true);

template <typename ContainerType>
	requires custom::is_specific_container<ContainerType, double>
ContainerType p_wilcox(const ContainerType& con, double m, double n, bool lower_tail = true) {

	const auto size = std::size(con);
	ContainerType ret(size);

	auto iter = std::cbegin(con);
	auto end = std::cend(con);

	auto to = std::begin(ret);

	for (; iter != end; ++iter, ++to) {
		*to = p_wilcox(*iter, m, n, lower_tail);
	}

	return ret;
}

template <typename ContainerType>
	requires custom::is_specific_container<ContainerType, double>
ContainerType p_normal(const ContainerType& con, double mu = 0.0, double sigma = 1.0, bool lower_tail = true) {

	boost::math::normal distribution(mu, sigma);

	const auto size = std::size(con);
	ContainerType ret(size);

	auto iter = std::cbegin(con);
	auto end = std::cend(con);

	auto to = std::begin(ret);

	for (; iter != end; ++iter, ++to) {
		double p = boost::math::cdf(distribution, *iter);
		*to = lower_tail ? p : (1.0 - p);
	}

	return ret;
}

template <typename ContainerType>
	requires custom::is_specific_container<ContainerType, double>
ContainerType p_students_t(const ContainerType& con, double n, bool lower_tail = true) {

	boost::math::students_t distribution(n);

	const auto size = std::size(con);
	ContainerType ret(size);

	auto iter = std::cbegin(con);
	auto end = std::cend(con);

	auto to = std::begin(ret);

	for (; iter != end; ++iter, ++to) {
		double p = boost::math::cdf(distribution, *iter);
		*to = lower_tail ? p : (1.0 - p);
	}

	return ret;
}

template <typename ContainerType>
	requires custom::is_specific_container<ContainerType, double>
ContainerType p_fisher_f(const ContainerType& con, double m, double n, bool lower_tail = true) {

	boost::math::fisher_f distribution(m, n);

	const auto size = std::size(con);
	ContainerType ret(size);

	auto iter = std::cbegin(con);
	auto end = std::cend(con);

	auto to = std::begin(ret);

	for (; iter != end; ++iter, ++to) {
		double p = boost::math::cdf(distribution, *iter);
		*to = lower_tail ? p : (1.0 - p);
	}

	return ret;
}

template <typename ContainerType>
requires custom::is_specific_container<ContainerType, double>
ContainerType p_logistic(const ContainerType& con, double location = 0.0, double scale = 1.0, bool lower_tail = true) {

	boost::math::logistic distribution(location, scale);

	const auto size = std::size(con);
	ContainerType ret(size);

	auto iter = std::cbegin(con);
	auto end = std::cend(con);

	auto to = std::begin(ret);

	for (; iter != end; ++iter, ++to) {
		double p = boost::math::cdf(distribution, *iter);
		*to = lower_tail ? p : (1.0 - p);
	}

	return ret;
}