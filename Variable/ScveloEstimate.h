#pragma once

#include "Identifier.h"

class ScveloEstimate
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(ScveloEstimate);

	template <
		typename GraphType,
		typename GraphNegType,
		typename SelfProbabilityType
		>
	ScveloEstimate(
		GraphType&& graph,
		GraphNegType&& graph_neg,
		SelfProbabilityType&& self_probability
	):
		graph_(std::forward<GraphType>(graph)),
		graph_neg_(std::forward<GraphNegType>(graph_neg)),
		self_probability_(std::forward<SelfProbabilityType>(self_probability))
	{}

	DataType data_type_{ DataType::Plain };

	Eigen::SparseMatrix<double> graph_;
	Eigen::SparseMatrix<double> graph_neg_;

	Eigen::ArrayXd self_probability_;
};

