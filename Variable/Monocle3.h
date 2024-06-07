#pragma once

#include "Identifier.h"

#include <igraph.h>

#include "Embedding.h"

class Monocle3
{
public:
	enum class DataType : int { Plain = 0 };

	Monocle3() = default;
	Monocle3(const Monocle3&);
	Monocle3& operator=(const Monocle3&);
	~Monocle3();

	DataType data_type_;

	igraph_t pr_graph_;
	igraph_t cell_graph_;

	Embedding original_embedding_;
	Eigen::ArrayX<bool> cell_included_;

	Eigen::MatrixXd cell_embedding_;
	Eigen::MatrixXd pr_embedding_;

	igraph_vector_t cell_graph_weights_;

	Eigen::ArrayXd pseudo_time_;

	G_SET_IDENTIFIER("Monocle3");
};

