#include "Monocle3.h"

#include "Custom.h"

Monocle3::Monocle3(const Monocle3& m):
	data_type_(m.data_type_),
	original_embedding_(m.original_embedding_),
	cell_included_(m.cell_included_),
	cell_embedding_(m.cell_embedding_),
	pr_embedding_(m.pr_embedding_),
	pseudo_time_(m.pseudo_time_)
{
	igraph_copy(&this->pr_graph_, &m.pr_graph_);
	igraph_copy(&this->cell_graph_, &m.cell_graph_);
	igraph_vector_init_copy(&this->cell_graph_weights_, &m.cell_graph_weights_);
};

Monocle3& Monocle3::operator=(const Monocle3& rhs) {
	this->data_type_ = rhs.data_type_;
	this->original_embedding_ = rhs.original_embedding_;
	this->cell_included_ = rhs.cell_included_;
	this->cell_embedding_ = rhs.cell_embedding_;
	this->pr_embedding_ = rhs.pr_embedding_;
	this->pseudo_time_ = rhs.pseudo_time_;
	igraph_copy(&this->pr_graph_, &rhs.pr_graph_);
	igraph_copy(&this->cell_graph_, &rhs.cell_graph_);
	custom::igraph_vector_copy(this->cell_graph_weights_, rhs.cell_graph_weights_);

	return *this;
};

Monocle3::~Monocle3() {
	igraph_destroy(&this->pr_graph_);
	igraph_destroy(&this->cell_graph_);

	if (this->cell_graph_weights_.stor_begin != NULL) {
		igraph_vector_destroy(&this->cell_graph_weights_);
	}
};