#include "RBERVertexPartition.h"

RBERVertexPartition::RBERVertexPartition(Graph* graph,
	vector<size_t> const& membership, double resolution_parameter) :
	LinearResolutionParameterVertexPartition(graph,
		membership, resolution_parameter)
{ }

RBERVertexPartition::RBERVertexPartition(Graph* graph,
	vector<size_t> const& membership) :
	LinearResolutionParameterVertexPartition(graph,
		membership)
{ }

RBERVertexPartition::RBERVertexPartition(Graph* graph,
	double resolution_parameter) :
	LinearResolutionParameterVertexPartition(graph, resolution_parameter)
{ }

RBERVertexPartition::RBERVertexPartition(Graph* graph) :
	LinearResolutionParameterVertexPartition(graph)
{ }

RBERVertexPartition* RBERVertexPartition::create(Graph* graph)
{
	return new RBERVertexPartition(graph, this->resolution_parameter);
}

RBERVertexPartition* RBERVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
	return new RBERVertexPartition(graph, membership, this->resolution_parameter);
}

RBERVertexPartition::~RBERVertexPartition()
{ }

/********************************************************************************
  RBER implementation of a vertex partition
  (which includes a resolution parameter).
 ********************************************************************************/
double RBERVertexPartition::diff_move(size_t v, size_t new_comm)
{
	size_t old_comm = this->membership(v);
	double diff = 0.0;
	if (new_comm != old_comm)
	{
		double w_to_old = this->weight_to_comm(v, old_comm);
		double w_to_new = this->weight_to_comm(v, new_comm);
		double w_from_old = this->weight_from_comm(v, old_comm);
		double w_from_new = this->weight_from_comm(v, new_comm);
		double nsize = this->graph->node_size(v);
		double csize_old = this->csize(old_comm);
		double csize_new = this->csize(new_comm);
		double self_weight = this->graph->node_self_weight(v);
		double possible_edge_difference_old = 0.0;
		if (this->graph->correct_self_loops())
			possible_edge_difference_old = nsize * (ptrdiff_t)(2.0 * csize_old - nsize);
		else
			possible_edge_difference_old = nsize * (ptrdiff_t)(2.0 * csize_old - nsize - 1.0);
		double diff_old = w_to_old + w_from_old -
			self_weight - this->resolution_parameter * this->graph->density() * possible_edge_difference_old;
		double possible_edge_difference_new = 0.0;
		if (this->graph->correct_self_loops())
			possible_edge_difference_new = nsize * (ptrdiff_t)(2.0 * csize_new + nsize);
		else
			possible_edge_difference_new = nsize * (ptrdiff_t)(2.0 * csize_new + nsize - 1.0);
		double diff_new = w_to_new + w_from_new + self_weight -
			this->resolution_parameter * this->graph->density() * possible_edge_difference_new;
		diff = diff_new - diff_old;
	}
	return diff;
}

double RBERVertexPartition::quality(double resolution_parameter)
{
	double mod = 0.0;
	for (size_t c = 0; c < this->n_communities(); c++)
	{
		double csize = this->csize(c);
		double w = this->total_weight_in_comm(c);
		double comm_possible_edges = this->graph->possible_edges(csize);

		mod += w - resolution_parameter * this->graph->density() * comm_possible_edges;
	}
	return (2.0 - this->graph->is_directed()) * mod;
}
