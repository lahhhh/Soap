#include "RBConfigurationVertexPartition.h"

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
	vector<size_t> const& membership, double resolution_parameter) :
	LinearResolutionParameterVertexPartition(graph,
		membership, resolution_parameter)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
	vector<size_t> const& membership) :
	LinearResolutionParameterVertexPartition(graph,
		membership)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph,
	double resolution_parameter) :
	LinearResolutionParameterVertexPartition(graph, resolution_parameter)
{ }

RBConfigurationVertexPartition::RBConfigurationVertexPartition(Graph* graph) :
	LinearResolutionParameterVertexPartition(graph)
{ }

RBConfigurationVertexPartition::~RBConfigurationVertexPartition()
{ }

RBConfigurationVertexPartition* RBConfigurationVertexPartition::create(Graph* graph)
{
	return new RBConfigurationVertexPartition(graph, this->resolution_parameter);
}

RBConfigurationVertexPartition* RBConfigurationVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
	return new RBConfigurationVertexPartition(graph, membership, this->resolution_parameter);
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
*****************************************************************************/
double RBConfigurationVertexPartition::diff_move(size_t v, size_t new_comm)
{
	size_t old_comm = this->_membership[v];
	double diff = 0.0;
	double total_weight = this->graph->total_weight() * (2.0 - this->graph->is_directed());
	if (total_weight == 0.0)
		return 0.0;
	if (new_comm != old_comm)
	{
		double w_to_old = this->weight_to_comm(v, old_comm);
		double w_from_old = this->weight_from_comm(v, old_comm);
		double w_to_new = this->weight_to_comm(v, new_comm);
		double w_from_new = this->weight_from_comm(v, new_comm);
		double k_out = this->graph->strength(v, IGRAPH_OUT);
		double k_in = this->graph->strength(v, IGRAPH_IN);
		double self_weight = this->graph->node_self_weight(v);
		double K_out_old = this->total_weight_from_comm(old_comm);
		double K_in_old = this->total_weight_to_comm(old_comm);
		double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
		double K_in_new = this->total_weight_to_comm(new_comm) + k_in;
		double diff_old = (w_to_old - this->resolution_parameter * k_out * K_in_old / total_weight) + \
			(w_from_old - this->resolution_parameter * k_in * K_out_old / total_weight);
		double diff_new = (w_to_new + self_weight - this->resolution_parameter * k_out * K_in_new / total_weight) + \
			(w_from_new + self_weight - this->resolution_parameter * k_in * K_out_new / total_weight);
		diff = diff_new - diff_old;
	}
	return diff;
}

/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/
double RBConfigurationVertexPartition::quality(double resolution_parameter)
{
	double mod = 0.0;

	double m;
	if (this->graph->is_directed())
		m = this->graph->total_weight();
	else
		m = 2 * this->graph->total_weight();

	if (m == 0)
		return 0.0;

	for (size_t c = 0; c < this->n_communities(); c++)
	{
		double w = this->total_weight_in_comm(c);
		double w_out = this->total_weight_from_comm(c);
		double w_in = this->total_weight_to_comm(c);
		mod += w - resolution_parameter * w_out * w_in / ((this->graph->is_directed() ? 1.0 : 4.0) * this->graph->total_weight());
	}
	double q = (2.0 - this->graph->is_directed()) * mod;
	return q;
}
