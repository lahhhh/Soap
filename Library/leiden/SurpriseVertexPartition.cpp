#include "SurpriseVertexPartition.h"

SurpriseVertexPartition::SurpriseVertexPartition(Graph* graph,
	vector<size_t> const& membership) :
	MutableVertexPartition(graph,
		membership)
{ }

SurpriseVertexPartition::SurpriseVertexPartition(Graph* graph) :
	MutableVertexPartition(graph)
{ }

SurpriseVertexPartition* SurpriseVertexPartition::create(Graph* graph)
{
	return new SurpriseVertexPartition(graph);
}

SurpriseVertexPartition* SurpriseVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
	return new  SurpriseVertexPartition(graph, membership);
}

SurpriseVertexPartition::~SurpriseVertexPartition()
{ }

double SurpriseVertexPartition::diff_move(size_t v, size_t new_comm)
{
	size_t old_comm = this->membership(v);
	size_t nsize = this->graph->node_size(v);
	double diff = 0.0;
	double m = this->graph->total_weight();

	if (m == 0)
		return 0.0;

	if (new_comm != old_comm)
	{
		double normalise = (2.0 - this->graph->is_directed());
		size_t n = this->graph->total_size();
		size_t n2 = this->graph->possible_edges(n);


		// Before move
		double mc = this->total_weight_in_all_comms();
		size_t nc2 = this->total_possible_edges_in_all_comms();

		// To old comm
		size_t n_old = this->csize(old_comm);
		double sw = this->graph->node_self_weight(v);
		double wtc = this->weight_to_comm(v, old_comm) - sw;
		double wfc = this->weight_from_comm(v, old_comm) - sw;
		double m_old = wtc / normalise + wfc / normalise + sw;

		// To new comm
		size_t n_new = this->csize(new_comm);
		wtc = this->weight_to_comm(v, new_comm);
		wfc = this->weight_from_comm(v, new_comm);
		sw = this->graph->node_self_weight(v);
		double m_new = wtc / normalise + wfc / normalise + sw;

		double q = mc / m;
		double s = (double)nc2 / (double)n2;
		double q_new = (mc - m_old + m_new) / m;
		double delta_nc2 = 2.0 * nsize * (ptrdiff_t)(n_new - n_old + nsize) / normalise;
		double s_new = (double)(nc2 + delta_nc2) / (double)n2;
		diff = m * (KLL(q_new, s_new) - KLL(q, s));
	}
	return diff;
}

double SurpriseVertexPartition::quality()
{

	double mc = this->total_weight_in_all_comms();
	size_t nc2 = this->total_possible_edges_in_all_comms();
	double m = this->graph->total_weight();
	size_t n = this->graph->total_size();

	if (m == 0)
		return 0.0;

	size_t n2 = this->graph->possible_edges(n);

	double q = mc / m;
	double s = (double)nc2 / (double)n2;
	double S = m * KLL(q, s);
	return S;
}
