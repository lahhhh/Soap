#include "SignificanceVertexPartition.h"


SignificanceVertexPartition::SignificanceVertexPartition(Graph* graph,
	vector<size_t> const& membership) :
	MutableVertexPartition(graph,
		membership)
{ }

SignificanceVertexPartition::SignificanceVertexPartition(Graph* graph) :
	MutableVertexPartition(graph)
{ }

SignificanceVertexPartition* SignificanceVertexPartition::create(Graph* graph)
{
	return new SignificanceVertexPartition(graph);
}

SignificanceVertexPartition* SignificanceVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
	return new SignificanceVertexPartition(graph, membership);
}

SignificanceVertexPartition::~SignificanceVertexPartition()
{ }

double SignificanceVertexPartition::diff_move(size_t v, size_t new_comm)
{
	size_t old_comm = this->membership(v);
	double nsize = this->graph->node_size(v);
	double diff = 0.0;
	if (new_comm != old_comm)
	{
		double normalise = (2.0 - this->graph->is_directed());
		double p = this->graph->density();

		//Old comm
		double n_old = this->csize(old_comm);
		double N_old = this->graph->possible_edges(n_old);
		double m_old = this->total_weight_in_comm(old_comm);
		double q_old = 0.0;
		if (N_old > 0)
			q_old = m_old / N_old;
		// Old comm after move
		double n_oldx = n_old - nsize; // It should not be possible that this becomes negative, so no need for ptrdiff_t here.
		double N_oldx = this->graph->possible_edges(n_oldx);
		double sw = this->graph->node_self_weight(v);
		// Be careful to exclude the self weight here, because this is include in the weight_to_comm function.
		double wtc = this->weight_to_comm(v, old_comm) - sw;
		double wfc = this->weight_from_comm(v, old_comm) - sw;
		double m_oldx = m_old - wtc / normalise - wfc / normalise - sw;
		double q_oldx = 0.0;
		if (N_oldx > 0)
			q_oldx = m_oldx / N_oldx;

		// New comm
		double n_new = this->csize(new_comm);
		double N_new = this->graph->possible_edges(n_new);
		double m_new = this->total_weight_in_comm(new_comm);
		double q_new = 0.0;
		if (N_new > 0)
			q_new = m_new / N_new;

		// New comm after move
		double n_newx = n_new + nsize;
		double N_newx = this->graph->possible_edges(n_newx);
		wtc = this->weight_to_comm(v, new_comm);
		wfc = this->weight_from_comm(v, new_comm);
		sw = this->graph->node_self_weight(v);
		double m_newx = m_new + wtc / normalise + wfc / normalise + sw;
		double q_newx = 0.0;
		if (N_newx > 0)
			q_newx = m_newx / N_newx;

		// Calculate actual diff

		if (N_oldx != N_new || q_oldx != q_new)
			diff += (double)N_oldx * KLL(q_oldx, p) - (double)N_new * KLL(q_new, p);

		if (N_newx != N_old || q_newx != q_old)
			diff += (double)N_newx * KLL(q_newx, p) - (double)N_old * KLL(q_old, p);

	}
	return diff;
}

/********************************************************************************
   Calculate the significance of the partition.
*********************************************************************************/
double SignificanceVertexPartition::quality()
{
	double S = 0.0;
	double p = this->graph->density();
	for (size_t c = 0; c < this->n_communities(); c++)
	{
		double n_c = this->csize(c);
		double m_c = this->total_weight_in_comm(c);
		double p_c = 0.0;
		size_t N_c = this->graph->possible_edges(n_c);
		if (N_c > 0)
			p_c = m_c / N_c;
		S += N_c * KLL(p_c, p);
	}
	return S;
}
