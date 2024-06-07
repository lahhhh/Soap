#include "Network.h"

#include "Custom.h"

Network::Network(int n_nodes, const std::vector< std::vector<int> >& edge)
	: Network(n_nodes, nullptr, edge, nullptr)
{}

Network::Network(
	int n_nodes,
	std::vector<double>* node_weight,
	const std::vector< std::vector<int> >& edge
)
	: Network(n_nodes, node_weight, edge, nullptr)
{}

Network::Network(
	int n_nodes,
	const std::vector< std::vector<int> >& edge,
	std::vector<double>* edge_weight
)
	: Network(n_nodes, nullptr, edge, edge_weight)
{}

Network::Network(
	int n_nodes,
	std::vector<double>* node_weight,
	const std::vector< std::vector<int> >& edge,
	std::vector<double>* edge_weight
) :
	n_nodes_(n_nodes),
	n_edges_(0),
	total_edge_weight_self_links_(0.0),
	first_neighbor_index_(n_nodes + 1, 0)
{
	const std::size_t edge_length = edge[0].size();
	std::vector<double> edge_weight2(edge_length, 0.0);

	std::vector<int> neighbor(edge_length, 0);

	std::size_t i = 1;
	for (std::size_t j = 0; j < edge_length; ++j) {

		if (edge[0][j] != edge[1][j])
		{
			if (edge[0][j] >= i) {
				for (; i <= edge[0][j]; ++i) {
					this->first_neighbor_index_[i] = this->n_edges_;
				}
			}
			neighbor[this->n_edges_] = edge[1][j];
			edge_weight2[this->n_edges_] = (edge_weight != nullptr) ? (*edge_weight)[j] : 1;
			++this->n_edges_;
		}
		else {
			this->total_edge_weight_self_links_ += (edge_weight != nullptr) ? (*edge_weight)[j] : 1;
		}
	}
	for (; i <= n_nodes; ++i) {
		this->first_neighbor_index_[i] = this->n_edges_;
	}

	this->neighbor_.resize(this->n_edges_);
	std::copy(neighbor.begin(), neighbor.begin() + this->n_edges_, this->neighbor_.begin());

	this->edge_weight_.resize(this->n_edges_);
	std::copy(edge_weight2.begin(), edge_weight2.begin() + this->n_edges_, this->edge_weight_.begin());

	this->node_weight_ = (node_weight != nullptr) ? *node_weight : this->get_total_edge_weight_per_node();
}

Network::Network(
	int n_nodes,
	std::vector<int>* first_neighbor_index,
	std::vector<int>* neighbor
)
	: Network(n_nodes, nullptr, first_neighbor_index, neighbor, nullptr)
{}

Network::Network(
	int n_nodes,
	std::vector<double>* node_weight,
	std::vector<int>* first_neighbor_index,
	std::vector<int>* neighbor
)
	: Network(n_nodes, node_weight, first_neighbor_index, neighbor, nullptr)
{}

Network::Network(
	int n_nodes,
	std::vector<int>* first_neighbor_index,
	std::vector<int>* neighbor,
	std::vector<double>* edge_weight
)
	: Network(n_nodes, nullptr, first_neighbor_index, neighbor, edge_weight)
{}

void Network::set_from(int n_nodes,
	const std::vector<int>& first_neighbor_index,
	const std::vector<int>& neighbor,
	const std::vector<double>& edge_weight
)
{
	this->n_nodes_ = n_nodes;
	this->n_edges_ = neighbor.size();
	this->total_edge_weight_self_links_ = 0;
	this->first_neighbor_index_ = first_neighbor_index;
	this->neighbor_ = neighbor;
	this->edge_weight_ = edge_weight;
	this->node_weight_ = this->get_total_edge_weight_per_node();
};

Network::Network(
	int n_nodes,
	std::vector<double>* node_weight,
	std::vector<int>* first_neighbor_index,
	std::vector<int>* neighbor,
	std::vector<double>* edge_weight
) :
	n_nodes_(n_nodes),
	n_edges_(neighbor->size()),
	total_edge_weight_self_links_(0.0),
	first_neighbor_index_(*first_neighbor_index),
	neighbor_(*neighbor)
{
	if (edge_weight != nullptr)
		this->edge_weight_ = *edge_weight;
	else
	{
		this->edge_weight_.resize(this->n_edges_, 1.0);
	}

	this->node_weight_ = (node_weight != nullptr) ? *node_weight : get_total_edge_weight_per_node();
}

void Network::set_from(int n_nodes,
	const std::vector<double>& node_weight,
	const std::vector<int>& first_neighbor_index,
	const std::vector<int>& neighbor,
	const std::vector<double>& edge_weight
)
{
	this->n_nodes_ = n_nodes;
	this->n_edges_ = neighbor.size();
	this->total_edge_weight_self_links_ = 0;
	this->first_neighbor_index_ = first_neighbor_index;
	this->neighbor_ = neighbor;
	this->edge_weight_ = edge_weight;
	this->node_weight_ = node_weight;
};

Network::Network(
	int n_nodes,
	int n_edges,
	double total_edge_weight_self_links,
	const std::vector<double>& node_weight,
	const std::vector<double>& edge_weight,
	const std::vector<int>& first_neighbor_index,
	const std::vector<int>& neighbor
) :
	n_nodes_(n_nodes),
	n_edges_(n_edges_),
	total_edge_weight_self_links_(total_edge_weight_self_links),
	node_weight_(node_weight),
	edge_weight_(edge_weight),
	first_neighbor_index_(first_neighbor_index),
	neighbor_(neighbor)
{}


int Network::get_n_nodes() const
{
	return this->n_nodes_;
}

double Network::get_total_node_weight() const
{
	return _Cs sum(this->node_weight_);
}

std::vector<double> Network::get_node_weights() const
{
	return this->node_weight_;
}

double Network::get_node_weight(int node) const
{
	return this->node_weight_[node];
}

int Network::get_n_edges() const
{
	return this->n_edges_ / 2;
}

int Network::get_n_edges(int node) const
{
	return this->first_neighbor_index_[node + 1] - this->first_neighbor_index_[node];
}

std::vector<int> Network::get_n_edges_per_node() const
{
	std::vector<int> n_edges_per_node(this->n_nodes_, 0);

	for (int i = 0; i < this->n_nodes_; ++i)
		n_edges_per_node[i] = this->first_neighbor_index_[i + 1] - this->first_neighbor_index_[i];
	return n_edges_per_node;
}

std::vector< std::vector<int> > Network::get_edges() const
{
	std::vector< std::vector<int> > edge(2);

	edge[0].resize(this->n_edges_, 0);
	for (std::size_t i = 0; i < this->n_nodes_; ++i) {
		std::fill(edge[0].begin() + this->first_neighbor_index_[i], edge[0].begin() + this->first_neighbor_index_[i + 1], i);
	}
	edge[1] = this->neighbor_;
	return edge;
}

std::vector<int> Network::get_edges(int node) const
{
	return std::vector<int>(this->neighbor_.cbegin() + this->first_neighbor_index_[node],
		this->neighbor_.cbegin() + this->first_neighbor_index_[node + 1]);
}

std::vector< std::vector<int> > Network::get_edges_per_node() const
{
	std::vector< std::vector<int> > edges_per_node(this->n_nodes_);

	for (std::size_t i = 0; i < this->n_nodes_; ++i) {
		edges_per_node[i] = std::vector<int>(this->neighbor_.begin() + this->first_neighbor_index_[i],
			this->neighbor_.begin() + this->first_neighbor_index_[i + 1]);
	}
	return edges_per_node;
}

double Network::get_total_edge_weight() const
{

	return _Cs sum(this->edge_weight_) / 2;
}

double Network::get_total_edge_weight(int node) const
{
	return std::accumulate(this->edge_weight_.cbegin() + this->first_neighbor_index_[node],
		this->edge_weight_.cbegin() + this->first_neighbor_index_[node + 1],
		0.0);
}

std::vector<double> Network::get_total_edge_weight_per_node() const
{
	std::vector<double> total_edge_weight_per_node(this->n_nodes_);

	for (std::size_t i = 0; i < this->n_nodes_; ++i) {
		total_edge_weight_per_node[i] = std::accumulate(this->edge_weight_.cbegin() + this->first_neighbor_index_[i],
			this->edge_weight_.cbegin() + this->first_neighbor_index_[i + 1],
			0.0);
	}
	return total_edge_weight_per_node;
}

std::vector<double> Network::get_edge_weights() const
{
	return this->edge_weight_;
}

std::vector<double> Network::get_edge_weights(int node) const
{
	return std::vector<double>(this->edge_weight_.cbegin() + this->first_neighbor_index_[node],
		this->edge_weight_.cbegin() + this->first_neighbor_index_[node + 1]);
}

std::vector< std::vector<double> > Network::get_edge_weights_per_node() const
{
	std::vector< std::vector<double> > edge_weight_per_node(this->n_nodes_);

	for (std::size_t i = 0; i < this->n_nodes_; ++i) {
		edge_weight_per_node[i] = std::vector<double>(this->edge_weight_.cbegin() + this->first_neighbor_index_[i],
			this->edge_weight_.cbegin() + this->first_neighbor_index_[i + 1]);
	}
	return edge_weight_per_node;
}

double Network::get_total_edge_weight_self_links() const
{
	return this->total_edge_weight_self_links_;
}

Network Network::create_network_without_node_weights() const
{

	return Network(this->n_nodes_,
		this->n_edges_,
		this->total_edge_weight_self_links_,
		std::vector<double>(this->n_nodes_, 1),
		this->edge_weight_,
		this->first_neighbor_index_,
		this->neighbor_);
}

Network Network::create_network_without_edge_weights() const
{
	return Network(this->n_nodes_,
		this->n_edges_,
		0,
		this->node_weight_,
		std::vector<double>(this->n_edges_, 1),
		this->first_neighbor_index_,
		this->neighbor_);
}

Network Network::create_network_without_node_and_edge_weights() const
{
	return Network(this->n_nodes_,
		this->n_edges_,
		0,
		std::vector<double>(this->n_nodes_, 1),
		std::vector<double>(this->n_edges_, 1),
		this->first_neighbor_index_,
		this->neighbor_);
}

Network Network::create_normalized_network1() const
{
	Network normalized_network(this->n_nodes_,
		this->n_edges_,
		0,
		std::vector<double>(this->n_nodes_, 1),
		std::vector<double>(this->n_edges_, 0),
		this->first_neighbor_index_,
		this->neighbor_
	);

	double total_node_weight = get_total_node_weight();
	for (int i = 0; i < this->n_nodes_; ++i) {
		for (int j = this->first_neighbor_index_[i]; j < this->first_neighbor_index_[i + 1]; ++j) {
			normalized_network.edge_weight_[j] = this->edge_weight_[j] / ((this->node_weight_[i] * this->node_weight_[this->neighbor_[j]]) / total_node_weight);
		}
	}

	normalized_network.total_edge_weight_self_links_ = 0;

	return normalized_network;
}

Network Network::create_normalized_network2() const
{
	Network normalized_network(this->n_nodes_,
		this->n_edges_,
		0,
		std::vector<double>(this->n_nodes_, 1),
		std::vector<double>(this->n_edges_, 0),
		this->first_neighbor_index_,
		this->neighbor_
	);

	for (int i = 0; i < this->n_nodes_; ++i) {
		for (int j = this->first_neighbor_index_[i]; j < this->first_neighbor_index_[i + 1]; ++j) {
			normalized_network.edge_weight_[j] = this->edge_weight_[j] / (2 / (this->n_nodes_ / this->node_weight_[i] + this->n_nodes_ / this->node_weight_[this->neighbor_[j]]));
		}
	}

	return normalized_network;
}

Network Network::create_pruned_network(int n_edges) const
{
	std::default_random_engine re;
	re.seed(1997);
	return create_pruned_network(n_edges, re);
}

Network Network::create_pruned_network(int n_edges, std::default_random_engine& re) const
{

	int i = 0;

	n_edges *= 2;

	if (n_edges >= this->n_edges_)
		return *this;

	std::vector<double> edge_weight(this->n_edges_ / 2);
	for (int j = 0; j < this->n_nodes_; ++j) {
		for (int k = this->first_neighbor_index_[j]; k < this->first_neighbor_index_[j + 1]; ++k) {
			if (this->neighbor_[k] < j)
			{
				edge_weight[i] = this->edge_weight_[k];
				++i;
			}
		}
	}

	std::ranges::sort(edge_weight);

	double edge_weight_threshold = edge_weight[(this->n_edges_ - n_edges) / 2];

	int n_edges_above_threshold = 0;
	while (edge_weight[this->n_edges_ / 2 - n_edges_above_threshold - 1] > edge_weight_threshold) {
		n_edges_above_threshold++;
	}
	int n_edgesAt_threshold = 0;
	while ((n_edges_above_threshold + n_edgesAt_threshold < this->n_edges_ / 2) && (edge_weight[this->n_edges_ / 2 - n_edges_above_threshold - n_edgesAt_threshold - 1] == edge_weight_threshold))
		n_edgesAt_threshold++;

	std::vector<int> node_permutation = _Cs generate_random_permutation(this->n_nodes_, re);

	std::vector<double> random_number(n_edgesAt_threshold, 0);
	i = 0;
	for (int j = 0; j < this->n_nodes_; ++j) {
		for (int k = this->first_neighbor_index_[j]; k < this->first_neighbor_index_[j + 1]; ++k) {
			if ((this->neighbor_[k] < j) && (this->edge_weight_[k] == edge_weight_threshold))
			{
				random_number[i] = generate_random_number(j, this->neighbor_[k], node_permutation);
				++i;
			}
		}
	}
	std::ranges::sort(random_number);

	double random_number_threshold = random_number[n_edges_above_threshold + n_edgesAt_threshold - n_edges / 2];

	Network pruned_network(
		this->n_nodes_,
		n_edges,
		this->total_edge_weight_self_links_,
		this->node_weight_,
		std::vector<double>(n_edges, 0),
		std::vector<int>(this->n_nodes_ + 1, 0),
		std::vector<int>(n_edges, 0)
	);

	pruned_network.edge_weight_.resize(n_edges);
	i = 0;
	for (int j = 0; j < this->n_nodes_; ++j)
	{
		for (int k = this->first_neighbor_index_[j]; k < this->first_neighbor_index_[j + 1]; ++k) {
			if ((this->edge_weight_[k] > edge_weight_threshold) || ((this->edge_weight_[k] == edge_weight_threshold) && (generate_random_number(j, this->neighbor_[k], node_permutation) >= random_number_threshold)))
			{
				pruned_network.neighbor_[i] = this->neighbor_[k];
				pruned_network.edge_weight_[i] = this->edge_weight_[k];
				++i;
			}
		}
		pruned_network.first_neighbor_index_[j + 1] = i;
	}

	return pruned_network;
}

Network Network::create_subnetwork(const std::vector<int>& node) const
{
	Network subnetwork;

	subnetwork.n_nodes_ = node.size();

	if (subnetwork.n_nodes_ == 1)
	{
		subnetwork.n_edges_ = 0;
		subnetwork.node_weight_.resize(1, this->node_weight_[node[0]]);
		subnetwork.first_neighbor_index_.resize(2, 0);
	}
	else
	{
		std::vector<int> subnetwork_node(this->n_nodes_, -1);
		for (int i = 0; i < node.size(); ++i) {
			subnetwork_node[node[i]] = i;
		}

		subnetwork.n_edges_ = 0;
		subnetwork.node_weight_.resize(subnetwork.n_nodes_, 0);
		subnetwork.first_neighbor_index_.resize(subnetwork.n_nodes_ + 1);
		std::vector<int> subnetwork_neighbor(this->n_edges_, 0);
		std::vector<double> subnetwork_edge_weight(this->n_edges_, 0.0);
		for (int i = 0; i < subnetwork.n_nodes_; ++i)
		{
			int j = node[i];
			subnetwork.node_weight_[i] = this->node_weight_[j];
			for (int k = this->first_neighbor_index_[j]; k < this->first_neighbor_index_[j + 1]; ++k)
				if (subnetwork_node[this->neighbor_[k]] >= 0)
				{
					subnetwork_neighbor[subnetwork.n_edges_] = subnetwork_node[this->neighbor_[k]];
					subnetwork_edge_weight[subnetwork.n_edges_] = this->edge_weight_[k];
					++subnetwork.n_edges_;
				}
			subnetwork.first_neighbor_index_[i + 1] = subnetwork.n_edges_;
		}

		subnetwork.neighbor_ = std::vector<int>(subnetwork_neighbor.cbegin(), subnetwork_neighbor.cbegin() + subnetwork.n_edges_);
		subnetwork.edge_weight_ = std::vector<double>(subnetwork_edge_weight.cbegin(), subnetwork_edge_weight.cbegin() + subnetwork.n_edges_);
	}

	subnetwork.total_edge_weight_self_links_ = 0;

	return subnetwork;
}

Network Network::create_subnetwork(const QVector<bool>& node_in_subnetwork) const
{
	int i = 0;

	for (int j = 0; j < this->n_nodes_; ++j) {
		if (node_in_subnetwork[j]) {
			++i;
		}
	}

	std::vector<int> node(i, 0);
	i = 0;
	for (int j = 0; j < this->n_nodes_; ++j) {
		if (node_in_subnetwork[j])
		{
			node[i] = j;
			++i;
		}
	}
	return create_subnetwork(node);
}

Network Network::create_subnetwork(const Clustering& clustering, int cluster) const
{

	std::vector< std::vector<int> > node_per_cluster = clustering.get_nodes_per_cluster();
	std::vector<int> subnetwork_node(this->n_nodes_, 0);
	std::vector<int> subnetwork_neighbor(this->n_edges_, 0);
	std::vector<double> subnetwork_edge_weight(this->n_edges_, 0);

	Network subnetwork = create_subnetwork(clustering, cluster, node_per_cluster[cluster], subnetwork_node, subnetwork_neighbor, subnetwork_edge_weight);

	return subnetwork;
}

std::vector<Network> Network::create_subnetworks(const Clustering& clustering) const
{
	std::vector<double> subnetwork_edge_weight(this->n_edges_, 0.0);
	std::vector<int> subnetwork_neighbor(this->n_edges_, 0), subnetwork_node(this->n_nodes_, 0);
	std::vector< std::vector<int> > node_per_cluster = clustering.get_nodes_per_cluster();;
	std::vector<Network> subnetwork(clustering.n_clusters_);

	for (int i = 0; i < clustering.n_clusters_; ++i) {
		subnetwork[i] = create_subnetwork(clustering, i, node_per_cluster[i], subnetwork_node, subnetwork_neighbor, subnetwork_edge_weight);
	}
	return subnetwork;
}

Network Network::create_subnetwork_largest_component() const
{
	return create_subnetwork(identify_components(), 0);
}

Network Network::create_reduced_network(const Clustering& clustering) const
{
	Network reduced_network;

	reduced_network.n_nodes_ = clustering.n_clusters_;

	reduced_network.n_edges_ = 0;
	reduced_network.node_weight_.resize(clustering.n_clusters_, 0);
	reduced_network.first_neighbor_index_.resize(clustering.n_clusters_ + 1, 0);
	reduced_network.total_edge_weight_self_links_ = this->total_edge_weight_self_links_;

	std::vector<int> reduced_network_neighbor1(this->n_edges_, 0);
	std::vector<double> reduced_network_edge_weight1(this->n_edges_, 0);

	std::vector<int> reduced_network_neighbor2(clustering.n_clusters_ - 1, 0);
	std::vector<double> reduced_network_edge_weight2(clustering.n_clusters_, 0);
	std::vector< std::vector<int> > node_per_cluster = clustering.get_nodes_per_cluster();
	for (int i = 0; i < clustering.n_clusters_; ++i)
	{
		int j = 0;
		for (int k = 0; k < node_per_cluster[i].size(); ++k)
		{
			int l = node_per_cluster[i][k];

			reduced_network.node_weight_[i] += this->node_weight_[l];

			for (int m = first_neighbor_index_[l]; m < first_neighbor_index_[l + 1]; ++m)
			{
				int n = clustering.cluster_[this->neighbor_[m]];
				if (n != i)
				{
					if (reduced_network_edge_weight2[n] == 0)
					{
						reduced_network_neighbor2[j] = n;
						++j;
					}
					reduced_network_edge_weight2[n] += this->edge_weight_[m];
				}
				else
					reduced_network.total_edge_weight_self_links_ += this->edge_weight_[m];
			}
		}

		for (int k = 0; k < j; ++k)
		{
			reduced_network_neighbor1[reduced_network.n_edges_ + k] = reduced_network_neighbor2[k];
			reduced_network_edge_weight1[reduced_network.n_edges_ + k] = reduced_network_edge_weight2[reduced_network_neighbor2[k]];
			reduced_network_edge_weight2[reduced_network_neighbor2[k]] = 0;
		}
		reduced_network.n_edges_ += j;
		reduced_network.first_neighbor_index_[i + 1] = reduced_network.n_edges_;
	}
	reduced_network.neighbor_ = std::vector<int>(reduced_network_neighbor1.cbegin(), reduced_network_neighbor1.cbegin() + reduced_network.n_edges_);
	reduced_network.edge_weight_ = std::vector<double>(reduced_network_edge_weight1.cbegin(), reduced_network_edge_weight1.cbegin() + reduced_network.n_edges_);

	return reduced_network;
}

Clustering Network::identify_components() const
{
	QVector<bool> node_visited(this->n_nodes_);
	Clustering clustering(this->n_nodes_);
	std::vector<int> node(this->n_nodes_, 0);

	clustering.n_clusters_ = 0;

	for (int i = 0; i < this->n_nodes_; ++i)
		if (!node_visited[i])
		{
			clustering.cluster_[i] = clustering.n_clusters_;
			node_visited[i] = true;
			node[0] = i;
			int j = 1, k = 0;
			do
			{
				for (int l = this->first_neighbor_index_[node[k]]; l < this->first_neighbor_index_[node[k] + 1]; l++)
					if (!node_visited[this->neighbor_[l]])
					{
						clustering.cluster_[this->neighbor_[l]] = clustering.n_clusters_;
						node_visited[this->neighbor_[l]] = true;
						node[j] = this->neighbor_[l];
						++j;
					}
				++k;
			} while (k < j);

			clustering.n_clusters_++;
		}

	clustering.order_clusters_by_n_nodes();

	return clustering;
}

double Network::generate_random_number(int node1, int node2, const std::vector<int>& node_permutation) const
{
	int i, j;
	std::default_random_engine e;
	std::uniform_real_distribution<double> u(0, 1);


	if (node1 < node2)
	{
		i = node1;
		j = node2;
	}
	else
	{
		i = node2;
		j = node1;
	}
	e.seed(node_permutation[i] * this->n_nodes_ + node_permutation[j]);

	return u(e);
}

Network Network::create_subnetwork(
	const Clustering& clustering,
	int cluster,
	const std::vector<int>& node,
	std::vector<int> subnetwork_node,
	std::vector<int> subnetwork_neighbor,
	std::vector<double> subnetwork_edge_weight)
	const
{
	Network subnetwork;

	subnetwork.n_nodes_ = node.size();

	if (subnetwork.n_nodes_ == 1)
	{
		subnetwork.n_edges_ = 0;
		subnetwork.node_weight_.resize(1, this->node_weight_[node[0]]);
		subnetwork.first_neighbor_index_.resize(2, 0);
	}
	else
	{
		for (int i = 0; i < node.size(); ++i) {
			subnetwork_node[node[i]] = i;
		}

		subnetwork.n_edges_ = 0;
		subnetwork.node_weight_.resize(subnetwork.n_nodes_, 0);
		subnetwork.first_neighbor_index_.resize(subnetwork.n_nodes_ + 1, 0);
		for (int i = 0; i < subnetwork.n_nodes_; ++i)
		{
			int j = node[i];
			subnetwork.node_weight_[i] = this->node_weight_[j];
			for (int k = this->first_neighbor_index_[j]; k < this->first_neighbor_index_[j + 1]; ++k)
				if (clustering.cluster_[this->neighbor_[k]] == cluster)
				{
					subnetwork_neighbor[subnetwork.n_edges_] = subnetwork_node[this->neighbor_[k]];
					subnetwork_edge_weight[subnetwork.n_edges_] = this->edge_weight_[k];
					++subnetwork.n_edges_;
				}
			subnetwork.first_neighbor_index_[i + 1] = subnetwork.n_edges_;
		}
		subnetwork.neighbor_ = std::vector<int>(subnetwork_neighbor.cbegin(), subnetwork_neighbor.cbegin() + subnetwork.n_edges_);
		subnetwork.edge_weight_ = std::vector<double>(subnetwork_edge_weight.cbegin(), subnetwork_edge_weight.cbegin() + subnetwork.n_edges_);
	}

	subnetwork.total_edge_weight_self_links_ = 0;

	return subnetwork;
}

