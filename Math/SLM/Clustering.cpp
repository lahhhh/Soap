#include "Clustering.h"

#include <algorithm>

#include "custom.h"

Clustering::Clustering(int n_nodes) :
	n_nodes_(n_nodes),
	n_clusters_(1),
	cluster_(n_nodes, 0)
{}

Clustering::Clustering(const std::vector<int>& cluster) :
	n_nodes_(cluster.size()),
	n_clusters_(std::ranges::max(cluster) + 1),
	cluster_(cluster)
{}


int Clustering::get_n_nodes() const
{
	return this->n_nodes_;
}

int Clustering::get_n_clusters() const
{
	return this->n_clusters_;
}

std::vector<int> Clustering::get_clusters() const
{
	return this->cluster_;
}

int Clustering::get_cluster(int node) const
{
	return this->cluster_[node];
}

std::vector<int> Clustering::get_n_nodes_per_cluster() const
{
	std::vector<int> n_nodes_per_cluster(this->n_clusters_, 0);

	for (int i = 0; i < this->n_nodes_; ++i) {
		++n_nodes_per_cluster[this->cluster_[i]];
	}

	return n_nodes_per_cluster;
}

std::vector< std::vector<int> > Clustering::get_nodes_per_cluster() const
{
	std::vector<int> n_nodes_per_cluster = get_n_nodes_per_cluster();
	std::vector< std::vector<int> > node_per_cluster(this->n_clusters_);

	for (int i = 0; i < this->n_clusters_; ++i)
	{
		node_per_cluster[i].reserve(n_nodes_per_cluster[i]);
	}

	for (int i = 0; i < this->n_nodes_; ++i)
	{
		node_per_cluster[this->cluster_[i]].push_back(i);
	}

	return node_per_cluster;
}

void Clustering::set_cluster(int node, int cluster)
{
	this->cluster_[node] = cluster;
	this->n_clusters_ = std::max(this->n_clusters_, cluster + 1);
}

void Clustering::init_singleton_cluster()
{
	for (int i = 0; i < this->n_nodes_; ++i) {
		this->cluster_[i] = i;
	}
	this->n_clusters_ = this->n_nodes_;
}

void Clustering::order_clusters_by_n_nodes()
{
	class ClusterNNodes
	{
	public:
		int cluster_;
		int n_nodes_;

	public:

		ClusterNNodes() = default;

		ClusterNNodes(int cluster, int n_nodes) :
			cluster_(cluster),
			n_nodes_(n_nodes)
		{
		}

		bool operator<(const ClusterNNodes& cluster_n_nodes) const {
			return cluster_n_nodes.n_nodes_ < this->n_nodes_;
		}
	};

	std::vector<ClusterNNodes> cluster_n_nodes(this->n_clusters_);
	std::vector<int> n_nodes_per_cluster = get_n_nodes_per_cluster();

	for (int i = 0; i < this->n_clusters_; i++) {
		cluster_n_nodes[i] = ClusterNNodes(i, n_nodes_per_cluster[i]);
	}

	std::ranges::sort(cluster_n_nodes, std::less<ClusterNNodes>{});

	std::vector<int> new_cluster(this->n_clusters_, 0);

	int i = 0;

	do
	{
		new_cluster[cluster_n_nodes[i].cluster_] = i;
		++i;
	} while ((i < this->n_clusters_) && (cluster_n_nodes[i].n_nodes_ > 0));

	this->n_clusters_ = i;
	for (i = 0; i < this->n_nodes_; ++i) {
		this->cluster_[i] = new_cluster[this->cluster_[i]];
	}
}

void Clustering::merge_clusters(const Clustering& clustering)
{

	for (int i = 0; i < this->n_nodes_; ++i) {
		this->cluster_[i] = clustering.cluster_[this->cluster_[i]];
	}

	this->n_clusters_ = clustering.n_clusters_;
}



