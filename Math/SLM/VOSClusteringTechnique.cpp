#include "VOSClusteringTechnique.h"

#include "Custom.h"

VOSClusteringTechnique::VOSClusteringTechnique(const Network& network, double resolution) :
	network_(network),
	clustering_(network.n_nodes_),
	resolution_(resolution)
{
	this->clustering_.init_singleton_cluster();
}

VOSClusteringTechnique::VOSClusteringTechnique(
	const Network& network,
	const Clustering& clustering,
	double resolution
) :
	network_(network),
	clustering_(clustering),
	resolution_(resolution)
{
}

Network VOSClusteringTechnique::get_network() const
{
	return this->network_;
}

Clustering VOSClusteringTechnique::get_clustering() const
{
	return this->clustering_;
}

double VOSClusteringTechnique::get_resolution() const
{
	return this->resolution_;
}

void VOSClusteringTechnique::set_network(const Network& network)
{
	this->network_ = network;
}

void VOSClusteringTechnique::set_clustering(const Clustering& clustering)
{
	this->clustering_ = clustering;
}

void VOSClusteringTechnique::set_resolution(double resolution)
{
	this->resolution_ = resolution;
}

double VOSClusteringTechnique::calc_quality_function() const
{
	double qualityFunction = 0;

	for (int i = 0; i < this->network_.n_nodes_; ++i)
	{
		int j = this->clustering_.cluster_[i];
		for (int k = this->network_.first_neighbor_index_[i]; k < this->network_.first_neighbor_index_[i + 1]; ++k)
			if (this->clustering_.cluster_[this->network_.neighbor_[k]] == j)
				qualityFunction += this->network_.edge_weight_[k];
	}
	qualityFunction += this->network_.total_edge_weight_self_links_;

	std::vector<double> cluster_weight(this->clustering_.n_clusters_, 0.0);
	for (int i = 0; i < this->network_.n_nodes_; ++i) {
		cluster_weight[this->clustering_.cluster_[i]] += this->network_.node_weight_[i];
	}
	for (int i = 0; i < this->clustering_.n_clusters_; ++i) {
		qualityFunction -= cluster_weight[i] * cluster_weight[i] * this->resolution_;
	}

	qualityFunction /= 2 * this->network_.get_total_edge_weight() + this->network_.total_edge_weight_self_links_;

	return qualityFunction;
}


bool VOSClusteringTechnique::run_local_moving_algorithm(std::default_random_engine& dre)
{
	bool update = false;
	double qualityFunction;
	int i;

	if (this->network_.n_nodes_ == 1)
		return false;

	update = false;

	std::vector<double> cluster_weight(this->network_.n_nodes_, 0.0);
	std::vector<int> n_nodes_per_cluster(this->network_.n_nodes_, 0);
	for (i = 0; i < this->network_.n_nodes_; ++i)
	{
		cluster_weight[this->clustering_.cluster_[i]] += this->network_.node_weight_[i];
		n_nodes_per_cluster[this->clustering_.cluster_[i]]++;
	}

	int n_unused_clusters = 0;
	std::vector<int> unused_cluster(this->network_.n_nodes_, 0);
	for (i = 0; i < this->network_.n_nodes_; ++i)
		if (n_nodes_per_cluster[i] == 0)
		{
			unused_cluster[n_unused_clusters] = i;
			n_unused_clusters++;
		}

	std::vector<int> node_permutation = custom::generate_random_permutation(this->network_.n_nodes_, dre);

	std::vector<double> edge_weight_per_cluster(this->network_.n_nodes_, 0);
	std::vector<int> neighboring_cluster(this->network_.n_nodes_ - 1, 0);
	int n_stable_nodes = 0;
	i = 0;
	do
	{
		int j = node_permutation[i];

		int n_neighboring_clusters = 0;
		for (int k = this->network_.first_neighbor_index_[j]; k < this->network_.first_neighbor_index_[j + 1]; ++k)
		{
			int l = this->clustering_.cluster_[this->network_.neighbor_[k]];
			if (edge_weight_per_cluster[l] == 0)
			{
				neighboring_cluster[n_neighboring_clusters] = l;
				n_neighboring_clusters++;
			}
			edge_weight_per_cluster[l] += this->network_.edge_weight_[k];
		}

		cluster_weight[this->clustering_.cluster_[j]] -= this->network_.node_weight_[j];
		n_nodes_per_cluster[this->clustering_.cluster_[j]]--;
		if (n_nodes_per_cluster[this->clustering_.cluster_[j]] == 0)
		{
			unused_cluster[n_unused_clusters] = this->clustering_.cluster_[j];
			n_unused_clusters++;
		}

		int best_cluster = -1;
		double maxQualityFunction = 0;
		for (int k = 0; k < n_neighboring_clusters; ++k)
		{
			int l = neighboring_cluster[k];
			qualityFunction = edge_weight_per_cluster[l] - this->network_.node_weight_[j] * cluster_weight[l] * this->resolution_;
			if ((qualityFunction > maxQualityFunction) || ((qualityFunction == maxQualityFunction) && (l < best_cluster)))
			{
				best_cluster = l;
				maxQualityFunction = qualityFunction;
			}
			edge_weight_per_cluster[l] = 0;
		}
		if (maxQualityFunction == 0)
		{
			best_cluster = unused_cluster[n_unused_clusters - 1];
			n_unused_clusters--;
		}

		cluster_weight[best_cluster] += this->network_.node_weight_[j];
		n_nodes_per_cluster[best_cluster]++;
		if (best_cluster == this->clustering_.cluster_[j])
			n_stable_nodes++;
		else
		{
			this->clustering_.cluster_[j] = best_cluster;
			n_stable_nodes = 1;
			update = true;
		}

		i = (i < this->network_.n_nodes_ - 1) ? (i + 1) : 0;
	} while (n_stable_nodes < this->network_.n_nodes_);

	std::vector<int> new_cluster(this->network_.n_nodes_, 0);
	this->clustering_.n_clusters_ = 0;
	for (i = 0; i < this->network_.n_nodes_; ++i)
		if (n_nodes_per_cluster[i] > 0)
		{
			new_cluster[i] = this->clustering_.n_clusters_;
			this->clustering_.n_clusters_++;
		}
	for (i = 0; i < this->network_.n_nodes_; ++i)
		this->clustering_.cluster_[i] = new_cluster[this->clustering_.cluster_[i]];

	return update;
}


bool VOSClusteringTechnique::run_louvain_algorithm(std::default_random_engine& dre)
{

	if (this->network_.n_nodes_ == 1)
		return false;

	bool update = run_local_moving_algorithm(dre);

	if (this->clustering_.n_clusters_ < this->network_.n_nodes_)
	{
		VOSClusteringTechnique VOSClusteringTechnique(this->network_.create_reduced_network(this->clustering_), this->resolution_);

		bool update2 = VOSClusteringTechnique.run_louvain_algorithm(dre);

		if (update2)
		{
			update = true;

			this->clustering_.merge_clusters(VOSClusteringTechnique.clustering_);
		}
	}

	return update;
}


bool VOSClusteringTechnique::run_iterated_louvain_algorithm(int max_n_iterations, std::default_random_engine& dre)
{
	bool update;

	int i = 0;
	do
	{
		update = run_louvain_algorithm(dre);
		++i;
	} while ((i < max_n_iterations) && update);
	return ((i > 1) || update);
}

bool VOSClusteringTechnique::run_louvain_algorithm_with_multilevel_refinement(std::default_random_engine& dre)
{

	if (this->network_.n_nodes_ == 1)
		return false;

	bool update = run_local_moving_algorithm(dre);

	if (this->clustering_.n_clusters_ < this->network_.n_nodes_)
	{
		VOSClusteringTechnique VOSClusteringTechnique(this->network_.create_reduced_network(this->clustering_), this->resolution_);

		bool update2 = VOSClusteringTechnique.run_louvain_algorithm_with_multilevel_refinement(dre);

		if (update2)
		{
			update = true;

			this->clustering_.merge_clusters(VOSClusteringTechnique.clustering_);

			run_local_moving_algorithm(dre);
		}
	}

	return update;
}

bool VOSClusteringTechnique::run_iterated_louvain_algorithm_with_multilevel_refinement(int max_n_iterations, std::default_random_engine& dre)
{
	bool update;
	int i;

	i = 0;
	do
	{
		update = run_louvain_algorithm_with_multilevel_refinement(dre);
		++i;
	} while ((i < max_n_iterations) && update);
	return ((i > 1) || update);
}

bool VOSClusteringTechnique::run_smart_local_moving_algorithm(std::default_random_engine& dre)
{
	int i = 0;

	if (this->network_.n_nodes_ == 1)
		return false;

	bool update = run_local_moving_algorithm(dre);

	if (this->clustering_.n_clusters_ < this->network_.n_nodes_)
	{
		std::vector<Network> subnetwork = this->network_.create_subnetworks(this->clustering_);

		std::vector< std::vector<int> > node_per_cluster = this->clustering_.get_nodes_per_cluster();

		this->clustering_.n_clusters_ = 0;
		std::vector<int>  n_nodes_per_clusterReduced_network(subnetwork.size(), 0);
		for (i = 0; i < subnetwork.size(); ++i)
		{
			VOSClusteringTechnique VOSClusteringTechnique(subnetwork[i], this->resolution_);

			VOSClusteringTechnique.run_local_moving_algorithm(dre);

			for (int j = 0; j < subnetwork[i].n_nodes_; ++j)
				this->clustering_.cluster_[node_per_cluster[i][j]] = this->clustering_.n_clusters_ + VOSClusteringTechnique.clustering_.cluster_[j];
			this->clustering_.n_clusters_ += VOSClusteringTechnique.clustering_.n_clusters_;
			n_nodes_per_clusterReduced_network[i] = VOSClusteringTechnique.clustering_.n_clusters_;
		}

		VOSClusteringTechnique VOSClusteringTechnique(this->network_.create_reduced_network(this->clustering_), this->resolution_);

		i = 0;
		for (int j = 0; j < n_nodes_per_clusterReduced_network.size(); ++j)
			for (int k = 0; k < n_nodes_per_clusterReduced_network[j]; ++k)
			{
				VOSClusteringTechnique.clustering_.cluster_[i] = j;
				++i;
			}
		VOSClusteringTechnique.clustering_.n_clusters_ = n_nodes_per_clusterReduced_network.size();

		update |= VOSClusteringTechnique.run_smart_local_moving_algorithm(dre);

		this->clustering_.merge_clusters(VOSClusteringTechnique.clustering_);
	}

	return update;
}

bool VOSClusteringTechnique::run_iterated_smart_local_moving_algorithm(int nIterations, std::default_random_engine& dre)
{
	bool update;
	int i;

	update = false;
	for (i = 0; i < nIterations; ++i)
		update |= run_smart_local_moving_algorithm(dre);
	return update;
}

int VOSClusteringTechnique::remove_cluster(int cluster)
{
	int i, j;

	std::vector<double> cluster_weight(this->clustering_.n_clusters_, 0.0);
	std::vector<double> totalEdge_weight_per_cluster(this->clustering_.n_clusters_, 0.0);
	for (i = 0; i < this->network_.n_nodes_; ++i)
	{
		cluster_weight[this->clustering_.cluster_[i]] += this->network_.node_weight_[i];
		if (this->clustering_.cluster_[i] == cluster)
			for (j = this->network_.first_neighbor_index_[i]; j < this->network_.first_neighbor_index_[i + 1]; ++j)
				totalEdge_weight_per_cluster[this->clustering_.cluster_[this->network_.neighbor_[j]]] += this->network_.edge_weight_[j];
	}

	i = -1;
	double maxQualityFunction = 0;
	for (j = 0; j < this->clustering_.n_clusters_; ++j)
		if ((j != cluster) && (cluster_weight[j] > 0))
		{
			double qualityFunction = totalEdge_weight_per_cluster[j] / cluster_weight[j];
			if (qualityFunction > maxQualityFunction)
			{
				i = j;
				maxQualityFunction = qualityFunction;
			}
		}

	if (i >= 0)
	{
		for (j = 0; j < this->network_.n_nodes_; ++j)
			if (this->clustering_.cluster_[j] == cluster)
				this->clustering_.cluster_[j] = i;
		if (cluster == this->clustering_.n_clusters_ - 1)
			this->clustering_.n_clusters_ = std::ranges::max(this->clustering_.cluster_) + 1;
	}

	return i;
}

void VOSClusteringTechnique::remove_small_clusters(int min_n_nodes_per_cluster)
{
	int i, j;

	VOSClusteringTechnique VOSClusteringTechnique(this->network_.create_reduced_network(this->clustering_), this->resolution_);

	std::vector<int> n_nodes_per_cluster = this->clustering_.get_n_nodes_per_cluster();

	do
	{
		i = -1;
		j = min_n_nodes_per_cluster;
		for (int k = 0; k < VOSClusteringTechnique.clustering_.n_clusters_; ++k)
			if ((n_nodes_per_cluster[k] > 0) && (n_nodes_per_cluster[k] < j))
			{
				i = k;
				j = n_nodes_per_cluster[k];
			}

		if (i >= 0)
		{
			j = VOSClusteringTechnique.remove_cluster(i);
			if (j >= 0)
				n_nodes_per_cluster[j] += n_nodes_per_cluster[i];
			n_nodes_per_cluster[i] = 0;
		}
	} while (i >= 0);

	this->clustering_.merge_clusters(VOSClusteringTechnique.clustering_);
}