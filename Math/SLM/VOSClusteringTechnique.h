#pragma once
/**
* modified from :
 * VOSClusteringTechnique
 *
 * @author Ludo Waltman
 * @author Nees Jan van Eck
 * @version 1.3.1, 11/23/14
 */

#include "Network.h"
#include "Clustering.h"

class VOSClusteringTechnique
{
public:
	Network network_;
	Clustering clustering_;
	double resolution_;

	VOSClusteringTechnique(const Network& network, double resolution);

	VOSClusteringTechnique(const Network& network, const Clustering& clustering, double resolution);

	Network get_network() const;

	Clustering get_clustering() const;

	double get_resolution() const;

	void set_network(const Network& network);

	void set_clustering(const Clustering& clustering);

	void set_resolution(double resolution);

	double calc_quality_function() const;


	bool run_local_moving_algorithm(std::default_random_engine& dre);

	bool run_louvain_algorithm(std::default_random_engine& dre);

	bool run_iterated_louvain_algorithm(int max_n_iterations, std::default_random_engine& dre);

	bool run_louvain_algorithm_with_multilevel_refinement(std::default_random_engine& dre);

	bool run_iterated_louvain_algorithm_with_multilevel_refinement(int max_n_iterations, std::default_random_engine& dre);


	bool run_smart_local_moving_algorithm(std::default_random_engine& dre);

	bool run_iterated_smart_local_moving_algorithm(int nIterations, std::default_random_engine& dre);

	int remove_cluster(int cluster);

	void remove_small_clusters(int min_n_nodes_per_cluster);
};
