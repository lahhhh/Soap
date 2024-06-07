#pragma once

#include <vector>
#include <QVector>
#include <random>

#include "Clustering.h"

class Network
{
public:
	int n_nodes_ = 0;
	int n_edges_ = 0;

	double total_edge_weight_self_links_ = 0;

	std::vector<double> node_weight_;
	std::vector<double> edge_weight_;
	std::vector<int> first_neighbor_index_;
	std::vector<int> neighbor_;

public:

	Network() = default;

	Network(int n_nodes, const std::vector< std::vector<int> >& edge);

	Network(
		int n_nodes, 
		std::vector<double>* node_weight, 
		const std::vector< std::vector<int> >& edge
	);

	Network(
		int n_nodes, 
		const std::vector< std::vector<int> >& edge,
		std::vector<double>* edge_weight
	);

	Network(
		int n_nodes, 
		std::vector<double>* node_weight, 
		const std::vector< std::vector<int> >& edge, 
		std::vector<double>* edge_weight
	);

	Network(
		int n_nodes, 
		std::vector<int>* first_neighbor_index, 
		std::vector<int>* neighbor
	);

	Network(
		int n_nodes, 
		std::vector<double> *node_weight, 
		std::vector<int>* first_neighbor_index, 
		std::vector<int>* neighbor
	);

	Network(
		int n_nodes, 
		std::vector<int>* first_neighbor_index, 
		std::vector<int>* neighbor, 
		std::vector<double>* edge_weight
	);

	void set_from(int n_nodes,
		const std::vector<int>& first_neighbor_index,
		const std::vector<int>& neighbor,
		const std::vector<double>& edge_weight
	);

	Network(
		int n_nodes, 
		std::vector<double>* node_weight, 
		std::vector<int>* first_neighbor_index, 
		std::vector<int>* neighbor, 
		std::vector<double>* edge_weight
	);

	void set_from(int n_nodes,
		const std::vector<double>& node_weight,
		const std::vector<int>& first_neighbor_index,
		const std::vector<int>& neighbor,
		const std::vector<double>& edge_weight
	);

	Network(
		int n_nodes,
		int n_edges,
		double total_edge_weight_self_links,
		const std::vector<double>& node_weight,
		const std::vector<double>& edge_weight,
		const std::vector<int>& first_neighbor_index,
		const std::vector<int>& neighbor
	);

	int get_n_nodes() const;

	double get_total_node_weight() const;

	std::vector<double> get_node_weights() const;

	double get_node_weight(int node) const;

	int get_n_edges() const;

	int get_n_edges(int node) const;

	std::vector<int> get_n_edges_per_node() const;

	std::vector< std::vector<int> > get_edges() const;

	std::vector<int> get_edges(int node) const;

	std::vector< std::vector<int> > get_edges_per_node() const;

	double get_total_edge_weight() const;

	double get_total_edge_weight(int node) const;

	std::vector<double> get_total_edge_weight_per_node() const;

	std::vector<double> get_edge_weights() const;

	std::vector<double> get_edge_weights(int node) const;

	std::vector< std::vector<double> > get_edge_weights_per_node() const;

	double get_total_edge_weight_self_links() const;

	Network create_network_without_node_weights() const;

	Network create_network_without_edge_weights() const;

	Network create_network_without_node_and_edge_weights() const;

	Network create_normalized_network1() const;

	Network create_normalized_network2() const;

	Network create_pruned_network(int n_edges) const;

	Network create_pruned_network(int n_edges, std::default_random_engine& re) const;

	Network create_subnetwork(const std::vector<int>& node) const;

	Network create_subnetwork(const QVector<bool>& node_in_subnetwork) const;

	Network create_subnetwork(const Clustering& clustering, int cluster) const;

	std::vector<Network> create_subnetworks(const Clustering& clustering) const;

	Network create_subnetwork_largest_component() const;

	Network create_reduced_network(const Clustering& clustering) const;

	Clustering identify_components() const;

	double generate_random_number(int node1, int node2, const std::vector<int>& node_permutation) const;

	Network create_subnetwork(
		const Clustering& clustering, 
		int cluster, 
		const std::vector<int>& node, 
		std::vector<int> subnetwork_node, 
		std::vector<int> subnetwork_neighbor, 
		std::vector<double> subnetworkEdge_weight) const;
};

