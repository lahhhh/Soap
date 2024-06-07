#pragma once

#include <vector>
#include <algorithm>

class Clustering
{

public:

    int n_nodes_;
    int n_clusters_;
    
    std::vector<int> cluster_;

public:

    Clustering() = default;
    Clustering(const Clustering&) = default;

    Clustering(int n_nodes);
    
    Clustering(const std::vector<int>& cluster);

    int get_n_nodes() const;

    int get_n_clusters() const;

    std::vector<int> get_clusters() const;

    int get_cluster(int node) const;

    std::vector<int> get_n_nodes_per_cluster() const;

    std::vector< std::vector<int> > get_nodes_per_cluster() const;

    void set_cluster(int node, int cluster);
    
    void init_singleton_cluster();
    
    void order_clusters_by_n_nodes();

    void merge_clusters(const Clustering& clustering);
};


