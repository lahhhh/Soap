#include "Louvain.h"

#include "Custom.h"

Eigen::SparseMatrix<double> create_shared_nearest_neighbors_matrix(
    const Eigen::MatrixXd& mat,
    const QString& nn_method,
    const int n_neighbors,
    const int n_trees
) {

    Eigen::MatrixXi knn;

    if (nn_method == "Euclidean") {
        knn = _Cs get_knn_mt<Euclidean>(mat, n_neighbors, n_trees);
    }
    else if (nn_method == "Angular") {
        knn = _Cs get_knn_mt<Angular>(mat, n_neighbors, n_trees);
    }
    else /* if (nn_method == "Manhattan") */ {
        knn = _Cs get_knn_mt<Manhattan>(mat, n_neighbors, n_trees);
    }

    return _Cs create_shared_nearest_neighbors_matrix(knn);
};

Network create_network(const Eigen::SparseMatrix<double>& snn,
    const int modularity_function) {

    const int n_nodes = snn.rows();

    std::vector<int> node1, node2;
    std::vector<double> edge_weight;

    for (Eigen::Index i = 0; i < n_nodes; ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
            if (it.row() < i) {
                node1.emplace_back(it.row());
                node2.emplace_back(i);
                edge_weight.emplace_back(it.value());
            }
        }
    }

    std::vector<int> n_neighbors(n_nodes), first_neighbor_index(n_nodes + 1);

    for (std::size_t i = 0; i < node1.size(); ++i) {
        ++n_neighbors[node1[i]];
        ++n_neighbors[node2[i]];
    }

    int n_edges = 0;
    for (int i = 0; i < n_nodes; ++i) {
        first_neighbor_index[i] = n_edges;
        n_edges += n_neighbors[i]; // n_edges = number of edges * 2
    }

    first_neighbor_index[n_nodes] = n_edges;

    std::vector<int> neighbor(n_edges);
    std::vector<double> edge_weight2(n_edges);
    std::fill(n_neighbors.begin(), n_neighbors.end(), 0);

    // resort node and edges
    for (std::size_t i = 0; i < node1.size(); ++i) {

        int node2_location = first_neighbor_index[node1[i]] + n_neighbors[node1[i]];
        neighbor[node2_location] = node2[i];
        edge_weight2[node2_location] = edge_weight[i];
        ++n_neighbors[node1[i]];

        int node1_location = first_neighbor_index[node2[i]] + n_neighbors[node2[i]];
        neighbor[node1_location] = node1[i];
        edge_weight2[node1_location] = edge_weight[i];
        ++n_neighbors[node2[i]];
    }

    Network network;

    if (modularity_function == 1) {
        network.set_from(n_nodes, first_neighbor_index, neighbor, edge_weight2);
    }
    else {
        network.set_from(n_nodes, std::vector<double>(n_nodes, 1.0), first_neighbor_index, neighbor, edge_weight2);
    }

    return network;
}

std::vector<int> louvain_cluster(
    const Eigen::MatrixXd& mat,
    const QString& method,
    const QString& nn_method,
    int modularity_function,
    int n_random_start,
    int n_iterations,
    int random_seed,
    int n_neighbors,
    int n_trees,
    double resolution
) {

    auto snn = create_shared_nearest_neighbors_matrix(mat, nn_method, n_neighbors, n_trees);

    auto network = create_network(snn, modularity_function);

    if (modularity_function == 1) {
        resolution /= (2 * network.get_total_edge_weight() + network.get_total_edge_weight_self_links());
    }

    double max_modularity = -std::numeric_limits<double>::infinity(), modularity = 0;

    std::default_random_engine dre;
    dre.seed(random_seed);

    Clustering clustering;

    for (int i = 0; i < n_random_start; ++i) {
        VOSClusteringTechnique vct(network, resolution);

        int j = 0;
        bool update = false;

        do {
            if (method == "Louvain") {
                update = vct.run_louvain_algorithm(dre);
            }
            else if (method == "Modified Louvain") {
                update = vct.run_louvain_algorithm_with_multilevel_refinement(dre);
            }
            else if (method == "SLM") {
                update = vct.run_smart_local_moving_algorithm(dre);
            }

            ++j;
            modularity = vct.calc_quality_function();

        } while (j < n_iterations && update);

        if (modularity > max_modularity) {
            clustering = vct.get_clustering();
            max_modularity = modularity;
        }
    }

    clustering.order_clusters_by_n_nodes();

    return clustering.cluster_;
};