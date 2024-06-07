#pragma once

#include "Identifier.h"

#include <QString>

#include "SLM/VOSClusteringTechnique.h"

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
);
