#pragma once

/*
modified from R package 'uwot'
*/

#include "Identifier.h"

#include <QString>
#include <mutex>

class UMAP
{
public:
    UMAP(
        Eigen::MatrixXd* mat, 
        int n_neighbors = 30, 
        const QString& metric = "Angular", 
        double learning_rate = 1.0, 
        const QString& init = "spectral",
        double minimum_distance = 0.3, 
        double spread = 1.0, 
        double set_op_mix_ratio = 1.0,
        double repulsion_strength = 1.0, 
        int negative_sample_rate = 5, 
        int random_state = 1997, 
        int n_trees = 50
    );

    void fit();

    void find_a_b_param(double spread_, double minimum_distance_);

    void set_disconnection_distance();

    void smooth_knn_distance(int n_iter = 64);

    void compute_membership_strengths();

    void fuzzy_simplicial_set();

    void fit_embed_data();

    void spectral_layout();

    void make_epoches_per_sample();

    void optimize_layout_euclidean();

    Eigen::MatrixXd* mat_{nullptr};

    QString metric_;
    QString init_;    
    
    double learning_rate_ = 1.0;
    double minimum_distance_ = 1.0;
    double spread_ = 1.0;
    double set_op_mix_ratio_ = 1.0;
    double local_connectivity_ = 1.0;
    double repulsion_strength_ = 1.0;
    double a_ = 0.0;
    double b_ = 0.0;
    double initial_alpha_ = 0.0;
    double disconnection_distance_ = 0.0;

    int n_neighbors_ = 15;
    int n_components_ = 2;
    int n_epoches_ = 0;
    int negative_sample_rate_ = 5;
    int random_state_ = 1997;
    int n_trees_ = 50;
    int n_vertices_ = 0;

    Eigen::MatrixXd knn_distance_;
    Eigen::MatrixXi knn_indices_;

    Eigen::ArrayXd smooth_sigmas_;
    Eigen::ArrayXd smooth_rho_;

    Eigen::SparseMatrix<double> strengths_;
    Eigen::ArrayXXd embedding_;
    Eigen::ArrayXd epochs_per_sample_;

    std::mutex mutex_;
};

