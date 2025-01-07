#pragma once

#include "Identifier.h"

#include "Embedding.h"

#include <armadillo>

/*
* Port from Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh PR, 
* Raychaudhuri S. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods.
* and https://github.com/immunogenomics/harmony
*/

class HarmonyWorker :
    public QObject
{
    Q_OBJECT
public:
    HarmonyWorker(
        const Eigen::MatrixXd& mat, 
        const QList<QStringList> & metadata_list
    );

    void setup(
        const arma::mat& __Z, 
        const arma::sp_mat& __Phi,
        const arma::vec& __sigma, 
        const arma::vec& __theta,
        const arma::vec& __lambda, 
        const float __alpha, 
        const float __epsilon_kmeans, 
        const float __epsilon_harmony,
        const int __K
    );

    bool harmonize();

    /* METHODS */
    void moe_correct_ridge_cpp();
    arma::cube moe_ridge_get_betas_cpp();
    int cluster_cpp();

    void init_cluster_cpp();
    void allocate_buffers();
    void compute_objective();
    int update_R();
    bool check_convergence(int type);

    Eigen::MatrixXd mat_;
    QList<QStringList> metadata_list;

    /* FIELDS */
    arma::mat R, Z_orig, Z_corr, Z_cos, Y;
    arma::sp_mat Phi, Phi_moe, Phi_moe_t, Phi_t, Rk;
    arma::vec Pr_b, theta, N_b, sigma, lambda;

    // auxilary data structures
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross, objective_harmony;
    vector<int> kmeans_rounds, B_vec; // OLD: Kb
    std::vector<arma::uvec>index;

    float block_size, epsilon_kmeans, epsilon_harmony, alpha;
    unsigned int N, K, B, d, max_iter_kmeans, max_iter_harmony, window_size;

    // buffers
    arma::mat W, _scale_dist, dist_mat, O, E, dir_prior; // N_k, N_kb, N_b, numerator, denominator, C;
    arma::uvec update_order, cells_update;


    // flags
    bool ran_setup, ran_init, lambda_estimation, verbose; // do_merge_R;

public:

    bool work();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_harmony_ready(Eigen::MatrixXd embedding);
    
};

