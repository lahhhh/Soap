#pragma once

#include "Identifier.h"
#include <vector>

#include "datapoint.h"

/*
 modified from Rtsne: https://github.com/jkrijthe/Rtsne
*/

static inline double sign_tsne(double x) { return (x == 0.0 ? 0.0 : (x < 0.0 ? -1.0 : 1.0)); }

class Tsne
{
public:

    Tsne(
        unsigned int random_state,
        double perplexity, 
        double theta, 
        int max_iter, 
        int stop_lying_iter,
        int mom_switch_iter, 
        double momentum, 
        double final_momentum, 
        double eta, 
        double exaggeration_factor
    );

    bool run(Eigen::MatrixXd* mat, Eigen::MatrixXd& embedding);

private:

    void symmetrize_matrix(unsigned int n_vertice);
    void train_iterations(int n_vertice, Eigen::MatrixXd& output);

    void compute_gradient(const Eigen::MatrixXd& output, Eigen::ArrayXXd& d_output);
    void compute_exact_gradient(Eigen::MatrixXd& output, int n_vertice, Eigen::ArrayXXd& d_output);
    void zero_mean(Eigen::MatrixXd& output);

    void compute_gaussian_perplexity(const Eigen::MatrixXd& data);
    void compute_gaussian_perplexity(const Eigen::MatrixXd& data, int K);
    void setup_approximate_memory(unsigned int n_vertice, int K);

    void compute_probabilities(const double perplexity, const int K, const double* distances, double* cur_P);
    void compute_squared_euclidean_distance(const Eigen::MatrixXd& output, Eigen::ArrayXXd& distance);

    Eigen::MatrixXd* mat_{ nullptr };

    unsigned int random_state_;

    // Member variables.
    double perplexity_;
    double theta_;
    double momentum_;
    double final_momentum_;
    double eta_;
    double exaggeration_factor_;

    int max_iter_;
    int stop_lying_iter_;
    int mom_switch_iter_;

    bool exact_;

    std::vector<unsigned int> row_P, col_P;
    std::vector<double> val_P;
    Eigen::MatrixXd P;
};

