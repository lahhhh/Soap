#pragma once

#include "Identifier.h"

/*
* modified from Scrublet package
* Wolock, Samuel L et al. ��Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data.�� 
* Cell systems vol. 8,4 (2019): 281-291.e9. doi:10.1016/j.cels.2018.11.005
*/

class ScrubletWorker 
    : public QObject
{
    Q_OBJECT
public:
    ScrubletWorker(
        const Eigen::SparseMatrix<int> & mat, 
        unsigned int random_state = 0, 
        double simulate_doublet_ratio = 2.0, 
        double expected_doublet_ratio = 0.06, 
        double stdev_doublet_ratio = 0.02, 
        double synthetic_doublet_umi_downsampling = 1.0
    );

    const Eigen::SparseMatrix<int> & mat_;

    int n_neighbors_;
    int n_cells_;
    int n_simulations_;
    int k_adjust_;

    unsigned int random_state_;

    double simulate_doublet_ratio_;
    double expected_doublet_ratio_;
    double stdev_doublet_ratio_;
    double synthetic_doublet_umi_downsampling_;

    Eigen::ArrayXd res1_, res2_;

public:

    bool work();

public slots:

    void run();

public:
    Eigen::SparseMatrix<double> simulate_doublet(const Eigen::SparseMatrix<double> & mat_);

    Eigen::SparseMatrix<double> normalize();

    void filter_gene(Eigen::SparseMatrix<double> &norm);

    Eigen::MatrixXd pca(Eigen::SparseMatrix<double> &norm, Eigen::SparseMatrix<double> &sim);

    void get_nearest_neighbors(Eigen::MatrixXi & idx, const Eigen::MatrixXd & pca__all);

    void calculate_score(Eigen::MatrixXi &idx);

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_scrublet_ready(Eigen::ArrayXd, Eigen::ArrayXd);
    
};
