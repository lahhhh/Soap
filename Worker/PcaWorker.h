#pragma once

#include "Identifier.h"

Eigen::MatrixXd pca_infercnv(
    const Eigen::MatrixXd& mat,
    int n_variable_feature = 2000,
    int nu = 50,
    int random_state = 0,
    double tol = 0.00001,
    int max_iter = 1000);

class PcaWorker
    : public QObject
{
    Q_OBJECT
public:
    PcaWorker(
        const Eigen::SparseMatrix<int>* mat,
        double feature_proportion, 
        int n_mat_u = 50, 
        int random_state = 0, 
        double tol = 0.00001, 
        int maximum_iteration = 1000
    ) :
        mat_(mat),
        feature_proportion_(feature_proportion),
        n_mat_u_(n_mat_u),
        random_state_(random_state),
        tol_(tol),
        maximum_iteration_(maximum_iteration)
    {};

    const Eigen::SparseMatrix<int>* mat_ = nullptr;
    double feature_proportion_ = 0.;

    int n_mat_u_;
    int random_state_;
    double tol_;
    int maximum_iteration_;

    Eigen::MatrixXd res_;
    QVector<double> sdev_;

public:

    bool work();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_pca_ready(Eigen::MatrixXd, QVector<double>);    
};
