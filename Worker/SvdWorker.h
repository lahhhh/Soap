#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"

class SvdWorker : public QObject
{
    Q_OBJECT
public:
    SvdWorker(
        const Eigen::SparseMatrix<double>& mat,
        int var_perc = 100,
        int n_mat_u = 50,
        int random_state = 0,
        double tol = 0.00001,
        int maximum_iteration = 1000
    );

    Eigen::SparseMatrix<double> mat_;
    int var_perc_{ 100 };
    int n_mat_u_{ 50 };
    int random_state_{ 0 };
    double tol_{ 0.00001 };
    int maximum_iteration_{ 1000 };

    Eigen::MatrixXd res_;
    Eigen::ArrayXd sdev_;

    bool find_variable_features();

public:

    bool work();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_svd_ready(Eigen::MatrixXd, QVector<double>);
    
};
