#pragma once

#include "Identifier.h"

class PcaWorker2
    : public QObject
{
    Q_OBJECT
public:
    PcaWorker2(
        const Eigen::MatrixXd& mat,
        bool use_variable_features = true,
        int n_variable_feature = 2000
    ) :
        mat_(mat),
        use_variable_features_(use_variable_features),
        n_variable_feature_(n_variable_feature)
    {};

    Eigen::MatrixXd mat_;

    bool use_variable_features_{ true };

    int n_variable_feature_{ 0 };

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_pca_ready(Eigen::MatrixXd, QVector<double>, QVector<double>);
};
