#pragma once

#include "Identifier.h"

#include "SparseInt.h"
#include "DenseInt.h"
#include "DenseDouble.h"

/*
* modified from SCENT package
* Chen, Weiyan, and Andrew E Teschendorff. ¡°Estimating Differentiation Potency of Single Cells Using Single-Cell Entropy (SCENT).¡± 
* Methods in molecular biology (Clifton, N.J.) vol. 1935 (2019): 125-139. doi:10.1007/978-1-4939-9057-3_9
*/

class ScentWorker
    : public QObject
{
    Q_OBJECT
public:

    ScentWorker(const SparseInt* counts) : counts_(counts) {};

    const SparseInt* counts_{ nullptr };
    DenseDouble expression_;

    DenseInt ppi_;

    double max_sr_{ 0.0 };
    bool local_{ false };

    Eigen::ArrayXd res_;

    bool read_ppi_database();

    bool do_integ_ppi();

    bool calculate_max_sr();

    bool calculate_sr();

public:

    bool work();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_sr_ready(Eigen::ArrayXd);
};

