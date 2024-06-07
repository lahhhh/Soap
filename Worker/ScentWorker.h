#pragma once

#include "Identifier.h"

#include "SparseInt.h"
#include "DenseInt.h"
#include "DenseDouble.h"

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

    bool read_ppi_database();

    bool do_integ_ppi();

    bool calculate_max_sr();

    bool calculate_sr();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_sr_ready(Eigen::ArrayXd);
};

