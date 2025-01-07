#pragma once

#include "Identifier.h"

#include "SignalEmitter.h"
#include "SingleCellMultiome.h"
#include "SingleCellAtac.h"

class TfidfWorker 
    : public QObject
{
    Q_OBJECT
public:

    TfidfWorker(
        const SparseInt* counts,
        double scale_factor = 10000.0
    ) : counts_(counts), scale_factor_(scale_factor) {};

    const SparseInt* counts_{ nullptr };

    double scale_factor_{ 10000.0 };

    std::unique_ptr<SparseDouble> res_{ nullptr };

public:

    bool work();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_tfidf_ready(SparseDouble*);    
};

