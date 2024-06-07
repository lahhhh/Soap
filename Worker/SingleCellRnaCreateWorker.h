#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class SingleCellRnaCreateWorker
	: public QObject
{
	Q_OBJECT
public:
    SingleCellRnaCreateWorker(const SparseInt* si): mode_(WorkMode::FromSparseInt), si_(si){}

    enum class WorkMode { FromSparseInt };

    WorkMode mode_{ WorkMode::FromSparseInt };

    const SparseInt* si_{ nullptr };

    void create_from_sparseint();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_single_cell_rna_created(SingleCellRna* data);
};

