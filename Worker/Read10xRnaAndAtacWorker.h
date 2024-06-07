#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"

class Read10XMultiomeWorker :
    public QObject
{
    Q_OBJECT
public:
    Read10XMultiomeWorker(const QString& path_10X);

    QString path_10X_;

    SingleCellMultiome* single_cell_multiome_ = nullptr;

    Eigen::SparseMatrix<int> counts_;

    QStringList barcodes_, gene_symbols_, peak_names_, feature_names_;


public slots:

    void run();


signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_data_create_soon(void* data, soap::VariableType type, QString name);
       

private:

    bool read_barcodes();

    bool read_features();

    bool read_matrix();

    void determine_species();

    void separate_counts();

    void calculate_metadata();
};

