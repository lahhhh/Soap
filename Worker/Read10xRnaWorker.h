#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class Read10XRnaWorker : public QObject
{
    Q_OBJECT
public:

    explicit Read10XRnaWorker(const QString& path_10X) : path_10X_(path_10X) {}

    QString path_10X_;

    SingleCellRna* single_cell_rna_{ nullptr };

    QStringList barcodes_;
    QStringList gene_symbols_;

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

    void calculate_metadata();
};
