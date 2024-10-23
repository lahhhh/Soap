#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class Read10XRnaWorker : public QObject
{
    Q_OBJECT
public:

	Read10XRnaWorker(
		const QString& barcodes_file_name,
		const QString& features_file_name,
		const QString& matrix_file_name) :
		barcodes_file_name_(barcodes_file_name),
		features_file_name_(features_file_name),
		matrix_file_name_(matrix_file_name)
	{}

	QString barcodes_file_name_;
	QString features_file_name_;
	QString matrix_file_name_;

    std::unique_ptr<SingleCellRna> single_cell_rna_;

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
