#pragma once

#include "Identifier.h"

#include "SingleCellAtac.h"

class LoadAtacFromFragmentsWorker
	: public QObject
{
	Q_OBJECT

public:
	explicit LoadAtacFromFragmentsWorker(const QString& file_name) : file_name_(file_name){}

	QString file_name_;

	std::unique_ptr<SingleCellAtac> single_cell_atac_;

	std::unordered_map<QString, int> barcode_index_;

	int n_barcode_{ 0 };

	GenomicRange peaks_;

	bool read_fragments();

	bool call_peak();

	bool calculate_matrix();

	void determine_species();

	void calculate_metadata();


public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_data_create_soon(void* data, soap::VariableType type, QString name);
};

