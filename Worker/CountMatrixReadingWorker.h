#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class CountMatrixReadingWorker 
	: public QObject
{
	Q_OBJECT
public:
	CountMatrixReadingWorker(const QString& file_path, const QString& delimiter);

	QString file_path_;
	QString delimiter_;

	int nrow_{ 0 };
	int ncol_{ 0 };

	QStringList colnames_;

	QVector<Eigen::Triplet<int>> triplets_;

	QMap<int, QString> row_name_map_;

	std::unique_ptr<SingleCellRna> single_cell_rna_{ nullptr };

	bool read_file();

	bool create_data();

public slots:

	void run();

	bool load();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_data_create_soon(void* data, soap::VariableType type, QString name);	

};

