#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class BatchLoading10XRnaWorker
	: public QObject
{
	Q_OBJECT
public:

	BatchLoading10XRnaWorker(const QString& dir) : dir_(dir) {};

	QString dir_;

	QList<SingleCellRna> objects_;

	std::unique_ptr<SingleCellRna> res_{ nullptr };

	bool load_single_object(
		SingleCellRna& object, 
		const QString& path,
		const QString& barcodes_file_name,
		const QString& features_file_name,
		const QString& matrix_file_name);

	void integrate_objects();

	QList<std::pair<QVector<int>, QVector<int> > >
		integrate_sparseint(
			SparseInt& to,
			QList<SparseInt const*> froms,
			bool distinguish_barcode
		);

	void integrate_metadata(Metadata& to, QList<Metadata*> froms);

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_data_create_soon(void* data, soap::VariableType type, QString name);

};
