#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class CountMatrixBatchLoadingWorker 
	: public QObject
{
	Q_OBJECT
public:
	CountMatrixBatchLoadingWorker(const QStringList& file_path, const QString& delimiter)
		:
		file_paths_(file_path),
		delimiter_(delimiter){}

	QStringList file_paths_;
	QString delimiter_;

	QList<SingleCellRna> objects_;

	std::unique_ptr<SingleCellRna> res_{ nullptr };

	bool load_single_object(SingleCellRna& object, const QString& file_path);

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

