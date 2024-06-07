#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"
#include "SparseInt.h"
#include "NumericMatrix.h"

class CellTypeAnnotationWorker : 
	public QObject
{
	Q_OBJECT

public:

	CellTypeAnnotationWorker(
		const SparseInt* counts,
		const QString& main_type_name,
		const QString& sub_type_name,
		bool annotate_by_cluster,
		QStringList cluster = {}
	) : 
		counts_(counts),
		main_type_name_(main_type_name),
		sub_type_name_(sub_type_name),
		annotate_by_cluster_(annotate_by_cluster),
		cluster_(cluster)
	{};

	const SparseInt* counts_;
	bool annotate_by_cluster_;
	QStringList cluster_;

	QString main_type_name_;
	QString sub_type_name_;

	Eigen::MatrixXd query_expression_;
	Eigen::MatrixXd celltype_expression_;

	QStringList main_type_database_;
	QStringList sub_type_database_;
	QStringList gene_name_database_;

	void load_database();

	bool filter_data();

	void assign_type();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_annotation_ready(
		QStringList main_type,
		QStringList sub_type,
		QString main_type_name,
		QString sub_type_name
	);
};

