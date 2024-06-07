#pragma once

#include "Identifier.h"

#include "CellChat.h"
#include "Custom.h"
#include "SparseDouble.h"

class CellChatWorker :
	public QObject
{
	Q_OBJECT
public:
	CellChatWorker(
		const QString& identity,
		const SparseDouble& normalized,
		soap::Species species,
		const QStringList& metadata,
		double minimum_percentage,
		const QString& p_adjust_method,
		int random_state,
		int n_boot,
		int minimum_cell_number,
		const QString& annotation_type);

	SparseDouble normalized_;

	QString identity_;
	QString p_adjust_method_;
	QString annotation_type_;

	soap::Species species_;

	QStringList metadata_;

	double minimum_percentage_;

	int random_state_;
	int n_boot_;
	int minimum_cell_number_;

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_cellchat_ready(CellChat);

};
