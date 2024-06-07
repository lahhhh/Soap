#pragma once

#include "Identifier.h"

#include "DenseInt.h"
#include "DenseDouble.h"

class RnaNormalizeWorker 
	: public QObject
{
	Q_OBJECT
public:

	RnaNormalizeWorker(
		const DenseInt* counts,
		const QString& method
	):
		counts_(counts),
		method_(method)
	{}

	const DenseInt* counts_{ nullptr };

	QString method_;

	std::pair<QStringList, QVector<int>> get_transcript_length();

	bool fpkm();
	bool tpm();
	bool normalize_1();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_normalize_ready(DenseDouble);
};

