#pragma once

#include "Identifier.h"

#include "SignalEmitter.h"
#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

class LogNormalizeWorker 
	: public QObject
{
	Q_OBJECT
public:
	LogNormalizeWorker(
		const SparseInt* data,
		double scale_factor = 10000.0
	) :
		data_(data),
		scale_factor_(scale_factor) 
	{};

	const SparseInt* data_{ nullptr };
	SparseDouble* normalized_{ nullptr };

	double scale_factor_{ 10000.0 };

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_log_normalize_ready(SparseDouble*);

};

