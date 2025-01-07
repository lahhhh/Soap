#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"

class TableReadingWorker 
	: public QObject
{
	Q_OBJECT
public:

	TableReadingWorker(const QString& file_name, bool fast) :
		file_name_(file_name),
		fast_(fast)
	{}

	QString file_name_;

	bool fast_{ false };

	std::unique_ptr<CustomMatrix> res_{ nullptr };

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_data_frame_ready(CustomMatrix*);
	
};
