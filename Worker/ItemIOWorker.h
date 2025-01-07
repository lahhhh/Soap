#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"
#include "BulkRna.h"
#include "DataFrame.h"

class ItemIOWorker 
	: public QObject
{
	Q_OBJECT
public:
	ItemIOWorker(
		const QString& file_path, 
		soap::VariableType data_type, 
		void* data, 
		const QString& item_name
	);

	ItemIOWorker(const QString& file_path);

	enum class WorkMode{Read, Write};

	WorkMode mode_ = WorkMode::Read;

	QString file_path_;
	QString item_name_;

	soap::VariableType data_type_;

	void* data_;

	void read();

	void write();

public:

	bool work();

public slots:

	void run();

signals:
	
	void x_message(QString, int);
	
	void x_results_ready();

	void x_data_create_soon(void* data, soap::VariableType type, QString name);
	
};

