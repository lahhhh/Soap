#pragma once

#include "Identifier.h"

class FileWritingWorker : public QObject
{
	Q_OBJECT
public:
	FileWritingWorker(
		void* data, 
		soap::VariableType data_type,
		const QString& file_name,
		soap::FileType file_type,
		const QStringList& settings = {}
	);

	void* data_ = nullptr;

	soap::VariableType data_type_;

	QString file_name_;

	soap::FileType file_type_;

	QStringList settings_;

public slots:

	void run();

signals:
	
	void x_message(QString, int);

	void x_results_ready();
	

};
