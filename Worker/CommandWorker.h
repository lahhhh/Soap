#pragma once

#include "Identifier.h"

#include <QProcess>

class CommandWorker
	: public QObject
{
	Q_OBJECT
public:

	CommandWorker(const QString& cmd) : command_(cmd) {};

	QString command_;

	QProcess* p_{ nullptr };

public slots:

	void run();

public:

	bool work();

	void output();

	void error();

	void finished(int exit_code, QProcess::ExitStatus exit_status);

signals:

	void x_message(QString, int);

	void x_results_ready();
};

