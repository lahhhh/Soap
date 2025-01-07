#pragma once

#include "Identifier.h"

#include <QProcess>

class MultiCommandWorker
	: public QObject
{
	Q_OBJECT
public:

	MultiCommandWorker(const QStringList& cmds) : commands_(cmds) {};

	QStringList commands_;

	QProcess* p_{ nullptr };

	int command_id_{ 0 };

public:

	bool work();

public slots:

	void run();

	void output();

	void error();

	void finished(int exit_code, QProcess::ExitStatus exit_status);

signals:

	void x_message(QString, int);

	void x_results_ready();
};

