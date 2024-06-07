#pragma once

#include "Identifier.h"

#include <QTextBrowser>

class InformationTextBrowser : 
	public QTextBrowser
{
	Q_OBJECT

public:
	InformationTextBrowser(QWidget* parant = nullptr);

	void log(const QString& info);
	void warn(const QString& info);
	void notice(const QString& info);

public slots:

	void s_receive_log(QString info);

	void s_receive_notice(QString info);

	void s_receive_warning(QString info);
};

