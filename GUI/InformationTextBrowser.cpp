#include "InformationTextBrowser.h"
#include <QDateTime>

InformationTextBrowser::InformationTextBrowser(QWidget* parent) : 
	QTextBrowser(parent)
{
	this->setFont(QFont("Arial", 13));
}

void InformationTextBrowser::log(const QString& info) {

	this->setTextColor(Qt::black);

	this->append(QDateTime::currentDateTime().toString("[yyyy-MM-dd ddd hh:mm:ss] ") + info);
}

void InformationTextBrowser::warn(const QString& info) {	

	this->setTextColor(Qt::red);

	this->append(QDateTime::currentDateTime().toString("[yyyy-MM-dd ddd hh:mm:ss] ") + info);
};

void InformationTextBrowser::notice(const QString& info) {

	this->setTextColor(Qt::darkCyan);

	this->append(QDateTime::currentDateTime().toString("[yyyy-MM-dd ddd hh:mm:ss] ") + info);
};

void InformationTextBrowser::s_receive_log(QString info) {

	this->log(info);
}

void InformationTextBrowser::s_receive_notice(QString info) {

	this->notice(info);
}

void InformationTextBrowser::s_receive_warning(QString info) {

	this->warn(info);
}