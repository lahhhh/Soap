#include "CommandWorker.h"

#include <QRegularExpression>

bool CommandWorker::work() {

    if (this->command_.isEmpty()) {
        G_TASK_WARN("Illegal Command.");
        return false;
    }

    this->p_ = new QProcess();

    connect(this->p_, &QProcess::readyReadStandardOutput, this, &CommandWorker::output);
    connect(this->p_, &QProcess::readyReadStandardError, this, &CommandWorker::error);
    connect(this->p_, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
        this, &CommandWorker::finished);

    this->p_->start("cmd.exe", QStringList{ "/C", this->command_ });

    return true;
};

void CommandWorker::run() {

    if (!this->work()) {
        G_TASK_END;
    }

    G_TASK_END;
};

void CommandWorker::finished(int exit_code, QProcess::ExitStatus exit_status) {

    G_TASK_LOG("Process [" + this->command_ + "] finished. Exit code : " + QString::number(exit_code) + ".");

    this->p_->close();
    this->p_->deleteLater();

    G_TASK_END;
}

void CommandWorker::output() {

    QByteArray data = this->p_->readAllStandardOutput();
    QString text = QString::fromLocal8Bit(data);
    text.remove(QRegularExpression("\x1B\\[[0-9;]*[A-Za-z]"));

    G_TASK_LOG(text);
};

void CommandWorker::error() {

    QByteArray data = this->p_->readAllStandardError();
    QString text = QString::fromLocal8Bit(data);
    text.remove(QRegularExpression("\x1B\\[[0-9;]*[A-Za-z]"));

    G_TASK_WARN(text);
};
