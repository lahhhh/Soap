#include "CommandWorker.h"

void CommandWorker::run() {

    if (this->command_.isEmpty()) {
        G_TASK_WARN("Illegal Command.");
        G_TASK_END;
    }

    QString program = this->command_[0];

    QStringList params;
    if (this->command_.size() > 1) {
        params = this->command_.sliced(1);
    }

    this->p_ = new QProcess();

    connect(this->p_, &QProcess::readyReadStandardOutput, this, &CommandWorker::output);
    connect(this->p_, &QProcess::readyReadStandardError, this, &CommandWorker::error);
    connect(this->p_, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
        this, &CommandWorker::finished);

    this->p_->start(program, params);
};

void CommandWorker::finished(int exit_code, QProcess::ExitStatus exit_status) {

    G_TASK_LOG("Process [" + this->command_.join(" ") + "] finished. Exit code : " + QString::number(exit_code) + ".");

    this->p_->close();
    this->p_->deleteLater();

    G_TASK_END;
}

void CommandWorker::output() {

    QByteArray data = this->p_->readAllStandardOutput();
    QString text = QString::fromLocal8Bit(data);

    G_TASK_LOG(text);
};

void CommandWorker::error() {

    QByteArray data = this->p_->readAllStandardError();
    QString text = QString::fromLocal8Bit(data);

    G_TASK_WARN(text);
};
