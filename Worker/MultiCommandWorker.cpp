#include "MultiCommandWorker.h"

bool MultiCommandWorker::work() {

    if (this->commands_.isEmpty()) {
        G_TASK_WARN("Empty Commands.");
        return false;
    }

    this->p_ = new QProcess();

    connect(this->p_, &QProcess::readyReadStandardOutput, this, &MultiCommandWorker::output);
    connect(this->p_, &QProcess::readyReadStandardError, this, &MultiCommandWorker::error);
    connect(this->p_, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
        this, &MultiCommandWorker::finished);

    QString cmd = this->commands_[this->command_id_];
    if (cmd.isEmpty()) {
        this->p_->close();

        this->p_->deleteLater();

        G_TASK_WARN("Empty Command.");
        return false;
    }

    this->p_->start("cmd.exe", { "/C", cmd });

    return true;
};

void MultiCommandWorker::run() {

    if (!this->work()) {
        G_TASK_END;
    }

    G_TASK_END;
};

void MultiCommandWorker::finished(int exit_code, QProcess::ExitStatus exit_status) {

    if (exit_status != QProcess::ExitStatus::NormalExit) {
        G_TASK_WARN("Meeting error in command " + QString::number(this->command_id_ + 1) + 
            ": [" + this->commands_[this->command_id_] + "].");
        this->p_->close();
        this->p_->deleteLater();
        G_TASK_END;
    }
    else {
        G_TASK_LOG("Command " + QString::number(this->command_id_ + 1) +
            ": [" + this->commands_[this->command_id_] + "] finished.");

        ++this->command_id_;

        if (this->command_id_ >= this->commands_.size()) {

            G_TASK_LOG("Script finished.");

            this->p_->close();
            this->p_->deleteLater();
            G_TASK_END;
        }
        else {
            QString cmd = this->commands_[this->command_id_];
            if (cmd.isEmpty()) {
                this->p_->close();

                this->p_->deleteLater();

                G_TASK_WARN("Empty Command " + QString::number(this->command_id_ + 1));
                G_TASK_END;
            }

            this->p_->start("cmd.exe", { "/C", cmd });
        }
    }
}

void MultiCommandWorker::output() {

    QByteArray data = this->p_->readAllStandardOutput();
    QString text = QString::fromLocal8Bit(data);

    G_TASK_LOG(text);
};

void MultiCommandWorker::error() {

    QByteArray data = this->p_->readAllStandardError();
    QString text = QString::fromLocal8Bit(data);

    G_TASK_WARN(text);
};
