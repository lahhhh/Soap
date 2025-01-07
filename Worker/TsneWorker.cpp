#include "TsneWorker.h"

#include "TSNE/Tsne.h"

TsneWorker::TsneWorker(
    const Eigen::MatrixXd& mat,
    unsigned int random_state,
    int max_iter,
    double perplexity,
    double theta,
    double momentum,
    double final_momentum,
    int stop_lying_iter,
    int momentum_switch_iter,
    double eta,
    double exaggeration_factor
):
    mat_(mat),
    random_state_(random_state),
    max_iter_(max_iter),
    perplexity_(perplexity),
    theta_(theta),
    momentum_(momentum),
    final_momentum_(final_momentum),
    stop_lying_iter_(stop_lying_iter),
    momentum_switch_iter_(momentum_switch_iter),
    eta_(eta),
    exaggeration_factor_(exaggeration_factor)
{};

bool TsneWorker::work() {

    this->mat_.transposeInPlace();

    Tsne tsne(
        this->random_state_,
        this->perplexity_,
        this->theta_,
        this->max_iter_,
        this->stop_lying_iter_,
        this->momentum_switch_iter_,
        this->momentum_,
        this->final_momentum_,
        this->eta_,
        this->exaggeration_factor_
    );

    int n_cell = this->mat_.cols();
    Eigen::MatrixXd embedding(2, n_cell);

    bool success = tsne.run(&this->mat_, embedding);

    if (!success) {
        G_TASK_WARN("tSNE failed.");
        return false;
    }

    this->res_ = embedding.transpose();

    return true;
};

void TsneWorker::run() {

    if (!this->work()) {
        G_TASK_END;
    }

    emit x_tsne_ready(this->res_);

    G_TASK_END;
}