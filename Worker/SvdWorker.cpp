#include "SvdWorker.h"

#include "TruncatedSvd.h"
#include "Custom.h"

SvdWorker::SvdWorker(
	const Eigen::SparseMatrix<double>& mat,
	int var_perc,
	int n_mat_u,
	int random_state,
	double tol,
	int maximum_iteration
) :
	mat_(mat.transpose()),
	var_perc_(var_perc),
	n_mat_u_(n_mat_u),
	random_state_(random_state),
	tol_(tol),
	maximum_iteration_(maximum_iteration)
{}

bool SvdWorker::find_variable_features() {

	if (this->var_perc_ == 100) {
		return true;
	}

	auto var = custom::col_var_mt(this->mat_);

	double threshold = custom::linear_percentile(var, this->var_perc_);

	Eigen::ArrayX<bool> filter = var >= threshold;

	if (filter.count() < 100) {
		G_TASK_WARN("Variable Percentage is two small.");
		return false;
	}

	this->mat_ = custom::col_sliced(this->mat_, filter);

	return true;
};

void SvdWorker::run() {

	this->find_variable_features();

	Eigen::ArrayX<bool> filter = custom::col_sum(this->mat_) > 0;
	if (filter.count() != this->mat_.cols()) {
		this->mat_ = custom::col_sliced(this->mat_, filter);
	}

	auto [U, S, V] = tsvd(&this->mat_, this->n_mat_u_, this->random_state_, this->tol_, this->maximum_iteration_);

	Eigen::ArrayXd sdev = S.array() / sqrt(this->mat_.rows() - 1.);
	this->mat_.resize(0, 0);

	Eigen::MatrixXd emb = U * S.asDiagonal();

	emit x_svd_ready(emb, custom::cast<QVector>(sdev));
	G_TASK_END;
};
