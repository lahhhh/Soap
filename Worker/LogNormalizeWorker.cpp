#include "LogNormalizeWorker.h"

#include "Custom.h"

bool LogNormalizeWorker::work() {

	this->res_.reset(new SparseDouble());

	this->res_->rownames_ = this->data_->rownames_;
	this->res_->colnames_ = this->data_->colnames_;
	this->res_->mat_ = this->data_->mat_.cast<double>();
	this->res_->data_type_ = SparseDouble::DataType::Normalized;

	auto&& normalized_data = this->res_->mat_;

	Eigen::ArrayXd column_sum = custom::col_sum_mt(normalized_data);
	int ncol = normalized_data.cols();

#pragma omp parallel for
	for (int i = 0; i < ncol; ++i) {
		double col_sum = column_sum[i];
		if (col_sum != 0) {
			double factor = this->scale_factor_ / col_sum;
			for (Eigen::SparseMatrix<double>::InnerIterator it(normalized_data, i); it; ++it) {
				it.valueRef() = log1p(double(it.value()) * factor);
			}
		}
	}

	return true;
};

void LogNormalizeWorker::run() {

	G_TASK_LOG("Start normalization...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_log_normalize_ready(this->res_.release());

	G_TASK_LOG("Normalization finished.");

	G_TASK_END;
}