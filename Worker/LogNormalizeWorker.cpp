#include "LogNormalizeWorker.h"

#include "Custom.h"

void LogNormalizeWorker::run() {

	this->normalized_ = new SparseDouble();

	this->normalized_->rownames_ = this->data_->rownames_;
	this->normalized_->colnames_ = this->data_->colnames_;
	this->normalized_->mat_ = this->data_->mat_.cast<double>();
	this->normalized_->data_type_ = SparseDouble::DataType::Normalized;

	auto&& normalized_data = this->normalized_->mat_;

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

	emit x_log_normalize_ready(this->normalized_);

	G_TASK_END;
}