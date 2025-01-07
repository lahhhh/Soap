#include "TfidfWorker.h"

#include "Custom.h"

bool TfidfWorker::work() {

	auto counts = this->counts_;

	this->res_.reset(new SparseDouble(
		SparseDouble::DataType::Normalized,
		counts->mat_.cast<double>(),
		counts->rownames_,
		counts->colnames_));

	Eigen::SparseMatrix<double>& normalized_data = this->res_->mat_;
	Eigen::ArrayXd column_sum = custom::col_sum_mt(normalized_data), idf = custom::row_sum(normalized_data);
	const int ncol = this->res_->cols();
	std::ranges::for_each(idf, [ncol](auto&& d) {if (d != 0) { d = ncol / d; }});

#pragma omp parallel for
	for (int i = 0; i < ncol; ++i) {
		double col_sum = column_sum[i];
		if (col_sum != 0) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(normalized_data, i); it; ++it) {
				it.valueRef() = log1p(it.value() / col_sum * idf[it.row()] * this->scale_factor_);
			}
		}
	};

	return true;
};

void TfidfWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_tfidf_ready(this->res_.release());

	G_TASK_END;
}

