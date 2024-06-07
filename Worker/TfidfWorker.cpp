#include "TfidfWorker.h"

#include "Custom.h"

void TfidfWorker::run() {

	auto counts = this->counts_;

	SparseDouble* normalized = new SparseDouble(
		SparseDouble::DataType::Normalized,
		counts->mat_.cast<double>(),
		counts->rownames_,
		counts->colnames_);

	Eigen::SparseMatrix<double>& normalized_data = normalized->mat_;
	Eigen::ArrayXd column_sum = _Cs col_sum_mt(normalized_data), idf = _Cs row_sum(normalized_data);
	const int ncol = normalized->cols();
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

	emit x_tfidf_ready(normalized);

	G_TASK_END;
}

