#include "CalculateCountsByGenomicRangeWorker.h"

#include "Custom.h"
#include "GenomeUtility.h"
#include <zlib.h>

CalculateCountsByGenomicRangeWorker::CalculateCountsByGenomicRangeWorker(
	const Fragments* fragments,
	const GenomicRange& genomic_range
) :
	fragments_(fragments),
	genomic_range_(genomic_range)
{};

SparseInt* CalculateCountsByGenomicRangeWorker::calculate_counts(const Fragments* fragments, const GenomicRange& genomic_range) {
	CalculateCountsByGenomicRangeWorker worker(fragments, genomic_range);
	worker.create_index();
	if (!worker.calculate_counts()) {
		return nullptr;
	}
	worker.build_matrix();
	return worker.counts_;
};

void CalculateCountsByGenomicRangeWorker::run() {

	G_TASK_LOG("Creating index...");

	this->create_index();
	G_TASK_LOG("Mapping fragments...");

	if (!this->calculate_counts()) {
		G_TASK_END;
	}
	G_TASK_LOG("Building counts matrix...");
	this->build_matrix();
	emit x_peak_counts_ready(this->counts_);

	G_TASK_LOG("Counts calculation finished.");

	G_TASK_END;
}

void CalculateCountsByGenomicRangeWorker::create_index() {

	this->genome_.set(this->genomic_range_);

	int n_range = this->genomic_range_.rows();
	const int n_cell = this->fragments_->cell_names_.size();

	this->dense_counts_.resize(n_range, n_cell);

};

void CalculateCountsByGenomicRangeWorker::build_matrix() {


	this->counts_ = new SparseInt();
	this->counts_->data_type_ = SparseInt::DataType::Counts;

	this->counts_->rownames_ = this->genomic_range_.get_range_names();


	this->counts_->colnames_ = this->fragments_->cell_names_;


	Eigen::SparseMatrix<int>& counts_matrix = this->counts_->mat_;

	counts_matrix = this->dense_counts_.sparseView();

	this->dense_counts_.resize(0, 0);

	auto filter = custom::greater_than(custom::row_sum(counts_matrix), 0);

	this->counts_->row_slice(filter);
}

void CalculateCountsByGenomicRangeWorker::find_row(const QString& seq_name, int cell_loc, int start, int end) {

	auto [r1, r2] = this->genome_.find_location(seq_name, start, end);

	if (this->genome_.success(r1)) {
		++this->dense_counts_(r1, cell_loc);
	}

	if (this->genome_.success(r2)) {
		++this->dense_counts_(r2, cell_loc);
	}
};

bool CalculateCountsByGenomicRangeWorker::calculate_counts() {

	for (const auto& [name, data] : this->fragments_->data_) {
		const int n_cell = data.size();
		for (int i = 0; i < n_cell; ++i) {
			const auto& cell_data = data[i];
			const std::size_t fragments_size = cell_data.first.size();
			for (size_t j = 0; j < fragments_size; ++j) {
				this->find_row(name, i, cell_data.first[j], cell_data.second[j]);
			}
		}
	}
	return true;
}
