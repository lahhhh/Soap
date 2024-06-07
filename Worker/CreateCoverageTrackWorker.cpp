#include "CreateCoverageTrackWorker.h"

#include <zlib.h>

#include "Custom.h"
#include "GenomeUtility.h"
#include "ItemDatabase.h"

void CreateCoverageTrackWorker::run() {
	this->build_index();

	if (!this->calculate_fragments_size()) {
		G_TASK_END;
	}

	if (!this->calculate_matrix()) {
		G_TASK_END;
	}
	this->create_track();
	if (!this->load_annotation()) {
		G_TASK_END;
	}

	emit x_coverage_track_ready(this->coverage_track_);

	G_TASK_END;
}

void CreateCoverageTrackWorker::build_index() {

	const CustomMatrix* metadata = &this->metadata_->mat_;

	if (metadata->data_type_.at(this->group_name_) == CustomMatrix::DataType::QStringFactor) {
		this->group_factors_ = metadata->string_factors_.at(this->group_name_);
	}
	else {
		this->group_factors_ = _Cs cast<QString>(metadata->integer_factors_.at(this->group_name_));
	}
	const int n_level = this->group_factors_.size();
	QStringList group = metadata->get_qstring(this->group_name_);
	this->group_distribution_ = _Cs table(group);
	const int n_cell = group.size();

	this->cell_index_.resize(n_cell, 0);
	for (int i = 0; i < n_cell; ++i) {
		this->cell_index_[i] = this->group_factors_.indexOf(group[i]);
	}
	this->cell_size_.resize(n_cell, 0.0);
};

bool CreateCoverageTrackWorker::load_annotation() {

	bool success = ItemDatabase::read_item(FILE_HUMAN_GENOME_GENOMIC_RANGE_GCS, this->coverage_track_->annotation_);

	if (!success) {
		G_TASK_WARN("Loading faied.");
		return false;
	}

	return true;
};

void CreateCoverageTrackWorker::create_track() {

	this->coverage_track_ = new CoverageTrack();
	this->coverage_track_->level_name_ = this->group_name_;
	this->coverage_track_->levels_ = this->group_factors_;

	QStringList sequence_names = this->temp_data_.keys();

	const int n_level = this->group_factors_.size();

	for (const auto& sequence : sequence_names) {
		auto& insertion_matrix = this->coverage_track_->insertion_matrix_[sequence];
		const auto& temp = this->temp_data_[sequence];
		insertion_matrix.resize(n_level);
		for (int i = 0; i < n_level; ++i) {
			insertion_matrix[i].set_from_pairs(temp[i]);
			insertion_matrix[i] *= 1.0 / this->group_distribution_[this->group_factors_[i]];
		}
		this->temp_data_.remove(sequence);
	}
};

bool CreateCoverageTrackWorker::calculate_matrix() {

	const int n_level = this->group_factors_.size();

	for (const auto& [name, data] : this->fragments_->data_) {

		auto& sequence_matrix = this->temp_data_[name];
		sequence_matrix.resize(n_level);

		const int n_cell = data.size();
		for (int i = 0; i < n_cell; ++i) {
			const auto& cell_data = data[i];
			const std::size_t n_fragments = cell_data.first.size();

			const auto& start_data = cell_data.first;
			const auto& end_data = cell_data.second;

			int row = this->cell_index_[i];
			float size_factor = this->cell_size_[i];

			for (std::size_t j = 0; j < n_fragments; ++j) {
				sequence_matrix[row].emplace_back((start_data[j] - 1) / 10, size_factor);
			}

			for (std::size_t j = 0; j < n_fragments; ++j) {
				sequence_matrix[row].emplace_back((end_data[j] - 1) / 10, size_factor);
			}
		}
	}
	return true;
};

bool CreateCoverageTrackWorker::calculate_fragments_size() {
	const int n_level = this->group_factors_.size();

	for (const auto& [name, data] : this->fragments_->data_) {

		const int n_cell = data.size();
		for (int i = 0; i < n_cell; ++i) {

			this->cell_size_[i] += data[i].first.size();
			this->cell_size_[i] += data[i].second.size();
		}
	}
	double mean_fragments_size = _Cs mean(this->cell_size_);

	for (auto&& size : this->cell_size_) {
		if (size != 0.0) {
			size = mean_fragments_size / size;
		}
	}

	return true;
}