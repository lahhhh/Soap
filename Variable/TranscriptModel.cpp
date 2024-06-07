#include "TranscriptModel.h"

#include "Custom.h"

void TranscriptModel::finalize() {
	auto size = this->size();
	auto n_interval = size / 2;

	std::vector<int> interval_starts(n_interval);

	for (std::size_t i = 0; i < n_interval; ++i) {
		interval_starts[i] = this->loc_[i * 2];
	}

	const auto order = _Cs order(interval_starts);

	std::vector<int> new_loc(size);

	for (std::size_t i = 0; i < n_interval; ++i) {
		new_loc[i * 2] = this->loc_[order[i] * 2];
		new_loc[i * 2 + 1] = this->loc_[order[i] * 2 + 1];
	}

	this->loc_ = new_loc;

	this->start_ = this->loc_.front();
	this->end_ = this->loc_.back();
};

int TranscriptModel::match(const std::vector<int>& segments) const {
	if (this->loc_.empty() || segments.empty()) {
		return 0;
	}

	int i{ 0 };
	int loc_size = this->loc_.size();
	int segment_size = segments.size();
	bool unspliced{ false };

	for (int j = 0; j < segment_size; j += 2) {

		if (segments[j] < this->start_ || segments[j + 1] > this->end_) {
			return 0;
		}

		while (segments[j] > this->loc_[i]) {
			++i;
		}

		if (i % 2 == 0 || segments[j + 1] > this->loc_[i]) {
			unspliced = true;
		}

	}

	if (unspliced) {
		return 1;
	}
	else {
		return 2;
	}

};