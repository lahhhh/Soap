#include "MotifPosition.h"

MotifPosition::MotifPosition(
	const std::size_t n_peak, 
	const std::size_t n_motif)
{
	this->motif_information_.resize(n_motif);
};

bool MotifPosition::contains(const QString& motif_name) const {
	return this->motifs_.contains(motif_name) && this->motif_names_.contains(motif_name);
};

Eigen::MatrixXi MotifPosition::get_motif_matrix() const {

	Eigen::MatrixXi ret = Eigen::MatrixXi::Zero(this->n_peak(), this->n_motif());

	int n_motif = this->n_motif();
	for (int i = 0; i < n_motif; ++i) {
		const auto& motif_data = this->motif_information_[i];

		for (auto&& match : motif_data) {
			++ret(match.peak_index, i);
		}
	}

	return ret;
};

QStringList MotifPosition::get_match_peak_position(const QString& motif_name) const {

	QStringList ret;
	qsizetype index = this->motif_names_.indexOf(motif_name);
	if (index == -1 || !this->motifs_.contains(motif_name)) {
		return ret;
	}

	auto&& motif_information = this->motif_information_[index];

	for (auto&& match : motif_information) {
		ret << this->peak_names_[match.peak_index];
	}

	return ret;
};

GenomicRange MotifPosition::get_position(const QString& motif_name) const {

	GenomicRange ret;
	qsizetype index = this->motif_names_.indexOf(motif_name);
	if (index == -1 || !this->motifs_.contains(motif_name)) {
		return ret;
	}

	const int motif_width = this->motifs_.at(motif_name).weight_.mat_.cols();
	auto&& motif_information = this->motif_information_[index];
	const std::size_t match_size = motif_information.size();

	for (std::size_t i = 0; i < match_size; ++i) {
		const auto& match = motif_information[i];
		std::size_t peak_index = match.peak_index;
		auto [sequence_name, start, width, strand] = this->peak_locations_[peak_index];
		ret.append(sequence_name, start + match.position, motif_width, match.strand);
	}

	ret.finalize();

	return ret;
};

QStringList MotifPosition::peak_to_motif(int peak_index) const {

	QStringList ret;

	int n_motif = this->n_motif();
	for (int i = 0; i < n_motif; ++i) {

		for (auto&& match : this->motif_information_[i]) {
			if (match.peak_index == peak_index) {
				ret << this->motif_names_[i];
				break;
			}
		}
	}

	return ret;
};

qsizetype MotifPosition::n_motif() const {

	return this->motif_names_.size();
};

qsizetype MotifPosition::n_peak() const {

	return this->peak_names_.size();
};

QStringList MotifPosition::peak_to_tf(int peak_index) const {

	QStringList ret;

	for (int i = 0; i < this->n_motif(); ++i) {

		for (auto&& match : this->motif_information_[i]) {
			if (match.peak_index == peak_index) {
				ret << this->motif_names_[i];
				break;
			}
		}
	}

	return ret;
};

MotifPosition& MotifPosition::append(
	const std::size_t peak_location,
	const std::size_t motif_location,
	const std::size_t position,
	const float score,
	const char strand)
{

	if (this->motif_information_.size() <= motif_location) {
		this->motif_information_.resize(motif_location + 1);
	}

	this->motif_information_[motif_location].emplace_back(peak_location, position, score, strand);

	return *this;
};

QVector<int> MotifPosition::get_motif_match_times(const QVector<int>& peak_index) const {

	int n_motif = this->n_motif();
	QVector<int> matches(n_motif, 0);
	for (int i = 0; i < n_motif; ++i) {
		const auto& motif_data = this->motif_information_[i];

		for (auto&& match : motif_data) {
			if (peak_index.contains(match.peak_index)) {
				++matches[i];
			}
		}
	}

	return matches;
};

int MotifPosition::get_match_count(const QVector<int>& peak_index, const int motif_index) const {

	const auto& motif_information = this->motif_information_[motif_index];
	const std::size_t match_size = motif_information.size();
	std::size_t count = 0;
	for (std::size_t i = 0; i < match_size; ++i) {
		if (peak_index.contains((int)motif_information[i].peak_index)) {
			++count;
		}
	}
	return count;
};

Eigen::ArrayXi MotifPosition::get_match_count(const int motif_index) const {

	const int n_peak = this->n_peak();

	Eigen::ArrayXi res = Eigen::ArrayXi::Zero(n_peak);

	if (motif_index >= this->n_motif() || motif_index < 0) {
		return res;
	}

	const auto& motif_information = this->motif_information_[motif_index];

	for (auto&& m : motif_information) {
		++res[m.peak_index];
	}

	return res;
};

int MotifPosition::get_match_count(const int peak_index, const int motif_index) const {
	if (peak_index < 0 || peak_index >= this->n_peak() || motif_index >= this->n_motif() || motif_index < 0) {
		return 0;
	}
	const auto& motif_information = this->motif_information_[motif_index];
	const std::size_t match_size = motif_information.size();
	std::size_t count = 0;
	for (std::size_t i = 0; i < match_size; ++i) {
		if (motif_information[i].peak_index == peak_index) {
			++count;
		}
	}
	return count;
};

int MotifPosition::get_match_count(const QString& peak_name, const QString& motif_name) const {

	return this->get_match_count(this->peak_names_.indexOf(peak_name), this->motif_names_.indexOf(motif_name));
};