#include "TranscriptionalFactorFootprintingWorker.h"

#include <QtConcurrent>

#include "GenomeUtility.h"
#include "Custom.h"

TranscriptionalFactorFootprintingWorker::TranscriptionalFactorFootprintingWorker(
	const Metadata* metadata,
	std::pair<std::vector<std::string>, std::vector<double>> bias,
	soap::Species species,
	const MotifPosition* motif_position,
	const QStringList& transcriptional_factor_names,
	const Fragments* fragments
) :
	mode_(WorkMode::Object),
	metadata_(metadata),
	bias_(bias),
	species_(species),
	motif_position_(motif_position),
	transcriptional_factor_names_(transcriptional_factor_names),
	fragments_object_(fragments)
{};

TranscriptionalFactorFootprintingWorker::TranscriptionalFactorFootprintingWorker(
	const Metadata* metadata,
	std::pair<std::vector<std::string>, std::vector<double>> bias,
	soap::Species species,
	const GenomicRange& genomic_range,
	const QString& file_name,
	const Fragments* fragments,
	const QString& factor_name,
	const QStringList& levels,
	const QStringList& factors,
	const GraphSettings& graph_settings,
	int height,
	int width
) :
	mode_(WorkMode::GenomicRange),
	metadata_(metadata),
	bias_(bias),
	species_(species),
	designated_location_(genomic_range),
	picture_name_(file_name),
	fragments_object_(fragments),
	factor_name_(factor_name),
	levels_(levels),
	factors_(factors),
	graph_settings_(graph_settings),
	height_(height),
	width_(width)
{};

TranscriptionalFactorFootprintingWorker::TranscriptionalFactorFootprintingWorker(
	const Metadata* metadata,
	std::pair<std::vector<std::string>, std::vector<double>> bias,
	soap::Species species,
	const MotifPosition* motif_position,
	const QStringList& transcriptional_factor_names,
	const Fragments* fragments,
	const QString& output_directory,
	const QString& factor_name,
	const QStringList& levels,
	const QStringList& factors,
	const GraphSettings& graph_settings,
	int height,
	int width
) :
	mode_(WorkMode::Batch),
	metadata_(metadata),
	bias_(bias),
	species_(species),
	motif_position_(motif_position),
	transcriptional_factor_names_(transcriptional_factor_names),
	output_directory_(output_directory),
	factor_name_(factor_name),
	fragments_object_(fragments),
	levels_(levels),
	factors_(factors),
	graph_settings_(graph_settings),
	height_(height),
	width_(width)
{};

std::tuple<bool, std::vector<std::string>, GenomicRange>
TranscriptionalFactorFootprintingWorker::get_motif_location_sequence(const QString& tf_name) {

	GenomicRange motif_location = this->motif_position_->get_position(tf_name);

	int n_motif = motif_location.rows();

	motif_location.extend(this->upstream_, this->downstream_, false);

	auto motif_sequence0 = this->genome_.get_std_sequence(motif_location);
	const std::size_t sequence_size = motif_sequence0.size();

	std::vector<std::string> motif_sequence;
	motif_sequence.reserve(sequence_size);

	for (std::size_t i = 0; i < sequence_size; ++i) {
		if (!motif_sequence0[i].empty()) {
			motif_sequence.push_back(motif_sequence0[i]);
			motif_location.append(motif_location[i]);
		}
	}

	if (motif_sequence.empty()) {
		return std::make_tuple(false, motif_sequence, motif_location);
	}
	else {

		motif_location.finalize();

		return std::make_tuple(true, motif_sequence, motif_location);
	}
};

// ignore possibility that motif locations overlapped because results will not change a lot
void TranscriptionalFactorFootprintingWorker::count_position(
	const GenomeIndex<IndexMode::Strand>& motif_index,
	Eigen::MatrixXi& insertion_matrix,
	const QString& sequence_name,
	const int start,
	const int end,
	const int cell_index) const
{

	auto [f1, f2, b1, b2] = motif_index.find_offset(sequence_name, start, end);

	if (motif_index.success(f1)) {

		++insertion_matrix(cell_index, f1);
	}

	if (motif_index.success(f2)) {

		++insertion_matrix(cell_index, f2);
	}

	if (motif_index.success(b1)) {

		++insertion_matrix(cell_index, b1);
	}

	if (motif_index.success(b2)) {

		++insertion_matrix(cell_index, b2);
	}

};

std::pair<bool, Eigen::MatrixXi> TranscriptionalFactorFootprintingWorker::get_insertion_matrix(
	const GenomeIndex<IndexMode::Strand>& motif_index, int insertion_width) const {

	const int n_cell = this->fragments_object_->cell_names_.size();
	Eigen::MatrixXi insertion_matrix = Eigen::MatrixXi::Zero(n_cell, insertion_width);

	Eigen::ArrayX<bool> cell_filter = Eigen::ArrayX<bool>::Constant(n_cell, false);
	bool filter_cell{ false };
	if (this->mode_ == WorkMode::Batch) {
		for (auto&& level : this->levels_) {
			cell_filter += custom::equal(this->factors_, level);
		}

		filter_cell = true;
	}

	QVector<std::vector<std::pair<std::vector<int>, std::vector<int>>> const*> chr_data_list;
	QVector<QString> chr_name_list;
	for (const auto& [name, data] : this->fragments_object_->data_) {
		chr_name_list << name;
		chr_data_list << &data;
	}

	int n_chr = chr_data_list.size();

#pragma omp parallel for
	for (int k = 0; k < n_chr; ++k) {

		Eigen::MatrixXi insertion_matrix_chr = Eigen::MatrixXi::Zero(n_cell, insertion_width);

		auto&& data = *chr_data_list[k];
		QString name = chr_name_list[k];

		for (int i = 0; i < n_cell; ++i) {

			if (filter_cell) {
				if (!cell_filter[i]) {
					continue;
				}
			}

			const auto& cell_data = data[i];
			const std::size_t n_fragments = cell_data.first.size();

			const auto& start_data = cell_data.first;
			const auto& end_data = cell_data.second;

			for (std::size_t j = 0; j < n_fragments; ++j) {

				this->count_position(motif_index, insertion_matrix_chr, name, start_data[j], end_data[j], i);
			}
		}

	#pragma omp critical
		{
			insertion_matrix += insertion_matrix_chr;
		}
	}

	return std::make_pair(true, insertion_matrix);
}

QVector<double>
TranscriptionalFactorFootprintingWorker::find_expected_insertions(const std::vector<std::string>& motif_sequence) const {
	int insertion_width = motif_sequence[0].size() - 6;
	Eigen::MatrixXd frequency = Eigen::MatrixXd::Zero(4096, insertion_width);
	const std::size_t sequence_size = motif_sequence.size();
	for (std::size_t i = 0; i < sequence_size; ++i) {
		const char* data = motif_sequence[i].data();
		for (std::size_t j = 0; j < insertion_width; ++j) {
			int loc = custom::get_hexamer_loc(data + j);
			if (loc > -1) {
				++frequency(loc, j);
			}
		}
	}
	auto bias = custom::cast<Eigen::ArrayX>(this->bias_.second);

	QVector<double> expected_insertions(insertion_width);

	for (int i = 0; i < insertion_width; ++i) {
		expected_insertions[i] = (frequency.col(i).array() * bias).sum();
	}
	double flank = (custom::mean(expected_insertions.sliced(0, 50)) + custom::mean(expected_insertions.sliced(insertion_width - 50))) / 2;
	expected_insertions = custom::divide(expected_insertions, flank);

	return expected_insertions;
};

bool TranscriptionalFactorFootprintingWorker::calculate_insertion_bias() {

	auto& bias = this->bias_;
	if (bias.first.size() == 4096 && bias.second.size() == 4096) {
		return true;
	}

	std::string	sequence = this->genome_.get_std_whole_sequence("1");

	auto genome_frequency = custom::get_6_nucleotide_frequency(sequence);

	auto insertion_frequency = custom::get_6_nucleotide_insertion_frequency("1", sequence, *this->fragments_object_);

	if (insertion_frequency.first.empty()) {
		G_TASK_WARN("Insertion Bias Calculation Failed.");
		return false;
	}

	auto insertion_bias = custom::partial_divide<double>(insertion_frequency.second, genome_frequency.second);
	bias.first = insertion_frequency.first;
	bias.second = insertion_bias;

	return true;
}

void TranscriptionalFactorFootprintingWorker::create_index() {
	this->cell_names_ = this->fragments_object_->cell_names_;

};

bool TranscriptionalFactorFootprintingWorker::load_genome() {

	if (this->species_ == soap::Species::Human) {

		this->genome_.set_sequence_file(FILE_HUMAN_GRCH38_2BIT);
	}
	else {

		G_TASK_WARN("Unsupported species for footprinting.");

		return false;
	}

	return true;
};

void TranscriptionalFactorFootprintingWorker::mode2() {

	if (!this->load_genome()) {
		return;
	}

	this->create_index();

	if (!this->calculate_insertion_bias()) {
		return;
	};


	GenomicRange motif_location = this->designated_location_;

	int n_motif = motif_location.rows();

	motif_location.extend(this->upstream_, this->downstream_, false);

	auto motif_sequence0 = this->genome_.get_std_sequence(motif_location);
	const std::size_t sequence_size = motif_sequence0.size();

	std::vector<std::string> motif_sequence;
	motif_sequence.reserve(sequence_size);

	for (std::size_t i = 0; i < sequence_size; ++i) {
		if (!motif_sequence0[i].empty()) {
			motif_sequence.push_back(motif_sequence0[i]);
			motif_location.append(motif_location[i]);
		}
	}

	if (motif_sequence.empty()) {

		G_TASK_WARN("Found no legal sequence");

		return;
	}
	else {

		motif_location.finalize();
	}

	auto expected_insertions = this->find_expected_insertions(motif_sequence);

	GenomeIndex<IndexMode::Strand> motif_index(motif_location.extended(-3, -3));
	int insertion_width = motif_sequence[0].size() - 6;

	auto [success2, insertion_matrix] = this->get_insertion_matrix(motif_index, insertion_width);
	if (!success2) {
		return;
	}

	Footprint fpt(
		PatternWeightMatrix{},
		insertion_matrix,
		this->cell_names_,
		custom::cast<QString>(custom::minus(custom::seq_n(1, insertion_width), insertion_width / 2)),
		expected_insertions,
		motif_location
	);

	bool success = fpt.draw(
		this->picture_name_,
		this->factor_name_,
		this->levels_,
		this->factors_,
		this->graph_settings_,
		this->height_,
		this->width_
	);

	G_TASK_LOG("Task Finished.");
};

void TranscriptionalFactorFootprintingWorker::run() {

	if (this->mode_ == WorkMode::GenomicRange) {
		this->mode2();
		G_TASK_END;
		return;
	}

	if (!this->load_genome()) {
		G_TASK_END;
	}

	this->create_index();

	if (!this->calculate_insertion_bias()) {
		G_TASK_END;
	};

	int n_tf = this->transcriptional_factor_names_.size();

	for (int i = 0; i < n_tf; ++i) {

		QString tf_name = this->transcriptional_factor_names_[i];

		auto [success, motif_sequence, motif_location] = this->get_motif_location_sequence(tf_name);

		if (!success) {
			G_TASK_WARN("Found no legal sequence in " + tf_name);
			continue;
		}
		auto expected_insertions = this->find_expected_insertions(motif_sequence);

		GenomeIndex<IndexMode::Strand> motif_index(motif_location.extended(-3, -3));
		int insertion_width = motif_sequence[0].size() - 6;

		auto [success2, insertion_matrix] = this->get_insertion_matrix(motif_index, insertion_width);
		if (!success2) {
			continue;
		}

		if (this->mode_ == WorkMode::Batch) {

			Footprint fpt(
				this->motif_position_->motifs_.at(tf_name),
				insertion_matrix,
				this->cell_names_,
				custom::cast<QString>(custom::minus(custom::seq_n(1, insertion_width), insertion_width / 2)),
				expected_insertions,
				motif_location
			);

			QString file_name = this->output_directory_ + "/" + custom::standardize_windows_file_name(fpt.motif_.motif_name_) + ".png";

			bool success = fpt.draw(
				file_name,
				this->factor_name_,
				this->levels_,
				this->factors_,
				this->graph_settings_,
				this->height_,
				this->width_
			);

			if (!success) {
				G_TASK_WARN(tf_name + " task is failed.");
			}
		}
		else {

			G_TASK_LOG("Calculation of " + tf_name + " finished.");

			emit x_footprint_ready(
				Footprint(
					this->motif_position_->motifs_.at(tf_name),
					insertion_matrix,
					this->cell_names_,
					custom::cast<QString>(custom::minus(custom::seq_n(1, insertion_width), insertion_width / 2)),
					expected_insertions,
					motif_location
				)
			);
		}
	}

	G_TASK_LOG("Footprinting Finished.");
	G_TASK_END;
}