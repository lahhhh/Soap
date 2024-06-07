#include "CoveragePlotWorker.h"

#include <zlib.h>

#include "Custom.h"
#include "GenomeUtility.h"
#include "ItemDatabase.h"

CoveragePlotWorker::CoveragePlotWorker(
	const SingleCellMultiome* single_cell_multiome,
	const Fragments* fragments,
	const QStringList& factors,
	const QStringList& levels,
	const QString& region,
	bool draw_gene,
	bool draw_link,
	double link_cutoff,
	bool draw_legend)
	:
	mode_(WorkMode::MultiomeFragmentsObject),
	single_cell_multiome_(single_cell_multiome),
	fragments_(fragments),
	factors_(factors),
	levels_(levels),
	region_(region),
	draw_gene_(draw_gene),
	draw_link_(draw_link),
	link_cutoff_(link_cutoff),
	draw_legend_(draw_legend)
{};

CoveragePlotWorker::CoveragePlotWorker(
	const SingleCellMultiome* single_cell_multiome,
	const Fragments* fragments,
	const QStringList& factors,
	const QStringList& levels,
	const QString& region,
	const QString& ccan_name,
	const QVector<std::pair<int, int>>& ccan_locs,
	bool draw_gene,
	bool draw_link,
	double link_cutoff,
	bool draw_legend)
	:
	mode_(WorkMode::MultiomeFragmentsObjectCcan),
	single_cell_multiome_(single_cell_multiome),
	fragments_(fragments),
	factors_(factors),
	levels_(levels),
	region_(region),
	ccan_name_(ccan_name),
	ccan_locs_(ccan_locs),
	draw_gene_(draw_gene),
	draw_link_(draw_link),
	link_cutoff_(link_cutoff),
	draw_legend_(draw_legend)
{};

CoveragePlotWorker::CoveragePlotWorker(
	const SingleCellAtac* single_cell_atac,
	const Fragments* fragments,
	const QStringList& factors,
	const QStringList& levels,
	const QString& region,
	bool draw_gene,
	bool draw_link,
	double link_cutoff,
	bool draw_legend)
	:
	mode_(WorkMode::AtacFragmentsObject),
	single_cell_atac_(single_cell_atac),
	fragments_(fragments),
	factors_(factors),
	levels_(levels),
	region_(region),
	draw_gene_(draw_gene),
	draw_link_(draw_link),
	link_cutoff_(link_cutoff),
	draw_legend_(draw_legend)
{};

CoveragePlotWorker::CoveragePlotWorker(
	const SingleCellAtac* single_cell_atac,
	const Fragments* fragments,
	const QStringList& factors,
	const QStringList& levels,
	const QString& region,
	const QString& ccan_name,
	const QVector<std::pair<int, int>>& ccan_locs,
	bool draw_gene,
	bool draw_link,
	double link_cutoff,
	bool draw_legend)
	:
	mode_(WorkMode::AtacFragmentsObjectCcan),
	single_cell_atac_(single_cell_atac),
	fragments_(fragments),
	factors_(factors),
	levels_(levels),
	region_(region),
	ccan_name_(ccan_name),
	ccan_locs_(ccan_locs),
	draw_gene_(draw_gene),
	draw_link_(draw_link),
	link_cutoff_(link_cutoff),
	draw_legend_(draw_legend)
{};

bool CoveragePlotWorker::find_peak_in_region() {
	QStringList peak_names;

	if (this->mode_ == WorkMode::MultiomeFragmentsObject ||
		this->mode_ == WorkMode::MultiomeFragmentsObjectCcan) {
		peak_names = this->single_cell_multiome_->atac_counts()->rownames_;
	}
	else {
		peak_names = this->single_cell_atac_->counts()->rownames_;
	}

	const qsizetype size = peak_names.size();
	for (qsizetype i = 0; i < size; ++i) {
		auto [sequence_name, start, end, success] = _Cs string_to_peak(peak_names[i]);
		if (success) {
			if (sequence_name == this->location_.sequence_name) {
				if (start < this->location_.end && end > this->location_.start) {
					this->peak_locations_ << std::make_pair(start, end);
					this->peak_index_ << i;
				}
			}
		}
	}
	return this->peak_locations_.size() != 0;
}

bool CoveragePlotWorker::find_link_in_region() {

	int n_peak = this->peak_index_.size();

	if (n_peak < 2) {
		return false;
	}

	const Cicero* cicero{ nullptr };
	if (this->mode_ == WorkMode::MultiomeFragmentsObject ||
		this->mode_ == WorkMode::MultiomeFragmentsObjectCcan) {
		cicero = this->single_cell_multiome_->cicero();
	}
	else {
		cicero = this->single_cell_atac_->cicero();
	}
	if (cicero == nullptr) {
		G_TASK_WARN("No Links Found.");
			return false;
	}

	std::ranges::sort(this->peak_index_);

	auto [min_peak_index, max_peak_index] = std::ranges::minmax(this->peak_index_);

	Eigen::MatrixXd links = Eigen::MatrixXd::Zero(n_peak, n_peak);

	auto&& connections = cicero->connections_;

	for (auto pd : this->peak_index_) {

		int p1 = this->peak_index_.indexOf(pd);

		for (Eigen::SparseMatrix<double>::InnerIterator it(connections, pd); it; ++it) {

			int row = it.row();

			if (row < min_peak_index) {
				continue;
			}

			if (row > max_peak_index) {
				break;
			}

			double val = it.value();

			if (val < this->link_cutoff_) {
				continue;
			}

			int p2 = this->peak_index_.indexOf(row);

			if (p2 != -1) {
				links(p1, p2) += val;
			}
		}
	}

	if (links.sum() == 0.0) {
		return false;
	}

	links += links.transpose().eval();
	links /= 2;

	for (int i = 0; i < n_peak; ++i) {
		for (int j = 0; j < i; ++j) {
			if (links(i, j) > 0) {
				this->peak_links_.emplace_back(i, j, links(i, j));
			}
		}
	}

	return true;
};

bool CoveragePlotWorker::find_gene_in_region() {

	auto&& genome_start = this->genome_.ranges_.start_;
	auto genome_end = this->genome_.get_sequence_end();
	QVector<char> strand = this->genome_.strand_.to_qvector();
	auto sequence_filter = this->genome_.sequence_names_ == this->location_.sequence_name;
	if (sequence_filter.count() == 0) {
		G_TASK_WARN("No Gene Found in region : " + this->location_.sequence_name + ":" + QString::number(this->location_.start) + "-" + QString::number(this->location_.end) + ".");
		return false;
	}

	auto start = _Cs sliced(genome_start, sequence_filter);
	auto end = _Cs sliced(genome_end, sequence_filter);
	strand = _Cs sliced(strand, sequence_filter);
	QStringList type = _Cs sliced(this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_TYPE), sequence_filter);
	QStringList gene_name = _Cs sliced(this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME), sequence_filter);

	const qsizetype sequence_size = start.size();
	int last_start = 0, last_end = 0;
	for (qsizetype i = 0; i < sequence_size; ++i) {
		int local_start = start[i], local_end = end[i];
		if (local_start < this->location_.end && local_end > this->location_.start) {
			if (local_start == last_start && local_end == last_end) {
				continue;
			}
			if (type[i] != "gap") {
				this->gene_structure_[gene_name[i]].first = strand[i];
			}
			last_start = local_start;
			last_end = local_end;
			this->gene_structure_[gene_name[i]].second.append(std::make_tuple(local_start, local_end, type[i]));
		}
	}
	if (_Cs sum(_Cs sapply(this->gene_structure_.values(), [](const QPair< char, QList< std::tuple<int, int, QString> > >& value) {return value.second.size(); })) == 0) {
		G_TASK_WARN("No Gene Found in region : " + this->location_.sequence_name + ":" + QString::number(this->location_.start) + "-" + QString::number(this->location_.end) + ".");
		return false;
	}
	for (const auto& name : this->gene_structure_.keys()) {
		if (this->gene_structure_[name].second.size() == 0) {
			this->gene_structure_.remove(name);
		}
	}
	return true;
}

bool CoveragePlotWorker::load_genome() {
	soap::Species species;

	if (this->mode_ == WorkMode::MultiomeFragmentsObject ||
		this->mode_ == WorkMode::MultiomeFragmentsObjectCcan) {
		species = this->single_cell_multiome_->species_;
	}
	else {
		species = this->single_cell_atac_->species_;
	}
	if (species == soap::Species::Human) {
		bool success = ItemDatabase::read_item(FILE_HUMAN_GENOME_GENOMIC_RANGE_GCS, this->genome_);

		if (!success) {
			G_TASK_WARN("Loading faied.");
			return false;
		}
		return true;
	}
	else {
		G_TASK_WARN("Coverage Plot Now Only Support Human Genome.");
		return false;
	}
};

void CoveragePlotWorker::smooth_matrix() {
	const Eigen::Index width = this->insertion_matrix_.cols();
	const Eigen::Index n_level = this->insertion_matrix_.rows();
	this->normalized_matrix_.resize(n_level, width);

	Eigen::MatrixXd average_matrix = this->insertion_matrix_ / 100;

	this->insertion_matrix_.resize(0, 0);

	for (Eigen::Index row = 0; row < n_level; ++row) {
		for (Eigen::Index col = 0; col < 51; ++col) {
			this->normalized_matrix_(row, col) = average_matrix.row(row).segment(0, col + 50).sum() / (col + 50) * 100;
		}
		double average = this->normalized_matrix_(row, 50);
		for (Eigen::Index col = 51; col < width - 49; ++col) {
			average += (average_matrix(row, col + 49) - average_matrix(row, col - 51));
			this->normalized_matrix_(row, col) = average;
		}
		for (Eigen::Index col = width - 49; col < width; ++col) {
			this->normalized_matrix_(row, col) = average_matrix.row(row).segment(col - 50, width - col + 50).sum() / (width - col + 50) * 100;
		}
	}
};

bool CoveragePlotWorker::calculate_matrix() {

	const int n_level = this->levels_.size();

	for (const auto& [name, data] : this->fragments_->data_) {

		if (name != this->location_.sequence_name) {
			continue;
		}

		const int n_cell = data.size();
		for (int i = 0; i < n_cell; ++i) {

			const auto& cell_data = data[i];
			const std::size_t n_fragments = cell_data.first.size();

			const auto& start_data = cell_data.first;
			const auto& end_data = cell_data.second;

			int row = this->cell_index_[i];

			if (row == -1) {
				continue;
			}

			double size_factor = this->cell_size_[i];

			auto lower_bound_it = std::ranges::lower_bound(start_data, this->location_.start);
			auto start_loc = std::distance(start_data.cbegin(), lower_bound_it);

			for (std::size_t j = start_loc; j < n_fragments; ++j) {
				int loc = start_data[j];

				if (loc >= this->location_.end) {
					break;
				}

				this->insertion_matrix_(row, loc - this->location_.start) += size_factor;
			}

			for (std::size_t j = start_loc; j < n_fragments; ++j) {
				int loc = end_data[j];

				if (loc >= this->location_.end) {
					break;
				}

				this->insertion_matrix_(row, loc - this->location_.start) += size_factor;
			}
		}
	}

	return true;
};

bool CoveragePlotWorker::prepare_draw_matrix() {

	const int n_level = this->levels_.size();
	for (int i = 0; i < n_level; ++i) {
		this->normalized_matrix_.row(i) /= (double)this->group_distribution_[this->levels_[i]];
	}

	double max_value = this->normalized_matrix_.maxCoeff();

	if (max_value == 0.0) {
		return false;
	}

	this->normalized_matrix_ /= max_value;
	return true;
};

bool CoveragePlotWorker::get_location_by_name() {

	auto [seq_name, start, end, strand, success] =
		_Cs find_gene_in_genome(this->region_, this->genome_);
	if (!success) {
		return false;
	}
	this->location_.sequence_name = seq_name;
	this->location_.start = start;
	this->location_.end = end;
	this->get_location_by_gene_name_ = true;
	return true;
};

void CoveragePlotWorker::expand_plot_region() {
	int extend_length = 0;
	if (this->get_location_by_gene_name_) {
		extend_length = (this->location_.end - this->location_.start) / 6;
		if (extend_length < 100) {
			extend_length = 100;
		}
	}
	else if (this->location_.end - this->location_.start < 200) {
		G_TASK_NOTICE("Given region is too narrow for visualization. Region is expanded.");
		extend_length = 100;
	}
	else {
		return;
	}
	this->location_.start -= extend_length;
	if (this->location_.start < 1) {
		this->location_.start = 1;
	}
	this->location_.end += extend_length;
};

void CoveragePlotWorker::build_index() {

	const int n_level = this->levels_.size();
	this->group_distribution_ = _Cs table(this->factors_);
	const int n_cell = this->factors_.size();

	const CustomMatrix* metadata{ nullptr };
	if (this->mode_ == WorkMode::MultiomeFragmentsObject ||
		this->mode_ == WorkMode::MultiomeFragmentsObjectCcan) {
		metadata = &this->single_cell_multiome_->metadata()->mat_;
	}
	else {
		metadata = &this->single_cell_atac_->metadata()->mat_;
	}

	this->cell_index_.resize(n_cell, 0);
	for (int i = 0; i < n_cell; ++i) {
		this->cell_index_[i] = this->levels_.indexOf(this->factors_[i]);
	}
	this->cell_size_.resize(n_cell, 0.0);


	const int width = this->location_.end - this->location_.start;
	this->insertion_matrix_ = Eigen::MatrixXd::Zero(n_level, width);
};

bool CoveragePlotWorker::calculate_fragments_size() {
	const int n_level = this->levels_.size();

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

void CoveragePlotWorker::run() {

	if (!this->load_genome()) {
		G_TASK_END;
	}

	if (!this->get_location()) {
		G_TASK_WARN("Can not determine the plot region: " + this->region_);
		G_TASK_END;
	}
	this->expand_plot_region();

	if (this->location_.end - this->location_.start > 1e7) {
		G_TASK_WARN("Region are two broad for visualization.");
		G_TASK_END;
	}

	this->build_index();

	if (!this->calculate_fragments_size()) {
		G_TASK_END;
	}
	if (!this->calculate_matrix()) {
		G_TASK_END;
	}
	this->smooth_matrix();


	if (this->draw_gene_) {
		this->find_gene_in_region();
	}

	this->find_peak_in_region();

	if (this->draw_link_) {
		this->find_link_in_region();
	}

	if (!this->prepare_draw_matrix()) {
		G_TASK_WARN("No insertion found in region.");
		G_TASK_END;
	}

	COVERAGE_PLOT_ELEMENTS res;
	res.region_name = this->region_;
	res.region = this->location_;
	res.group_factors = this->levels_;
	res.normalized_matrix = this->normalized_matrix_;
	res.gene_structure = this->gene_structure_;
	res.peak_locations = this->peak_locations_;
	res.peak_links = this->peak_links_;
	res.draw_legend = this->draw_legend_;
	res.ccan_locations = this->ccan_locs_;

	if (this->mode_ == WorkMode::MultiomeFragmentsObjectCcan || this->mode_ == WorkMode::AtacFragmentsObjectCcan) {
		res.plot_title = this->ccan_name_;
	}
	else {
		res.plot_title = this->region_;
	}

	emit x_plot_ready(res);
	G_TASK_END;
}

bool CoveragePlotWorker::get_location() {
	auto [sequence_name, start, end, success] = _Cs string_to_peak(this->region_);
	if (!success) {
		return get_location_by_name();
	}
	else {
		this->location_.sequence_name = sequence_name;
		this->location_.start = start;
		this->location_.end = end;
		return true;
	}
};