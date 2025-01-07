#include "FragmentsQualityViewWorker.h"

#include "ItemDatabase.h"
#include "Custom.h"
#include "GenomeUtility.h"

FragmentsQualityViewWorker::FragmentsQualityViewWorker(const SingleCellMultiome* single_cell_multiome, const Fragments* fragments) :
	mode_(WorkMode::MultiomeFragmentsObject),
	single_cell_multiome_(single_cell_multiome),
	fragments_(fragments),
	species_(single_cell_multiome->species_)
{};

FragmentsQualityViewWorker::FragmentsQualityViewWorker(const SingleCellAtac* single_cell_atac, const Fragments* fragments) :
	mode_(WorkMode::AtacFragmentsObject),
	single_cell_atac_(single_cell_atac),
	fragments_(fragments),
	species_(single_cell_atac->species_)
{};

void FragmentsQualityViewWorker::create_index() {

	int n_cell = 0;

	if (this->mode_ == WorkMode::AtacFragmentsObject) {
		n_cell = this->single_cell_atac_->counts()->cols();
	}
	else {
		n_cell = this->single_cell_multiome_->atac_counts()->cols();
	}

	this->length_distribution_.resize(1000, 0);

	this->flank_counts_.resize(n_cell, 0);
	this->center_counts_.resize(n_cell, 0);
	this->blacklist_counts_.resize(n_cell, 0);
	this->not_blacklist_counts_.resize(n_cell, 0);
	this->mono_nucleosome_count_.resize(n_cell, 0);
	this->nucleosome_free_count_.resize(n_cell, 0);
	this->fragments_in_peak_.resize(n_cell, 0);
	this->fragments_not_in_peak_.resize(n_cell, 0);
	this->cell_index_.reserve(n_cell);
};

void FragmentsQualityViewWorker::get_results() {

	auto flank_mean = custom::cast<double>(custom::multiply(this->flank_counts_, 5));
	custom::assign(flank_mean, custom::mean(flank_mean), custom::equal(this->flank_counts_, 0));
	auto center_norm = custom::divide(this->center_counts_, flank_mean);

	this->res_[METADATA_TSS_ENRICHMENT] = center_norm;
	this->res_[METADATA_TSS_PERCENTILE] = custom::empirical_cumulative_distribution(center_norm);
	this->res_[METADATA_NUCLEOSOME_SIGNAL] = custom::partial_divide<double>(this->mono_nucleosome_count_, this->nucleosome_free_count_);
	this->res_[METADATA_BLACKLIST_RATIO] = custom::partial_divide<double>(this->blacklist_counts_, this->not_blacklist_counts_);
	this->res_[METADATA_FRIP_SCORE] = custom::partial_divide<double>(this->fragments_in_peak_, custom::add(this->fragments_in_peak_, this->fragments_not_in_peak_));
	this->res_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION] =	custom::partial_divide<double>(this->length_distribution_, std::ranges::max(this->length_distribution_));

};

bool FragmentsQualityViewWorker::create_blacklist_index(const QString& bed_file) {
	QFile file(bed_file);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		G_TASK_WARN("Unexpected error when reading blacklist file.");
		return false;
	}

	QTextStream in(&file);

	QString seq;
	int start{ 0 }, end{ 0 };
	GenomicRange blacklist_region;
	while (!in.atEnd()) {
		in >> seq >> start >> end;
		blacklist_region.append(custom::standardize_chromosome_name(seq), start, end - start, '*');
	}

	blacklist_region.finalize();

	if (blacklist_region.is_empty()) {
		return false;
	}

	this->blacklist_index_.set(blacklist_region);

	return true;
};

bool FragmentsQualityViewWorker::work() {

	if (this->species_ == soap::Species::Human) {

		bool success = ItemDatabase::read_item(FILE_HUMAN_GENOME_GENOMIC_RANGE_SIF, this->genome_);

		if (this->genome_.is_empty()) {
			G_TASK_WARN("Loading failed.");
			return false;
		}
	}
	else {
		G_TASK_WARN("Now only human genome is supported.");
		return false;
	}

	this->get_tss_position();

	GenomicRange upstream_flank = this->genome_.extended(1000, -900, true);
	GenomicRange downstream_flank = this->genome_.extended(-900, 1000, true);
	GenomicRange centers = this->genome_.extended(500, 500, true);

	upstream_flank.append(downstream_flank);

	this->flank_index_.set(upstream_flank);
	this->center_index_.set(centers);

	if (this->mode_ == WorkMode::AtacFragmentsObject) {
		this->peak_index_.set(this->single_cell_atac_->counts()->rownames_);
	}
	else {
		this->peak_index_.set(this->single_cell_multiome_->atac_counts()->rownames_);
	}

	this->get_blacklist_region();

	this->create_index();

	this->calculate_metric();

	this->get_results();

	return true;
};

void FragmentsQualityViewWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_qc_ready(this->res_);

	G_TASK_LOG("Fragments Quality Check Finished.");

	G_TASK_END;
}

bool FragmentsQualityViewWorker::calculate_metric() {
	const auto& cell_names = this->fragments_->cell_names_;
	for (const auto& [name, data] : this->fragments_->data_) {

		const int n_cell = data.size();
		for (int i = 0; i < n_cell; ++i) {
			const auto& cell_data = data[i];
			const std::size_t n_fragments = cell_data.first.size();

			const auto& start_data = cell_data.first;
			const auto& end_data = cell_data.second;

			for (std::size_t j = 0; j < n_fragments; ++j) {

				int start = start_data[j];
				int end = end_data[j];
				int width = end - start;

				bool is_nucleosome_free = width <= 147;
				bool is_mono_nucleosome = (width > 147) && (width <= 294);

				if (width < 1000 && width >= 0) {
					++length_distribution_[width];
				}

				if (is_mono_nucleosome) {
					++this->mono_nucleosome_count_[i];
				}
				if (is_nucleosome_free) {
					++this->nucleosome_free_count_[i];
				}

				int peak_count = this->peak_index_.count_insertion(name, start, end);
				this->fragments_in_peak_[i] += peak_count;
				this->fragments_not_in_peak_[i] += (2 - peak_count);

				int center_count = this->center_index_.count_insertion(name, start, end);
				if (center_count > 0) {
					this->center_counts_[i] += center_count;
				}

				int flank_count = this->flank_index_.count_insertion(name, start, end);
				if (flank_count > 0) {
					this->flank_counts_[i] += flank_count;
				}

				int blacklist_count = this->blacklist_index_.count_insertion(name, start, end);
				this->blacklist_counts_[i] += blacklist_count;
				this->not_blacklist_counts_[i] += (2 - blacklist_count);

			}
		}
	}

	return true;
};

bool FragmentsQualityViewWorker::get_blacklist_region() {
	if (this->species_ == soap::Species::Human) {
		return this->create_blacklist_index(FILE_HUMAN_GRCH38_BLACKLIST_BED);
	}
	else {
		G_TASK_WARN("Now only human genome is supported.");
		return false;
	}
};

void FragmentsQualityViewWorker::get_tss_position() {

	auto filter = custom::equal(this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_BIOTYPE), QString("protein_coding"));
	filter *= !custom::test(this->genome_.sequence_names_,
		[](const QString& seq_name)->bool {return seq_name == "chrM" || seq_name == "MT" || seq_name == "Mt"; });

	this->genome_.row_slice(filter);
	this->genome_.strand_.assign('+', this->genome_.strand_ == '*');

	const QStringList& gene_id = this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_ID);
	QStringList unique_gene_ids = custom::unique(gene_id);

	QVector<int> index_first = custom::index_of(unique_gene_ids, gene_id);
	QVector<int> index_last = custom::last_index_of(unique_gene_ids, gene_id);

	auto end = this->genome_.get_sequence_end();

	this->genome_.sequence_names_.reorder(index_first);
	this->genome_.strand_.reorder(index_first);
	this->genome_.ranges_.start_ = custom::ifelse<QVector<int>>(
		this->genome_.strand_ == '+',
		custom::reordered(this->genome_.ranges_.start_, index_first),
		custom::reordered(end, index_last)
		);
	this->genome_.ranges_.width_.fill(1);

	this->genome_.metadata_.row_reorder(index_first);
	this->genome_.metadata_.update(METADATA_GENOMIC_RANGE_GENE_NAME, custom::make_unique(this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME)));
};