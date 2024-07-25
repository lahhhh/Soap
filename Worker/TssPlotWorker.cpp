#include "TssPlotWorker.h"

#include <zlib.h>

#include "ItemDatabase.h"
#include "Custom.h"
#include "GenomeUtility.h"

TssPlotWorker::TssPlotWorker(soap::Species species, const Fragments* fragments) :
	fragments_(fragments),
	species_(species)
{};

void TssPlotWorker::run() {

	if (this->species_ == soap::Species::Human) {
		bool success = ItemDatabase::read_item(FILE_HUMAN_GENOME_GENOMIC_RANGE_GCS, this->genome_);

		if (!success) {
			G_TASK_WARN("Loading faied.");
			G_TASK_END;
		}
	}
	else {
		G_TASK_WARN("Now only human genome is supported.");
		G_TASK_END;
	}

	this->get_tss_position();

	this->create_index();

	this->calculate_tss();

	this->get_results();

	emit x_tss_ready(this->tss_matrix_, this->tss_vector_);
	G_TASK_END;
}

void TssPlotWorker::create_index() {

	this->tss_index_.set(this->genome_);

	int n_range = this->genome_.size();
	this->tss_matrix_.resize(n_range, 6000);
};

void TssPlotWorker::get_results() {

	Eigen::ArrayXd row_mean = this->tss_matrix_.rowwise().mean();
	auto order = _Cs order(row_mean);
	this->tss_matrix_ = this->tss_matrix_(order, Eigen::all);

	int n_gene = this->tss_matrix_.rows();

	int n_row = n_gene;
	int compress_fold = n_row / 500;

	if (compress_fold > 1) {
		this->tss_matrix_ = _Cs row_compressed(this->tss_matrix_, compress_fold);
	}

	this->tss_matrix_ = _Cs column_compressed(this->tss_matrix_, 10);

	this->tss_vector_ = _Cs cast<QVector>(this->tss_matrix_.colwise().sum());
	this->tss_vector_ = _Cs divide(this->tss_vector_, this->tss_vector_.first() + this->tss_vector_.last() / 2.0);
};

void TssPlotWorker::find_insertion_location_in_tss(
	const QString& seq_name,
	int start,
	int end)
{
	auto [pid1, poff1, pid2, poff2, nid1, noff1, nid2, noff2] = this->tss_index_.find_offset(seq_name, start, end);



	if (this->tss_index_.success(pid1)) {
		++this->tss_matrix_(pid1, poff1);
	}

	if (this->tss_index_.success(pid2)) {
		++this->tss_matrix_(pid2, poff2);
	}

	if (this->tss_index_.success(nid1)) {
		++this->tss_matrix_(nid1, noff1);
	}

	if (this->tss_index_.success(nid2)) {
		++this->tss_matrix_(nid2, noff2);
	}
};

void TssPlotWorker::calculate_tss() {
	const auto& cell_names = this->fragments_->cell_names_;

	for (const auto& [name, data] : this->fragments_->data_) {
		const int n_cell = data.size();

		for (int i = 0; i < n_cell; ++i) {
			const auto& cell_data = data[i];
			const std::size_t n_fragments = cell_data.first.size();

			const auto& start_data = cell_data.first;
			const auto& end_data = cell_data.second;

			for (std::size_t j = 0; j < n_fragments; ++j) {

				this->find_insertion_location_in_tss(name, start_data[j], end_data[j]);

			}
		}
	}
};

void TssPlotWorker::get_tss_position() {
	auto filter = _Cs equal(this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_BIOTYPE), QString("protein_coding"));
	filter *= !_Cs test(this->genome_.sequence_names_,
		[](const QString& seq_name)->bool {return seq_name == "chrM" || seq_name == "MT" || seq_name == "Mt"; });

	this->genome_.row_slice(filter);
	this->genome_.strand_.assign('+', this->genome_.strand_ == '*');

	const QStringList& gene_id = this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_ID);
	QStringList unique_gene_ids = _Cs unique(gene_id);

	QVector<int> index_first = _Cs index_of(unique_gene_ids, gene_id);
	QVector<int> index_last = _Cs last_index_of(unique_gene_ids, gene_id);

	auto end = this->genome_.get_sequence_end();

	this->genome_.sequence_names_.reorder(index_first);
	this->genome_.strand_.reorder(index_first);
	this->genome_.ranges_.start_ = _Cs ifelse<QVector<int>>(
		this->genome_.strand_ == '+',
		_Cs reordered(this->genome_.ranges_.start_, index_first),
		_Cs reordered(_Cs minus(end, 1), index_last)
		);
	this->genome_.ranges_.width_.resize(index_first.size());
	this->genome_.ranges_.width_.fill(1);

	this->genome_.metadata_.row_reorder(index_first);
	this->genome_.metadata_.update(METADATA_GENOMIC_RANGE_GENE_NAME,
		_Cs make_unique(this->genome_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME)));

	this->genome_.extend(2999, 3000);
};