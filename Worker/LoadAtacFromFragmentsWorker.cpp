#include "LoadAtacFromFragmentsWorker.h"

#include "MacsCallPeakWorker.h"
#include "CalculateCountsByGenomicRangeWorker.h"

#include <zlib.h>

bool LoadAtacFromFragmentsWorker::read_fragments() {

	this->single_cell_atac_.reset(new SingleCellAtac());
	auto fragments = &SUBMODULES(*this->single_cell_atac_, Fragments)[VARIABLE_FRAGMENTS];

	gzFile file = gzopen(this->file_name_.toUtf8().data(), "rb");

	if (file == NULL) {
		G_TASK_WARN("File : " + this->file_name_ + " is broken.");
		return false;
	}

	int buffer_length = 256;
	char buffer[256];

	std::size_t process = 0;

	bool not_end_of_file = false;
	while (not_end_of_file = gzgets(file, buffer, buffer_length) != 0)
	{
		if (buffer[0] != '#')break;
	}
	if (!not_end_of_file) {
		gzclose(file);
		G_TASK_WARN("File : " + this->file_name_ + " is broken.");
		return false;
	}

	QString barcode;

	do {
		const char* c = buffer, * d, * seq_end;
		while (*c != '\t') {
			++c;
		}
		seq_end = c;

		int start = custom::atoi_specialized(&c);
		int end = custom::atoi_specialized(&c);

		if (start >= end) {
			G_TASK_WARN("Invalid Fragments File: " + QString::fromUtf8(buffer));
			gzclose(file);
			return false;
		}

		++c;
		d = c;
		while (*c != '\t') {
			++c;
		}
		++c;
		barcode = QString::fromUtf8(d, c - d - 1);

		QString seq_name = custom::standardize_chromosome_name(QString::fromUtf8(buffer, seq_end - buffer));
		auto& sequence_matrix = fragments->data_[seq_name];

		auto iter = this->barcode_index_.find(barcode);
		if (iter != this->barcode_index_.end()) {

			int index = iter->second;

			if (sequence_matrix.empty()) {
				sequence_matrix.resize(this->n_barcode_);
			}

			auto& barcode_fragments = sequence_matrix[index];
			barcode_fragments.first.emplace_back(start);
			barcode_fragments.second.emplace_back(end);
		}
		else {

			int index = this->n_barcode_++;

			this->barcode_index_[barcode] = index;

			sequence_matrix.resize(this->n_barcode_);

			auto& barcode_fragments = sequence_matrix[index];
			barcode_fragments.first.emplace_back(start);
			barcode_fragments.second.emplace_back(end);
		}

		if (process % 50000000 == 0) {
			G_TASK_LOG(QString::number(process) + " sequences has been processed.");
		}

		++process;
	} while (not_end_of_file = gzgets(file, buffer, buffer_length) != 0);

	gzclose(file);

	QList<std::size_t> insertion_size(this->n_barcode_, 0);

	for (auto&& [name, insertion] : fragments->data_) {
		if (insertion.size() < this->n_barcode_) {
			insertion.resize(this->n_barcode_);
		}

		auto chr_insertion_size = custom::sapply(insertion, [](auto&& i) {return i.first.size(); });

		insertion_size = custom::add(insertion_size, chr_insertion_size);
	}

	fragments->cell_names_.resize(this->n_barcode_);
	for (auto&& [name, ind] : this->barcode_index_) {
		fragments->cell_names_[ind] = name;
	}

	auto filter = custom::sapply(insertion_size, [](auto d) {return d > 1000; });

	fragments->slice(filter);

	return true;
};

bool LoadAtacFromFragmentsWorker::call_peak() {

	this->peaks_ = MacsCallPeakWorker::call_peak({this->single_cell_atac_->fragments()});

	return !this->peaks_.is_empty();
};

bool LoadAtacFromFragmentsWorker::calculate_matrix() {

	auto counts = CalculateCountsByGenomicRangeWorker::calculate_counts(this->single_cell_atac_->fragments(), this->peaks_);

	if (counts == nullptr) {
		return false;
	}

	SUBMODULES(*this->single_cell_atac_, SparseInt)[VARIABLE_COUNTS] = std::move(*counts);
	delete counts;

	return true;
};

void LoadAtacFromFragmentsWorker::determine_species() {

	if (this->single_cell_atac_->fragments()->data_.contains("20")) {
		this->single_cell_atac_->species_ = soap::Species::Human;
	}
	else {
		this->single_cell_atac_->species_ = soap::Species::Mouse;
	}
};

void LoadAtacFromFragmentsWorker::calculate_metadata() {

	SparseInt& counts = *this->single_cell_atac_->counts();

	counts.data_type_ = SparseInt::DataType::Counts;

	Eigen::ArrayXi col_count = custom::col_sum_mt(counts.mat_);
	const int ncol = counts.mat_.cols();
	Eigen::ArrayXi col_peak(ncol);

	for (int i = 0; i < ncol; ++i) {
		col_peak[i] = counts.mat_.outerIndexPtr()[i + 1] - counts.mat_.outerIndexPtr()[i];
	}

	Metadata& metadata = SUBMODULES(*this->single_cell_atac_, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);
	metadata.mat_.update(METADATA_ATAC_UMI_NUMBER, QVector<int>(col_count.begin(), col_count.end()));
	metadata.mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, QVector<int>(col_peak.begin(), col_peak.end()));
	metadata.mat_.update(METADATA_BARCODES, this->single_cell_atac_->fragments()->cell_names_);
};

void LoadAtacFromFragmentsWorker::run() {

	if (!this->read_fragments()) {
		G_TASK_END;
	}

	if (!this->call_peak()) {
		G_TASK_WARN("Peak Calling Failed.");
		G_TASK_END;
	}

	if (!this->calculate_matrix()) {
		G_TASK_WARN("Counts Calculation Failed.");
		G_TASK_END;
	}

	this->determine_species();

	this->calculate_metadata();

	emit x_data_create_soon(this->single_cell_atac_.release(), soap::VariableType::SingleCellAtac, "Single Cell ATAC Data");

};