#include "Read10xRnaWorker.h"

#include <zlib.h>
#include <QFile>
#include "Custom.h"

bool Read10XRnaWorker::read_barcodes() {

	QString barcodes_file_path = this->path_10X_ + "/" + BARCODES_FILE_NAME_10X;
	auto file = barcodes_file_path.toUtf8();

	char* buffer = new char[256];

	gzFile barcodes = gzopen(file.data(), "rb");

	if (barcodes == NULL) {
		G_TASK_WARN("Please check if the path contains illegal characters.");
		return false;
	}

	while ((gzgets(barcodes, buffer, 256)) != 0)
	{
		this->barcodes_ << QString::fromUtf8(buffer, custom::line_length(buffer));
	}
	gzclose(barcodes);

	delete[] buffer;

	return true;
};

bool Read10XRnaWorker::read_features() {

	QString features_file_path = this->path_10X_ + "/" + FEATURE_FILE_NAME_10X;
	auto file = features_file_path.toUtf8();

	char* buffer = new char[1024];

	gzFile features = gzopen(file.data(), "rb");

	if (features == NULL) {
		return false;
	}

	char* c;

	while ((gzgets(features, buffer, 1024)) != 0)
	{
		c = buffer;
		while (*c != '\t') {
			++c;
		}
		++c;
		char* feature_start = c;
		while (*c != '\t') {
			++c;
		}
		this->gene_symbols_ << QString::fromUtf8(feature_start, c - feature_start);
	}
	gzclose(features);

	delete[] buffer;

	return true;
};

void Read10XRnaWorker::determine_species() {
	for (auto&& i : this->gene_symbols_) {
		if (i.startsWith("MT-")) {
			this->single_cell_rna_->species_ = soap::Species::Human;
			break;
		}
		if (i.startsWith("mt-")) {
			this->single_cell_rna_->species_ = soap::Species::Mouse;
			break;
		}
	}
};

bool Read10XRnaWorker::read_matrix() {
	QString matrix_file_path = this->path_10X_ + "/" + MATRIX_FILE_NAME_10X;
	auto file = matrix_file_path.toUtf8();

	char* buffer = new char[256];

	gzFile matrix = gzopen(file.data(), "rb");

	if (matrix == NULL) {
		return false;
	}

	while (gzgets(matrix, buffer, 256) != 0) {
		if (buffer[0] != '%') {
			break;
		}
	}

	const char* c = buffer;

	int n_row = custom::atoi_specialized(&c);
	int n_column = custom::atoi_specialized(&c);
	int n_counts = custom::atoi_specialized(&c);

	if (n_row != this->gene_symbols_.size() || n_column != this->barcodes_.length()) {
		return false;
	}

	SparseInt& counts = SUBMODULES(*this->single_cell_rna_, SparseInt)[VARIABLE_COUNTS];
	counts.data_type_ = SparseInt::DataType::Counts;

	Eigen::SparseMatrix<int, Eigen::RowMajor> trmat(n_row, n_column);
	typename Eigen::SparseMatrix<int>::IndexVector wi(n_row);
	wi.setZero();

	int* data = new int[n_counts * 3];
	for (int i = 0; i < n_counts; ++i) {
		gzgets(matrix, buffer, 256);
		c = buffer;
		data[3 * i] = custom::atoi_specialized(&c) - 1;
		data[3 * i + 1] = custom::atoi_specialized(&c) - 1;
		data[3 * i + 2] = custom::atoi_specialized(&c);
		++wi(data[3 * i]);
	}

	gzclose(matrix);
	trmat.reserve(wi);
	for (int i = 0; i < n_counts; ++i) {
		trmat.insertBackUncompressed(data[3 * i], data[3 * i + 1]) = data[3 * i + 2];
	}
	delete[]data;
	trmat.collapseDuplicates(Eigen::internal::scalar_sum_op<int, int>());

	counts.mat_ = trmat;
	counts.rownames_ = custom::make_unique(this->gene_symbols_);
	counts.colnames_ = custom::make_unique(this->barcodes_);
	delete[] buffer;

	return true;
};

void Read10XRnaWorker::calculate_metadata() {
	SparseInt& counts = SUBMODULES(*this->single_cell_rna_, SparseInt)[VARIABLE_COUNTS];
	
	Eigen::ArrayX<bool> gene_detected = custom::row_sum(counts.mat_) > 0;
	counts.row_slice(gene_detected);

	Eigen::ArrayXi col_count = custom::col_sum_mt(counts.mat_);
	const int ncol = counts.mat_.cols();
	Eigen::ArrayXi col_gene(ncol);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(ncol);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(ncol);

	for (int i = 0; i < ncol; ++i) {
		col_gene[i] = counts.mat_.outerIndexPtr()[i + 1] - counts.mat_.outerIndexPtr()[i];
	}

	QList<int> mitochondrial_location, ribosomal_location;

	auto& gene_symbols = counts.rownames_;
	if (this->single_cell_rna_->species_ == soap::Species::Human) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (this->single_cell_rna_->species_ == soap::Species::Mouse) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("mt-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("Rps") || gene_symbols.at(i).startsWith("Rpl")) {
				ribosomal_location << i;
			}
		}
	}

	if (mitochondrial_location.length() > 0) {
		for (int& i : mitochondrial_location) {
			mitochondrial_content += counts.mat_.row(i).cast<double>();
		}
		mitochondrial_content /= col_count.cast<double>();
	}

	if (ribosomal_location.length() > 0) {
		for (int& i : ribosomal_location) {
			ribosomal_content += counts.mat_.row(i).cast<double>();
		}
		ribosomal_content /= col_count.cast<double>();
	}

	custom::remove_na(mitochondrial_content);
	custom::remove_na(ribosomal_content);

	Metadata& metadata = SUBMODULES(*this->single_cell_rna_, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, QVector<int>(col_count.begin(), col_count.end()));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, QVector<int>(col_gene.begin(), col_gene.end()));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content)); 
	metadata.mat_.update(METADATA_BARCODES, this->barcodes_);
};

void Read10XRnaWorker::run() {

	this->single_cell_rna_ = new SingleCellRna();

	if (!this->read_barcodes()) {
		delete this->single_cell_rna_;
		G_TASK_WARN("Barcodes loading failed.");
		G_TASK_END;
	}

	if (!this->read_features()) {
		delete this->single_cell_rna_;
		G_TASK_WARN("Features loading failed.");
		G_TASK_END;
	}

	if (!this->read_matrix()) {
		delete this->single_cell_rna_;
		G_TASK_WARN("Matrix loading failed.");
		G_TASK_END;
	}

	this->determine_species();
	this->calculate_metadata();

	emit x_data_create_soon(this->single_cell_rna_, soap::VariableType::SingleCellRna, "Single Cell RNA Data");
	G_TASK_END;
}