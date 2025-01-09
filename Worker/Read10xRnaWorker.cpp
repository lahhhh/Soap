#include "Read10xRnaWorker.h"

#include <zlib.h>
#include <QFile>
#include "Custom.h"

bool Read10XRnaWorker::read_barcodes() {

	std::unique_ptr<char[]> bu(new char[1024]);

	char* buffer = bu.get();

	gzFile barcodes = gzopen_w((const wchar_t*)this->barcodes_file_name_.utf16(), "rb");

	if (barcodes == NULL) {
		G_TASK_WARN("Please check if the path contains illegal characters.");
		return false;
	}

	while ((gzgets(barcodes, buffer, 1024)) != 0)
	{
		this->barcodes_ << QString::fromUtf8(buffer, custom::line_length(buffer));
	}

	gzclose(barcodes);

	return true;
};

bool Read10XRnaWorker::read_features() {

	std::unique_ptr<char[]> bu(new char[1024]);

	char* buffer = bu.get();

	gzFile features = gzopen_w((const wchar_t*)this->features_file_name_.utf16(), "rb");

	if (features == NULL) {
		return false;
	}

	char* c;

	while ((gzgets(features, buffer, 1024)) != 0)
	{
		c = buffer;
		int tab_count{ 0 };
		while (*c != '\0') {
			if (*c == '\t') {
				++tab_count;
			}
			++c;
		}

		if (tab_count >= 2) {

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
		else { // someone edited the original feature file
			c = buffer;
			char* feature_start = c;
			while (*c != '\t' && *c != '\0') {
				++c;
			}
			if (c - feature_start < 2) {
				return false;
			}
			this->gene_symbols_ << QString::fromUtf8(feature_start, c - feature_start - 1);
		}
	}

	gzclose(features);

	return true;
};

void Read10XRnaWorker::determine_species() {
	for (auto&& i : this->gene_symbols_) {
		if (i.startsWith("MT-")) {
			this->res_->species_ = soap::Species::Human;
			break;
		}
		if (i.startsWith("mt-")) {
			this->res_->species_ = soap::Species::Mouse;
			break;
		}
	}
};

bool Read10XRnaWorker::read_matrix() {

	std::unique_ptr<char[]> bu(new char[1024]);

	char* buffer = bu.get();

	gzFile matrix = gzopen_w((const wchar_t*)this->matrix_file_name_.utf16(), "rb");

	if (matrix == NULL) {
		return false;
	}

	while (gzgets(matrix, buffer, 1024) != 0) {
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

	SparseInt& counts = SUBMODULES(*this->res_, SparseInt)[VARIABLE_COUNTS];
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

	return true;
};

void Read10XRnaWorker::calculate_metadata() {
	SparseInt& counts = SUBMODULES(*this->res_, SparseInt)[VARIABLE_COUNTS];
	
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
	if (this->res_->species_ == soap::Species::Human) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (this->res_->species_ == soap::Species::Mouse) {
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

	Metadata& metadata = SUBMODULES(*this->res_, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, QVector<int>(col_count.begin(), col_count.end()));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, QVector<int>(col_gene.begin(), col_gene.end()));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content)); 
	metadata.mat_.update(METADATA_BARCODES, this->barcodes_);
};

bool Read10XRnaWorker::work() {

	this->res_.reset(new SingleCellRna());

	try {
		if (!this->read_barcodes()) {
			G_TASK_WARN("Barcodes loading failed.");
			return false;
		}
	}
	catch (...) {
		G_TASK_WARN("Barcodes loading failed.");
		return false;
	}

	try {
		if (!this->read_features()) {
			G_TASK_WARN("Features loading failed.");
			return false;
		}
	}
	catch (...) {
		G_TASK_WARN("Features loading failed.");
		return false;
	}

	try {
		if (!this->read_matrix()) {
			G_TASK_WARN("Matrix loading failed.");
			return false;
		}
	}
	catch (...) {
		G_TASK_WARN("Matrix loading failed.");
		return false;
	}

	this->determine_species();

	this->calculate_metadata();

	return true;
};

void Read10XRnaWorker::run() {

	G_TASK_LOG("Start loading 10X scRNA...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_data_create_soon(this->res_.release(), soap::VariableType::SingleCellRna, "Single Cell RNA Data");

	G_TASK_LOG("Loading finished.");
	
	G_TASK_END;
}