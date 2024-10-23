#include "Read10xRnaAndAtacWorker.h"
#include <zlib.h>
#include "Custom.h"
#include "Identifier.h"

void Read10XMultiomeWorker::run() {
	this->single_cell_multiome_.reset(new SingleCellMultiome());

	this->single_cell_multiome_->create_field("RNA", DataField::DataType::Rna);
	this->single_cell_multiome_->create_field("ATAC", DataField::DataType::Atac);

	if (!this->read_barcodes()) {
		G_TASK_WARN("Barcodes loading failed.");
		G_TASK_END;
	}

	if (!this->read_features()) {
		G_TASK_WARN("Features loading failed.");
		G_TASK_END;
	}

	if (!this->read_matrix()) {
		G_TASK_WARN("Matrix loading failed.");
		G_TASK_END;
	}
	determine_species();
	separate_counts();
	calculate_metadata();
	emit x_data_create_soon(this->single_cell_multiome_.release(), soap::VariableType::SingleCellMultiome, "SingleCellMultiome");
	G_TASK_END;
};

bool Read10XMultiomeWorker::read_barcodes() {

	auto file = this->barcodes_file_name_.toStdString();

	std::unique_ptr<char[]> bu(new char[256]);

	char* buffer = bu.get();

	gzFile barcodes = gzopen(file.data(), "rb");

	if (barcodes == NULL) {
		return false;
	}

	while ((gzgets(barcodes, buffer, 1024)) != 0)
	{
		this->barcodes_ << QString::fromUtf8(buffer, custom::line_length(buffer));;
	}
	gzclose(barcodes);

	return true;
};

bool Read10XMultiomeWorker::read_features() {

	auto file = this->features_file_name_.toStdString();

	std::unique_ptr<char[]> bu(new char[1024]);

	char* buffer = bu.get();

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
		char* name1 = c;
		while (*c != '\t') {
			++c;
		}
		char* name2 = c;
		++c;
		while (*c != '\t') {
			++c;
		}
		int type_length = c - name2 - 1;
		this->feature_names_ << QString::fromUtf8(name2 + 1, c - name2 - 1);
		if (type_length == 15) {
			this->gene_symbols_ << QString::fromUtf8(name1, name2 - name1);
		}
		else {
			this->peak_names_ << QString::fromUtf8(name1, name2 - name1);
		}
	}
	gzclose(features);

	return true;
};

void Read10XMultiomeWorker::determine_species() {
	for (QString& i : this->gene_symbols_) {
		if (i.startsWith("MT-")) {
			this->single_cell_multiome_->species_ = soap::Species::Human;
			break;
		}
		if (i.startsWith("mt-")) {
			this->single_cell_multiome_->species_ = soap::Species::Mouse;
			break;
		}
	}
};

bool Read10XMultiomeWorker::read_matrix() {

	auto file = this->matrix_file_name_.toStdString();

	std::unique_ptr<char[]> bu(new char[1024]);

	char* buffer = bu.get();

	gzFile matrix = gzopen(file.data(), "rb");

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
	int n_number = custom::atoi_specialized(&c);

	if (n_row != this->gene_symbols_.size() + this->peak_names_.size() || n_column != this->barcodes_.length()) {
		return false;
	}

	Eigen::SparseMatrix<int, Eigen::RowMajor> trmat(n_row, n_column);
	typename Eigen::SparseMatrix<int>::IndexVector wi(n_row);
	wi.setZero();

	int* data = new int[n_number * 3];
	for (int i = 0; i < n_number; ++i) {
		gzgets(matrix, buffer, 1024);
		c = buffer;
		data[3 * i] = custom::atoi_specialized(&c) - 1;
		data[3 * i + 1] = custom::atoi_specialized(&c) - 1;
		data[3 * i + 2] = custom::atoi_specialized(&c);
		++wi(data[3 * i]);
	}

	gzclose(matrix);

	trmat.reserve(wi);
	for (int i = 0; i < n_number; ++i) {
		trmat.insertBackUncompressed(data[3 * i], data[3 * i + 1]) = data[3 * i + 2];
	}
	delete[]data;
	trmat.collapseDuplicates(Eigen::internal::scalar_sum_op<int, int>());

	this->counts_ = trmat;

	return true;
};

void Read10XMultiomeWorker::separate_counts() {

	SparseInt& rna_counts = SUBMODULES(*this->single_cell_multiome_->rna_field(), SparseInt)[VARIABLE_RNA_COUNTS];
	rna_counts.data_type_ = SparseInt::DataType::Counts;

	SparseInt& atac_counts = SUBMODULES(*this->single_cell_multiome_->atac_field(), SparseInt)[VARIABLE_ATAC_COUNTS];
	atac_counts.data_type_ = SparseInt::DataType::Counts;

	QStringList barcodes = custom::make_unique(this->barcodes_);

	rna_counts.mat_ = custom::row_sliced(this->counts_, custom::equal(this->feature_names_, QString("Gene Expression")));
	rna_counts.rownames_ = custom::make_unique(this->gene_symbols_);
	rna_counts.colnames_ = barcodes;

	atac_counts.mat_ = custom::row_sliced(this->counts_, custom::equal(this->feature_names_, QString("Peaks")));
	atac_counts.rownames_ = custom::make_unique(this->peak_names_);
	atac_counts.colnames_ = barcodes;
}

void Read10XMultiomeWorker::calculate_metadata() {

	SparseInt& rna_counts = SUBMODULES(*this->single_cell_multiome_->rna_field(), SparseInt)[VARIABLE_RNA_COUNTS];
	Eigen::ArrayX<bool> gene_detected = custom::row_sum(rna_counts.mat_) > 0;
	rna_counts.row_slice(gene_detected);

	SparseInt& atac_counts = SUBMODULES(*this->single_cell_multiome_->atac_field(), SparseInt)[VARIABLE_ATAC_COUNTS];
	Eigen::ArrayX<bool> peak_detected = custom::row_sum(atac_counts.mat_) > 0;
	atac_counts.row_slice(peak_detected);

	Eigen::ArrayXi rna_count = custom::col_sum_mt(rna_counts.mat_);
	const int n_cell = rna_counts.mat_.cols();
	Eigen::ArrayXi rna_gene(n_cell);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(n_cell);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(n_cell);

	for (std::size_t i = 0; i < n_cell; ++i) {
		rna_gene[i] = rna_counts.mat_.outerIndexPtr()[i + 1] - rna_counts.mat_.outerIndexPtr()[i];
	}

	QVector<int> mitochondrial_location, ribosomal_location;

	auto& gene_symbols = rna_counts.rownames_;
	if (this->single_cell_multiome_->species_ == soap::Species::Human) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (this->single_cell_multiome_->species_ == soap::Species::Mouse) {
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
			mitochondrial_content += rna_counts.mat_.row(i).cast<double>();
		}
		mitochondrial_content /= rna_count.cast<double>();
	}

	if (ribosomal_location.length() > 0) {
		for (int& i : ribosomal_location) {
			ribosomal_content += rna_counts.mat_.row(i).cast<double>();
		}
		ribosomal_content /= rna_count.cast<double>();
	}

	custom::remove_na(mitochondrial_content);
	custom::remove_na(ribosomal_content);

	Metadata& metadata = SUBMODULES(*this->single_cell_multiome_, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(rna_counts.colnames_);
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, QVector<int>(rna_count.begin(), rna_count.end()));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, QVector<int>(rna_gene.begin(), rna_gene.end()));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content)); 
	metadata.mat_.update(METADATA_BARCODES, this->barcodes_);

	Eigen::ArrayXi peak_count = custom::col_sum_mt(atac_counts.mat_);
	const int ncol_atac = atac_counts.mat_.cols();
	Eigen::ArrayXi peak_gene(ncol_atac);

	for (std::size_t i = 0; i < ncol_atac; ++i) {
		peak_gene[i] = atac_counts.mat_.outerIndexPtr()[i + 1] - atac_counts.mat_.outerIndexPtr()[i];
	}

	metadata.mat_.update(METADATA_ATAC_UMI_NUMBER, custom::cast<QVector>(peak_count));
	metadata.mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, custom::cast<QVector>(peak_gene));
};