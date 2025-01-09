#include "SingleCellRnaCreateWorker.h"

void SingleCellRnaCreateWorker::run() {

	G_TASK_LOG("Start creating object...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_single_cell_rna_created(this->res_.release());

	G_TASK_LOG("Object created.");

	G_TASK_END;
};

bool SingleCellRnaCreateWorker::work() {

	switch (this->mode_)
	{
	case WorkMode::FromSparseInt:
		return this->create_from_sparseint();
		break;
	default:
		break;
	}

	return false;
};

bool SingleCellRnaCreateWorker::create_from_sparseint() {

	if (this->si_->rows() < 100 || this->si_->cols() < 100 || this->si_->mat_.nonZeros() < 100) {
		G_TASK_WARN("Too Small Object for SingleCellRna");
		return false;
	}

	this->res_.reset(new SingleCellRna());

	SparseInt& counts = SUBMODULES(*this->res_, SparseInt)[VARIABLE_COUNTS];

	counts = *this->si_;
	counts.data_type_ = SparseInt::DataType::Counts;

	for (auto&& i : counts.rownames_) {
		if (i.startsWith("MT-")) {
			this->res_->species_ = soap::Species::Human;
			break;
		}
		if (i.startsWith("mt-")) {
			this->res_->species_ = soap::Species::Mouse;
			break;
		}
	}

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
	
	return true;
};