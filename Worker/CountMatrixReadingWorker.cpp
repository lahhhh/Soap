#include "CountMatrixReadingWorker.h"

#include "Custom.h"
#include <QtConcurrent>

CountMatrixReadingWorker::CountMatrixReadingWorker(const QString& file_path, const QString& delimiter) :
	file_path_(file_path),
	delimiter_(delimiter)
{}

bool CountMatrixReadingWorker::read_file() {

	QFile file(this->file_path_);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		G_TASK_WARN("Can not open " + this->file_path_);
		return false;
	}

	QTextStream in(&file);
	
	QString line = in.readLine();

	if (line.isNull()) {
		
		G_TASK_WARN("File is empty.");
		return false;
	}

	QString line2 = in.readLine();

	if (line2.isNull()) {
		
		G_TASK_WARN("File is broken.");
		return false;
	}

	QChar delimiter = this->delimiter_[0];

	if (this->delimiter_ == "auto-detect") {
		delimiter = custom::detect_delimiter(line, line2);
	}

	this->colnames_ = custom::digest(line, delimiter);
	this->ncol_ = this->colnames_.size();

	// common dataset should contains more than 100 cells
	if (this->ncol_ < 100) {
		G_TASK_WARN("Count Matrix Parse Failed.");
		
		return false;
	}

	if (this->colnames_[0].isEmpty()) {
		--this->ncol_;
		this->colnames_ = this->colnames_.sliced(1, this->ncol_);
	}

	QStringList first_row = custom::digest(line2, delimiter);
	if (first_row.size() == this->ncol_) {// means the first element should be "Gene Name" or something instead of empty ""
		--this->ncol_;
		this->colnames_ = this->colnames_.sliced(1, this->ncol_);
	}
	else if (first_row.size() != this->ncol_ + 1) {
		G_TASK_WARN("File is not valid.");
		
		return false;
	}

	this->row_name_map_[0] = first_row[0];
	for (int i = 1; i < this->ncol_ + 1; ++i) {
		if (first_row[i] != "0") {
			this->triplets_.emplace_back(0, i - 1, first_row[i].toInt());
		}
	}
	this->nrow_ = 1;

	bool success = true;
	std::mutex mutex, mutex2;
	int count = 0;
	QStringList cache;

	auto fun = [&cache, this, &mutex, &mutex2, &success, delimiter](int row_index) {

		QStringList row = custom::digest(cache[row_index], (delimiter));

		if (row.size() != this->ncol_ + 1) {
			success = false;
			return;
		}

		QVector<Eigen::Triplet<int>> sub_triplets;

		mutex.lock();
		this->row_name_map_[row_index + this->nrow_] = row[0];
		mutex.unlock();

		for (int i = 1; i < this->ncol_ + 1; ++i) {
			if (row[i] != "0") {
				sub_triplets.emplace_back(row_index + this->nrow_, i - 1, row[i].toInt());
			}
		}

		mutex2.lock();
		this->triplets_ << sub_triplets;
		mutex2.unlock();
	};

	while (!in.atEnd()) {
		cache.clear();
		count = 0;
		while (!in.atEnd() && count < 2000) {
			cache << in.readLine();
			++count;
		}
		QVector<int> params = custom::seq_n(0, count);
		QFuture<void> f = QtConcurrent::map(params, fun);
		f.waitForFinished();
		if (!success) {
			G_TASK_WARN("File is broken.");
			
			return false;
		}
		this->nrow_ += count;
	}

	return success;
};

void CountMatrixReadingWorker::create_data() {

	SingleCellRna* sc = new SingleCellRna();
	SparseInt& counts = SUBMODULES(*sc, SparseInt)[VARIABLE_COUNTS];
	counts.data_type_ = SparseInt::DataType::Counts;
	counts.mat_.resize(this->nrow_, this->ncol_);
	counts.mat_.reserve(this->triplets_.size());
	counts.mat_.setFromTriplets(this->triplets_.cbegin(), this->triplets_.cend());

	this->triplets_.clear();
	QStringList gene_symbols(this->nrow_);
	for (int i = 0; i < this->nrow_; ++i) {
		gene_symbols[i] = this->row_name_map_[i];
	}

	Eigen::ArrayX<bool> gene_detected = custom::row_sum(counts.mat_) > 0;

	// gene number should be more than 1000
	if (gene_detected.count() < 1000) {
		G_TASK_WARN("Count Matrix Loading Failed.");
		delete sc;
		return;
	}

	counts.mat_ = custom::row_sliced(counts.mat_, gene_detected);
	gene_symbols = custom::sliced(gene_symbols, gene_detected);

	Eigen::ArrayXi col_count = custom::col_sum(counts.mat_);
	Eigen::ArrayXi col_gene(this->ncol_);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(this->ncol_);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(this->ncol_);

	for (int i = 0; i < this->ncol_; ++i) {
		col_gene[i] = counts.mat_.outerIndexPtr()[i + 1] - counts.mat_.outerIndexPtr()[i];
	}
	QList<int> mitochondrial_location, ribosomal_location;
	soap::Species species = soap::Species::Undefined;
	for (QString& i : gene_symbols) {
		if (i.startsWith("MT-")) {
			species = soap::Species::Human;
			break;
		}
		if (i.startsWith("mt-")) {
			species = soap::Species::Mouse;
			break;
		}
	}
	if (species == soap::Species::Human) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (species == soap::Species::Mouse) {
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

	gene_symbols = custom::make_unique(gene_symbols);
	this->colnames_ = custom::make_unique(this->colnames_);
	counts.rownames_ = gene_symbols;
	counts.colnames_ = this->colnames_;
	Metadata& metadata = SUBMODULES(*sc, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(this->colnames_);
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, custom::cast<QVector>(col_count));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, custom::cast<QVector>(col_gene));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content));

	sc->species_ = species;

	emit x_data_create_soon(sc, soap::VariableType::SingleCellRna, "SingleCellRna");
}

void CountMatrixReadingWorker::run() {

	if (this->read_file()) {

		this->create_data();
	}

	G_TASK_END;
}
