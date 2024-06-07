#include "CountMatrixBatchLoadingWorker.h"

#include <QFile>
#include <QtConcurrent/QtConcurrent>

void CountMatrixBatchLoadingWorker::run() {

	int n_object = this->file_paths_.size();

	this->objects_.resize(n_object);

	for (int i = 0; i < n_object; ++i) {
		if (!this->load_single_object(this->objects_[i], this->file_paths_[i])) {
			G_TASK_WARN("Trouble in loading " + this->file_paths_[i]);
			G_TASK_END;
		}
	}

	this->integrate_objects();

	G_TASK_END;
};

QList<std::pair<QVector<int>, QVector<int> > >
CountMatrixBatchLoadingWorker::integrate_sparseint(
	SparseInt& to,
	QList<SparseInt const*> froms,
	bool distinguish_barcode)
{

	QList<std::pair<QVector<int>, QVector<int> > > index_map;

	QStringList barcodes, gene_names;
	QList<int> locations;

	for (auto ptr : froms) {
		locations << barcodes.size();
		barcodes << ptr->colnames_;
		gene_names << ptr->rownames_;
	}
	gene_names = _Cs unique(gene_names);

	int nrow = gene_names.size();

	if (distinguish_barcode) {

		barcodes = _Cs make_unique(barcodes);
		int ncol = barcodes.size(), index = 0;
		to.mat_.resize(nrow, ncol);
		std::vector<Eigen::Triplet<int> > triplets;

		for (auto from : froms) {
			int loc = locations[index++];
			QVector<int> row_map = _Cs index_of(from->rownames_, gene_names);
			for (int k = 0; k < from->mat_.cols(); ++k) {
				for (Eigen::SparseMatrix<int>::InnerIterator it(from->mat_, k); it; ++it) {
					triplets.emplace_back(row_map[it.row()], k + loc, it.value());
				}
			}

			index_map.emplace_back(std::make_pair(row_map, _Cs seq_n(loc, from->colnames_.size())));
		}

		to.mat_.setFromTriplets(triplets.cbegin(), triplets.cend());
		triplets.clear();
		to.rownames_ = gene_names;
		to.colnames_ = barcodes;

		return index_map;
	}
	else {

		barcodes = _Cs unique(barcodes);
		int ncol = barcodes.size(), index = 0;
		to.mat_.resize(nrow, ncol);
		std::vector<Eigen::Triplet<int> > triplets;

		for (auto from : froms) {
			QVector<int> row_map = _Cs index_of(from->rownames_, gene_names);
			QVector<int> col_map = _Cs index_of(from->colnames_, barcodes);
			for (int k = 0; k < from->mat_.outerSize(); ++k) {
				for (Eigen::SparseMatrix<int>::InnerIterator it(from->mat_, k); it; ++it) {
					triplets.emplace_back(row_map[it.row()], col_map[k], it.value());
				}
			}

			index_map.emplace_back(std::make_pair(row_map, col_map));
		}

		to.mat_.setFromTriplets(triplets.cbegin(), triplets.cend());
		triplets.clear();
		to.rownames_ = gene_names;
		to.colnames_ = barcodes;

		return index_map;
	}
}


void CountMatrixBatchLoadingWorker::integrate_objects() {

	auto sp = _Cs unique(_Cs sapply(this->objects_, [](auto&& o) {return o.species_; }));

	SingleCellRna* single_cell_rna = new SingleCellRna();
	single_cell_rna->species_ = sp[0];
	single_cell_rna->data_type_ = SingleCellRna::DataType::Integrated;

	SparseInt& counts = SUBMODULES(*single_cell_rna, SparseInt)[VARIABLE_COUNTS];
	counts.data_type_ = SparseInt::DataType::Counts;

	QList<SparseInt const*> list_of_counts = _Cs sapply(this->objects_,
		[](auto&& data) {return static_cast<SparseInt const*>(data.counts()); });

	this->integrate_sparseint(counts, list_of_counts, true);

	Metadata& metadata = SUBMODULES(*single_cell_rna, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);

	QList<Metadata*> list_of_metadata = _Cs sapply(this->objects_, [](auto&& data) {return data.metadata(); });

	this->integrate_metadata(metadata, list_of_metadata);

	emit x_data_create_soon(single_cell_rna, soap::VariableType::SingleCellRna, "SingleCellRna");
};

void CountMatrixBatchLoadingWorker::integrate_metadata(Metadata& to, QList<Metadata*> froms) {
	QStringList metadata_names;
	QMap<QString, CustomMatrix::DataType > final_data_type;
	for (auto data : froms) {
		metadata_names << data->mat_.colnames_;
	}
	metadata_names = _Cs unique(metadata_names);

	for (const auto& name : metadata_names) {
		CustomMatrix::DataType data_type = CustomMatrix::DataType::NoType;
		bool accepted = true;
		for (auto data : froms) {
			if (data->mat_.contains(name)) {
				CustomMatrix::DataType type = data->mat_.data_type_[name];
				if (data_type != type) {
					if (data_type == CustomMatrix::DataType::NoType) {
						data_type = type;
					}
					else {
						accepted = false;
					}
				}
			}
		}
		if (accepted) {
			final_data_type[name] = data_type;
		}
	}
	for (const auto& name : final_data_type.keys()) {
		CustomMatrix::DataType data_type = final_data_type[name];
		if (data_type == CustomMatrix::DataType::IntegerNumeric) {
			QVector<int> final_data;
			for (auto data : froms) {
				if (data->mat_.contains(name)) {
					final_data << data->mat_.get_const_integer_reference(name);
				}
				else {
					final_data << QVector<int>(data->mat_.rows(), -1);
				}
			}
			to.mat_.update(name, final_data, CustomMatrix::DataType::IntegerNumeric);
		}
		else if (data_type == CustomMatrix::DataType::IntegerFactor) {
			QVector<int> final_data;
			for (auto data : froms) {
				if (data->mat_.contains(name)) {
					final_data << data->mat_.get_const_integer_reference(name);
				}
				else {
					final_data << QVector<int>(data->mat_.rows(), -1);
				}
			}
			to.mat_.update(name, final_data, CustomMatrix::DataType::IntegerFactor);
		}
		else if (data_type == CustomMatrix::DataType::DoubleNumeric) {
			QVector<double> final_data;
			for (auto data : froms) {
				if (data->mat_.contains(name)) {
					final_data << data->mat_.get_const_double_reference(name);
				}
				else {
					final_data << QVector<double>(data->mat_.rows(), -1);
				}
			}
			to.mat_.update(name, final_data, CustomMatrix::DataType::DoubleNumeric);
		}
		else if (data_type == CustomMatrix::DataType::QStringFactor) {
			QStringList final_data;
			for (auto data : froms) {
				if (data->mat_.contains(name)) {
					final_data << data->mat_.get_const_qstring_reference(name);
				}
				else {
					final_data << QStringList(data->mat_.rows(), "NA");
				}
			}
			to.mat_.update(name, final_data, CustomMatrix::DataType::QStringFactor);
		}
		else if (data_type == CustomMatrix::DataType::QString) {
			QStringList final_data;
			for (auto data : froms) {
				if (data->mat_.contains(name)) {
					final_data << data->mat_.get_const_qstring_reference(name);
				}
				else {
					final_data << QStringList(data->mat_.rows(), "NA");
				}
			}
			to.mat_.update(name, final_data, CustomMatrix::DataType::QString);
		}
	}
}

bool CountMatrixBatchLoadingWorker::load_single_object(SingleCellRna& object, const QString& file_path) {

	QFile file(file_path);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		return false;
	}

	QTextStream in(&file);

	QString line = in.readLine();

	if (line.isNull()) {
		
		return false;
	}

	QString line2 = in.readLine();

	if (line2.isNull()) {
		
		return false;
	}

	QChar delimiter = this->delimiter_[0];

	if (this->delimiter_ == "auto-detect") {
		delimiter = _Cs detect_delimiter(line, line2);
	}

	QStringList colnames = _Cs digest(line, delimiter);
	int ncol = colnames.size();

	// common dataset should contains more than 100 cells
	if (ncol < 100) {
		
		return false;
	}

	if (colnames[0].isEmpty()) {
		--ncol;
		colnames = colnames.sliced(1, ncol);
	}

	QStringList first_row = _Cs digest(line2, delimiter);
	if (first_row.size() == ncol) {
		--ncol;
		colnames = colnames.sliced(1, ncol);
	}
	else if (first_row.size() != ncol + 1) {
		return false;
	}

	QVector<Eigen::Triplet<int>> triplets;

	QMap<int, QString> row_name_map;

	row_name_map[0] = first_row[0];
	for (int i = 1; i < ncol + 1; ++i) {
		if (first_row[i] != "0") {
			triplets.emplace_back(0, i - 1, first_row[i].toInt());
		}
	}

	bool success = true;
	std::mutex mutex, mutex2;
	int count{ 0 };
	int nrow{ 0 };
	QStringList cache;

	auto fun = [&cache, &mutex, &mutex2, &success, delimiter, ncol, &row_name_map, &nrow, &triplets](int row_index) {

		QStringList row = _Cs digest(cache[row_index], delimiter);

		if (row.size() != ncol + 1) {
			success = false;
			return;
		}

		QVector<Eigen::Triplet<int>> sub_triplets;

		mutex.lock();
		row_name_map[row_index + nrow] = row[0];
		mutex.unlock();

		for (int i = 1; i < ncol + 1; ++i) {
			if (row[i] != "0") {
				sub_triplets.emplace_back(row_index + nrow, i - 1, row[i].toInt());
			}
		}

		mutex2.lock();
		triplets << sub_triplets;
		mutex2.unlock();
	};

	while (!in.atEnd()) {
		cache.clear();
		count = 0;
		while (!in.atEnd() && count < 2000) {
			cache << in.readLine();
			++count;
		}

		QVector<int> params = _Cs seq_n(0, count);
		QFuture<void> f = QtConcurrent::map(params, fun);
		f.waitForFinished();

		if (!success) {
			
			return false;
		}
		nrow += count;
	}


	SparseInt& counts = SUBMODULES(object, SparseInt)[VARIABLE_COUNTS];
	counts.data_type_ = SparseInt::DataType::Counts;
	counts.mat_.resize(nrow, ncol);
	counts.mat_.reserve(triplets.size());
	counts.mat_.setFromTriplets(triplets.cbegin(), triplets.cend());

	triplets.clear();
	QStringList gene_symbols(nrow);
	for (int i = 0; i < nrow; ++i) {
		gene_symbols[i] = row_name_map[i];
	}

	Eigen::ArrayX<bool> gene_detected = _Cs row_sum(counts.mat_) > 0;

	// gene number should be more than 1000
	if (gene_detected.count() < 1000) {
		return false;
	}

	counts.mat_ = _Cs row_sliced(counts.mat_, gene_detected);
	gene_symbols = _Cs sliced(gene_symbols, gene_detected);

	Eigen::ArrayXi col_count = _Cs col_sum(counts.mat_);
	Eigen::ArrayXi col_gene(ncol);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(ncol);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(ncol);

	for (int i = 0; i < ncol; ++i) {
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

	_Cs remove_na(mitochondrial_content);
	_Cs remove_na(ribosomal_content);

	gene_symbols = _Cs make_unique(gene_symbols);
	colnames = _Cs make_unique(colnames);
	counts.rownames_ = gene_symbols;
	counts.colnames_ = colnames;

	Metadata& metadata = SUBMODULES(object, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(colnames);
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, _Cs cast<QVector>(col_count));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, _Cs cast<QVector>(col_gene));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, _Cs cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, _Cs cast<QVector>(ribosomal_content));
	metadata.mat_.update("Source", QStringList(ncol, file_path), CustomMatrix::DataType::QStringFactor);

	object.species_ = species;

	return true;
};