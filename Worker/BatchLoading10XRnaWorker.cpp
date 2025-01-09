#include "BatchLoading10XRnaWorker.h"

#include "Read10xRnaWorker.h"

#include <QDir>
#include <QDirIterator>
#include <QFileInfo>

bool BatchLoading10XRnaWorker::work() {

	QDir dir(this->dir_);
	if (!dir.exists()) {
		G_TASK_WARN("Invalid Folder");
		return false;
	}

	this->objects_.reserve(1000);

	QDirIterator dir_it(this->dir_, QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::Subdirectories);
	while (dir_it.hasNext()) {
		dir_it.next();
		QDir current_dir(dir_it.filePath());

		QString barcodes_file;
		QString features_file;
		QString matrix_file;

		QDirIterator file_iterator(current_dir.absolutePath(), QDir::Files);
		while (file_iterator.hasNext()) {
			file_iterator.next();
			QFileInfo file_info(file_iterator.fileInfo());
			QString file_name = file_info.absoluteFilePath();

			if (file_name.endsWith("barcodes.tsv.gz")) {
				barcodes_file = file_name;
			}
			else if (file_name.endsWith("features.tsv.gz")) {
				features_file = file_name;
			}
			else if (file_name.endsWith("matrix.mtx.gz")) {
				matrix_file = file_name;
			}
		}

		if (!barcodes_file.isEmpty() && !features_file.isEmpty() && !matrix_file.isEmpty()) {
			this->objects_ << SingleCellRna();

			if (!this->load_single_object(this->objects_.last(), file_iterator.path(), barcodes_file, features_file, matrix_file)) {
				G_TASK_WARN("Trouble in loading file from " + file_iterator.path());
				return false;
			}
		}
	}

	this->integrate_objects();
}

void BatchLoading10XRnaWorker::run() {

	G_TASK_LOG("Start loading...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_data_create_soon(this->res_.release(), soap::VariableType::SingleCellRna, "SingleCellRna");

	G_TASK_LOG("Loading finished.");

	G_TASK_END;
};


QList<std::pair<QVector<int>, QVector<int> > >
BatchLoading10XRnaWorker::integrate_sparseint(
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
	gene_names = custom::unique(gene_names);

	int nrow = gene_names.size();

	if (distinguish_barcode) {

		barcodes = custom::make_unique(barcodes);
		int ncol = barcodes.size(), index = 0;
		to.mat_.resize(nrow, ncol);
		std::vector<Eigen::Triplet<int> > triplets;

		for (auto from : froms) {
			int loc = locations[index++];
			QVector<int> row_map = custom::index_of(from->rownames_, gene_names);
			for (int k = 0; k < from->mat_.cols(); ++k) {
				for (Eigen::SparseMatrix<int>::InnerIterator it(from->mat_, k); it; ++it) {
					triplets.emplace_back(row_map[it.row()], k + loc, it.value());
				}
			}

			index_map.emplace_back(std::make_pair(row_map, custom::seq_n(loc, from->colnames_.size())));
		}

		to.mat_.setFromTriplets(triplets.cbegin(), triplets.cend());
		triplets.clear();
		to.rownames_ = gene_names;
		to.colnames_ = barcodes;

		return index_map;
	}
	else {

		barcodes = custom::unique(barcodes);
		int ncol = barcodes.size(), index = 0;
		to.mat_.resize(nrow, ncol);
		std::vector<Eigen::Triplet<int> > triplets;

		for (auto from : froms) {
			QVector<int> row_map = custom::index_of(from->rownames_, gene_names);
			QVector<int> col_map = custom::index_of(from->colnames_, barcodes);
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


void BatchLoading10XRnaWorker::integrate_objects() {

	auto sp = custom::unique(custom::sapply(this->objects_, [](auto&& o) {return o.species_; }));

	this->res_.reset(new SingleCellRna());

	this->res_->species_ = sp[0];
	this->res_->data_type_ = SingleCellRna::DataType::Integrated;

	SparseInt& counts = SUBMODULES(*this->res_, SparseInt)[VARIABLE_COUNTS];
	counts.data_type_ = SparseInt::DataType::Counts;

	QList<SparseInt const*> list_of_counts = custom::sapply(this->objects_,
		[](auto&& data) {return static_cast<SparseInt const*>(data.counts()); });

	this->integrate_sparseint(counts, list_of_counts, true);

	Metadata& metadata = SUBMODULES(*this->res_, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);

	QList<Metadata*> list_of_metadata = custom::sapply(this->objects_, [](auto&& data) {return data.metadata(); });

	this->integrate_metadata(metadata, list_of_metadata);
};

void BatchLoading10XRnaWorker::integrate_metadata(Metadata& to, QList<Metadata*> froms) {
	QStringList metadata_names;
	QMap<QString, CustomMatrix::DataType > final_data_type;
	for (auto data : froms) {
		metadata_names << data->mat_.colnames_;
	}
	metadata_names = custom::unique(metadata_names);

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

bool BatchLoading10XRnaWorker::load_single_object(
	SingleCellRna& object,
	const QString& path,
	const QString& barcodes_file_name,
	const QString& features_file_name,
	const QString& matrix_file_name) {

	Read10XRnaWorker worker(barcodes_file_name, features_file_name, matrix_file_name);

	if (!worker.work()) {
		return false;
	}

	SingleCellRna* ptr = worker.res_.release();

	int n_cell = ptr->counts()->mat_.cols();

	ptr->metadata()->mat_.update("Source", QStringList(n_cell, path), CustomMatrix::DataType::QStringFactor);

	object = std::move(*ptr);
	delete ptr;

	return true;
};