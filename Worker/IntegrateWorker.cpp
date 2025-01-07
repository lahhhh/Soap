#include "IntegrateWorker.h"

#include "MacsCallPeakWorker.h"
#include "CalculateCountsByGenomicRangeWorker.h"

std::pair<QStringList, QList<QVector<int>>>
IntegrateWorker::get_index(const QList<QStringList>& components, bool distinguish) {

	QStringList concatenated = custom::unroll(components);

	QList<QVector<int>> index;

	if (!distinguish) {

		concatenated = custom::unique(concatenated);

		for (auto&& comp : components) {
			index << custom::index_of(comp, concatenated);
		}
	}
	else {

		concatenated = custom::make_unique(concatenated);

		int loc{ 0 };
		for (auto&& comp : components) {
			int size = comp.size();
			index << custom::seq_n(loc, size);
			loc += size;
		}
	}

	return { concatenated, index };
};

void IntegrateWorker::map_dense_int_value(
	DenseInt& to,
	QList<DenseInt const*> froms,
	const QList<QVector<int>>& row_maps,
	const QList<QVector<int>>& col_maps,
	const QStringList& rownames,
	const QStringList& colnames) {

	to.rownames_ = rownames;
	to.colnames_ = colnames;
	int nrow = to.rownames_.size();
	int ncol = to.colnames_.size();
	to.mat_ = Eigen::MatrixXi::Zero(nrow, ncol);

	int n_data = froms.size();
	auto&& m = to.mat_;
	for (int i = 0; i < n_data; ++i) {
		auto&& mat = froms[i]->mat_;

		int n_sub_col = mat.cols();
		int n_sub_row = mat.rows();

		auto&& row_map = row_maps[i];
		auto&& col_map = col_maps[i];

		for (int k = 0; k < n_sub_col; ++k) {
			for (int j = 0; j < n_sub_row; ++j) {
				m(row_map[j], col_map[k]) += mat(j, k);
			}
		}
	}
};

void IntegrateWorker::map_sparse_int_value(
	SparseInt& to,
	QList<SparseInt const*> froms,
	const QList<QVector<int>>& row_maps,
	const QList<QVector<int>>& col_maps,
	const QStringList& rownames,
	const QStringList& colnames)
{
	to.rownames_ = rownames;
	to.colnames_ = colnames;
	int nrow = to.rownames_.size();
	int ncol = to.colnames_.size();
	to.mat_.resize(nrow, ncol);

	std::vector<Eigen::Triplet<int> > triplets;

	int n_data = froms.size();

	for (int i = 0; i < n_data; ++i) {
		auto&& mat = froms[i]->mat_;
		int n_sub_col = mat.cols();
		auto&& row_map = row_maps[i];
		auto&& col_map = col_maps[i];

		for (int k = 0; k < n_sub_col; ++k) {
			for (Eigen::SparseMatrix<int>::InnerIterator it(mat, k); it; ++it) {
				triplets.emplace_back(row_map[it.row()], col_map[k], it.value());
			}
		}
	}

	to.mat_.setFromTriplets(triplets.cbegin(), triplets.cend());
	triplets.clear();
}

void IntegrateWorker::integrate_velocyto(VelocytoBase& to, QList<const VelocytoBase*> froms, bool distinguish_barcode) {

	auto& spliced = SUBMODULES(to, SparseInt)[VARIABLE_RNA_SPLICED];
	spliced.data_type_ = SparseInt::DataType::Spliced;
	auto spliced_froms = custom::sapply(froms, [](auto* f) {return f->get_spliced(); });

	auto [counts_rownames, row_map] = this->get_index(
		custom::sapply(spliced_froms, [](auto&& t) {return t->rownames_; }),
		false);

	auto [counts_colnames, col_map] = this->get_index(
		custom::sapply(spliced_froms, [](auto&& t) {return t->colnames_; }),
		this->distinguish_barcode_);

	this->map_sparse_int_value(spliced, spliced_froms, row_map, col_map, counts_rownames, counts_colnames);


	auto& unspliced = SUBMODULES(to, SparseInt)[VARIABLE_RNA_UNSPLICED];
	unspliced.data_type_ = SparseInt::DataType::Unspliced;
	auto unspliced_froms = custom::sapply(froms, [](auto* f) {return f->get_unspliced(); });

	auto [counts_rownames2, row_map2] = this->get_index(
		custom::sapply(unspliced_froms, [](auto&& t) {return t->rownames_; }),
		false);

	auto [counts_colnames2, col_map2] = this->get_index(
		custom::sapply(unspliced_froms, [](auto&& t) {return t->colnames_; }),
		this->distinguish_barcode_);

	this->map_sparse_int_value(unspliced, unspliced_froms, row_map2, col_map2, counts_rownames2, counts_colnames2);

};

void IntegrateWorker::integrate_metadata(Metadata& to, QList<const Metadata*> froms) {
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
				CustomMatrix::DataType type = data->mat_.data_type_.at(name);
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

void IntegrateWorker::multiome_quality_control(SingleCellMultiome& single_cell_multiome) {
	SparseInt& rna_counts = SUBMODULES(*single_cell_multiome.rna_field(), SparseInt)[VARIABLE_RNA_COUNTS];
	Eigen::ArrayX<bool> gene_detected = custom::row_sum(rna_counts.mat_) > 0;
	rna_counts.row_slice(gene_detected);

	SparseInt& atac_counts = SUBMODULES(*single_cell_multiome.atac_field(), SparseInt)[VARIABLE_ATAC_COUNTS];
	Eigen::ArrayX<bool> peak_detected = custom::row_sum(atac_counts.mat_) > 0;
	atac_counts.row_slice(peak_detected);

	Eigen::ArrayXi rna_count = custom::col_sum_mt(rna_counts.mat_);
	const int ncol_rna = rna_counts.mat_.cols();
	Eigen::ArrayXi rna_gene(ncol_rna);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(ncol_rna);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(ncol_rna);

	for (std::size_t i = 0; i < ncol_rna; ++i) {
		rna_gene[i] = rna_counts.mat_.outerIndexPtr()[i + 1] - rna_counts.mat_.outerIndexPtr()[i];
	}

	QVector<int> mitochondrial_location, ribosomal_location;

	auto& gene_symbols = rna_counts.rownames_;
	if (single_cell_multiome.species_ == soap::Species::Human) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (single_cell_multiome.species_ == soap::Species::Mouse) {
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

	Metadata& metadata = SUBMODULES(single_cell_multiome, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(rna_counts.colnames_);
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, QVector<int>(rna_count.begin(), rna_count.end()));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, QVector<int>(rna_gene.begin(), rna_gene.end()));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content));

	Eigen::ArrayXi peak_count = custom::col_sum_mt(atac_counts.mat_);
	const int ncol_atac = atac_counts.mat_.cols();
	Eigen::ArrayXi peak_gene(ncol_atac);

	for (std::size_t i = 0; i < ncol_atac; ++i) {
		peak_gene[i] = atac_counts.mat_.outerIndexPtr()[i + 1] - atac_counts.mat_.outerIndexPtr()[i];
	}

	metadata.mat_.update(METADATA_ATAC_UMI_NUMBER, QVector<int>(peak_count.begin(), peak_count.end()));
	metadata.mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, QVector<int>(peak_gene.begin(), peak_gene.end()));
};

void IntegrateWorker::atac_quality_control(SingleCellAtac& single_cell_atac) {

	SparseInt& counts = *single_cell_atac.counts();

	Eigen::ArrayXi col_count = custom::col_sum_mt(counts.mat_);
	const int ncol = counts.mat_.cols();
	Eigen::ArrayXi col_peak(ncol);

	for (int i = 0; i < ncol; ++i) {
		col_peak[i] = counts.mat_.outerIndexPtr()[i + 1] - counts.mat_.outerIndexPtr()[i];
	}

	Metadata& metadata = SUBMODULES(single_cell_atac, Metadata)[VARIABLE_METADATA];
	metadata.mat_.update(METADATA_ATAC_UMI_NUMBER, custom::cast<QVector>(col_count));
	metadata.mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, custom::cast<QVector>(col_peak));
	metadata.mat_.update(METADATA_BARCODES, single_cell_atac.fragments()->cell_names_);
};

void IntegrateWorker::rna_quality_control(SingleCellRna& single_cell_rna) {
	Eigen::SparseMatrix<int>& counts = SUBMODULES(single_cell_rna, SparseInt)[VARIABLE_COUNTS].mat_;
	Eigen::ArrayXi col_count = custom::col_sum(counts);
	const int ncol = counts.cols();
	Eigen::ArrayXi col_gene(ncol);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(ncol);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(ncol);

	for (int i = 0; i < ncol; ++i) {
		col_gene[i] = counts.outerIndexPtr()[i + 1] - counts.outerIndexPtr()[i];
	}
	QList<int> mitochondrial_location, ribosomal_location;
	QStringList& gene_names = SUBMODULES(single_cell_rna, SparseInt)[VARIABLE_COUNTS].rownames_;
	if (single_cell_rna.species_ == soap::Species::Human) {
		for (int i = 0; i < gene_names.length(); ++i) {
			if (gene_names.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_names.at(i).startsWith("RPS") || gene_names.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (single_cell_rna.species_ == soap::Species::Mouse) {
		for (int i = 0; i < gene_names.length(); ++i) {
			if (gene_names.at(i).startsWith("mt-")) {
				mitochondrial_location << i;
			}

			if (gene_names.at(i).startsWith("Rps") || gene_names.at(i).startsWith("Rpl")) {
				ribosomal_location << i;
			}
		}
	}

	if (mitochondrial_location.length() > 0) {
		for (int& i : mitochondrial_location) {
			mitochondrial_content += counts.row(i).cast<double>();
		}
		mitochondrial_content /= col_count.cast<double>();
	}

	if (ribosomal_location.length() > 0) {
		for (int& i : ribosomal_location) {
			ribosomal_content += counts.row(i).cast<double>();
		}
		ribosomal_content /= col_count.cast<double>();
	}

	custom::remove_na(mitochondrial_content);
	custom::remove_na(ribosomal_content);

	CustomMatrix& metadata = SUBMODULES(single_cell_rna, Metadata)[VARIABLE_METADATA].mat_;

	metadata.update(METADATA_RNA_UMI_NUMBER, custom::cast<QVector>(col_count));
	metadata.update(METADATA_RNA_UNIQUE_GENE_NUMBER, custom::cast<QVector>(col_gene));
	metadata.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content));
}

void IntegrateWorker::bulkrna_mode() {
	BulkRna* bulk_rna = new BulkRna();
	bulk_rna->species_ = this->species_;
	bulk_rna->data_type_ = BulkRna::DataType::Integrated;

	DenseInt& counts = SUBMODULES(*bulk_rna, DenseInt)[VARIABLE_COUNTS];
	counts.data_type_ = DenseInt::DataType::Counts;

	QList<DenseInt const*> list_of_counts = custom::sapply(this->bulkrna_data_,
		[](auto* data) {return static_cast<DenseInt const*>(data->counts()); });

	auto [counts_rownames, row_map] = this->get_index(
		custom::sapply(list_of_counts, [](auto&& t) {return t->rownames_; }),
		false);

	auto [counts_colnames, col_map] = this->get_index(
		custom::sapply(list_of_counts, [](auto&& t) {return t->colnames_; }),
		this->distinguish_barcode_);

	this->map_dense_int_value(counts, list_of_counts, row_map, col_map, counts_rownames, counts_colnames);

	Metadata& metadata = SUBMODULES(*bulk_rna, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);

	QList<const Metadata*> list_of_metadata = custom::sapply(this->bulkrna_data_, [](auto* data) {return data->metadata(); });

	this->integrate_metadata(metadata, list_of_metadata);

	emit x_bulkrna_ready(bulk_rna, this->bulkrna_data_);
};

void IntegrateWorker::scrna_mode() {

	SingleCellRna* single_cell_rna = new SingleCellRna();
	single_cell_rna->species_ = this->species_;
	single_cell_rna->data_type_ = SingleCellRna::DataType::Integrated;

	SparseInt& counts = SUBMODULES(*single_cell_rna, SparseInt)[VARIABLE_COUNTS];
	counts.data_type_ = SparseInt::DataType::Counts;

	QList<SparseInt const*> list_of_counts = custom::sapply(this->scrna_data_,
		[](auto* data) {return static_cast<SparseInt const*>(data->counts()); });

	auto [counts_rownames, row_map] = this->get_index(
		custom::sapply(list_of_counts, [](auto&& t) {return t->rownames_; }),
		false);

	auto [counts_colnames, col_map] = this->get_index(
		custom::sapply(list_of_counts, [](auto&& t) {return t->colnames_; }),
		this->distinguish_barcode_);

	this->map_sparse_int_value(counts, list_of_counts, row_map, col_map, counts_rownames, counts_colnames);

	QList<const VelocytoBase*> list_of_velocyto = custom::sapply(this->scrna_data_,
		[](auto* data) {return data->velocyto_base(); });

	if (!list_of_velocyto.contains(nullptr)) {

		VelocytoBase& vb = SUBMODULES(*single_cell_rna, VelocytoBase)[VARIABLE_VELOCYTO_BASE];

		this->integrate_velocyto(vb, list_of_velocyto, this->distinguish_barcode_);
	}

	Metadata& metadata = SUBMODULES(*single_cell_rna, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(counts.colnames_);

	if (this->distinguish_barcode_) {

		QList<const Metadata*> list_of_metadata = custom::sapply(this->scrna_data_, [](auto* data) {return data->metadata(); });

		this->integrate_metadata(metadata, list_of_metadata);
	}

	this->rna_quality_control(*single_cell_rna);
	emit x_scrna_ready(single_cell_rna, this->scrna_data_);
};

bool IntegrateWorker::integrate_fragments_object(
	Fragments& to,
	const QList< QVector<int> >& col_map,
	const QList<const Fragments* >& fragmentss
)
{
	const int n_cell_all = to.cell_names_.size();

	for (int i = 0; i < fragmentss.size(); ++i) {
		const auto* fragments = fragmentss[i];
		const int n_cell = fragments->cell_names_.size();
		auto&& map = col_map[i];

		for (const auto& [sub_sequence_name, sub_sequence_data] : fragments->data_) {

			auto& all_sequence_data = to.data_[sub_sequence_name];

			if (all_sequence_data.empty()) {
				all_sequence_data.resize(n_cell_all);
			}

			for (int j = 0; j < n_cell; ++j) {
				auto& [start_all, end_all] = all_sequence_data[map[j]];
				auto& [start, end] = sub_sequence_data[j];

				start_all.insert(start_all.end(), start.begin(), start.end());
				end_all.insert(end_all.end(), end.begin(), end.end());
			}
		}
	}
	return true;
};

void IntegrateWorker::scatac_mode() {

	std::unique_ptr<SingleCellAtac> single_cell_atac(new SingleCellAtac());
	single_cell_atac->species_ = this->species_;
	single_cell_atac->data_type_ = SingleCellAtac::DataType::Integrated;

	// Fragments Integrate
	Fragments& fragments = SUBMODULES(*single_cell_atac, Fragments)[VARIABLE_FRAGMENTS];

	auto [cell_names, col_map] = this->get_index(
		custom::sapply(this->scatac_data_, [](auto&& t) {return t->counts()->colnames_; }),
		this->distinguish_barcode_
	);

	fragments.cell_names_ = cell_names;
	fragments.adjust_length_by_cell_name_length();

	if (!this->integrate_fragments_object(
		fragments,
		col_map,
		custom::sapply(this->scatac_data_, [](auto&& d) {return d->fragments(); })))
	{
		G_TASK_WARN("Fragments integration failed due to incomplete fragments files.");
		return;
	}

	MacsCallPeakWorker worker1({ &fragments });

	worker1.work();

	GenomicRange new_peak = worker1.res_;

	// ATAC Counts Integrate

	CalculateCountsByGenomicRangeWorker worker2(&fragments, new_peak);

	worker2.work();

	SUBMODULES(*single_cell_atac, SparseInt)[VARIABLE_COUNTS] = std::move(*worker2.res_);
	SUBMODULES(*single_cell_atac, SparseInt)[VARIABLE_COUNTS].data_type_ = SparseInt::DataType::Counts;

	worker2.res_.reset();
	// metadata

	Metadata& metadata = SUBMODULES(*single_cell_atac, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(cell_names);
	if (this->distinguish_barcode_) {

		QList<const Metadata*> list_of_metadata = custom::sapply(this->scatac_data_, [](auto* data) {return data->metadata(); });

		this->integrate_metadata(metadata, list_of_metadata);
	}

	this->atac_quality_control(*single_cell_atac);

	emit x_scatac_ready(single_cell_atac.release(), this->scatac_data_);
}

void IntegrateWorker::scmultiome_mode() {

	std::unique_ptr<SingleCellMultiome> single_cell_multiome(new SingleCellMultiome());
	single_cell_multiome->species_ = this->species_;
	single_cell_multiome->data_type_ = SingleCellMultiome::DataType::Integrated;

	auto& rna_field = single_cell_multiome->create_field("RNA", DataField::DataType::Rna);
	auto& atac_field = single_cell_multiome->create_field("ATAC", DataField::DataType::Atac);

	// RNA Counts Integrate

	SparseInt& rna_counts = SUBMODULES(rna_field, SparseInt)[VARIABLE_RNA_COUNTS];
	rna_counts.data_type_ = SparseInt::DataType::Counts;

	QList<SparseInt const*> list_of_counts = custom::sapply(this->scmultiome_data_,
		[](const SingleCellMultiome* data) {return static_cast<SparseInt const*>(data->rna_counts()); });

	auto [counts_rownames, row_map] = this->get_index(
		custom::sapply(list_of_counts, [](auto&& t) {return t->rownames_; }),
		false);

	auto [counts_colnames, col_map] = this->get_index(
		custom::sapply(list_of_counts, [](auto&& t) {return t->colnames_; }),
		this->distinguish_barcode_);

	this->map_sparse_int_value(rna_counts, list_of_counts, row_map, col_map, counts_rownames, counts_colnames);

	QList<const VelocytoBase*> list_of_velocyto = custom::sapply(this->scmultiome_data_,
		[](auto* data) {return data->velocyto_base(); });

	if (!list_of_velocyto.contains(nullptr)) {

		VelocytoBase& vb = SUBMODULES(*single_cell_multiome, VelocytoBase)[VARIABLE_VELOCYTO_BASE];

		this->integrate_velocyto(vb, list_of_velocyto, this->distinguish_barcode_);
	}

	// Fragments Integrate
	Fragments& fragments = SUBMODULES(*single_cell_multiome, Fragments)[VARIABLE_FRAGMENTS];
	fragments.cell_names_ = rna_counts.colnames_;
	fragments.adjust_length_by_cell_name_length();

	if (!this->integrate_fragments_object(
		fragments,
		col_map,
		custom::sapply(this->scmultiome_data_, [](auto&& d) {return d->fragments(); })))
	{
		G_TASK_WARN("Fragments integration failed due to incomplete fragments files.");
		return;
	}

	MacsCallPeakWorker worker1({ &fragments });

	worker1.work();

	GenomicRange new_peak = worker1.res_;

	// ATAC Counts Integrate

	CalculateCountsByGenomicRangeWorker worker2(&fragments, new_peak);

	worker2.work();

	SUBMODULES(atac_field, SparseInt)[VARIABLE_ATAC_COUNTS] = std::move(*worker2.res_);
	SUBMODULES(atac_field, SparseInt)[VARIABLE_ATAC_COUNTS].data_type_ = SparseInt::DataType::Counts;
	worker2.res_.reset();

	// Metadata Integrate

	Metadata& metadata = SUBMODULES(*single_cell_multiome, Metadata)[VARIABLE_METADATA];
	metadata.mat_.set_rownames(rna_counts.colnames_);
	if (this->distinguish_barcode_) {

		QList<const Metadata*> list_of_metadata = custom::sapply(this->scmultiome_data_, [](auto* data) {return data->metadata(); });

		this->integrate_metadata(metadata, list_of_metadata);
	}

	this->multiome_quality_control(*single_cell_multiome);

	emit x_scmultiome_ready(single_cell_multiome.release(), this->scmultiome_data_);
};


bool IntegrateWorker::work() {

	switch (this->mode_)
	{
	case WorkMode::BulkRna:
		bulkrna_mode();
		break;
	case WorkMode::SingleCellRna:
		scrna_mode();
		break;
	case WorkMode::SingleCellAtac:
		scatac_mode();
		break;
	case WorkMode::SingleCellMultiome:
		scmultiome_mode();
		break;
	default:
		break;
	}

	return true;
}

void IntegrateWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	G_TASK_END;
}
