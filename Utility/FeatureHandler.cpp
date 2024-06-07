#include "FeatureHandler.h"

soap::Species FeatureHandler::get_species() {

	switch (this->type_)
	{
	case DataType::BulkRna:
		return static_cast<BulkRna*>(this->data_)->species_;
		break;
	case DataType::SingleCellRna:
		return static_cast<SingleCellRna*>(this->data_)->species_;
		break;
	case DataType::SingleCellAtac:
		return static_cast<SingleCellAtac*>(this->data_)->species_;
		break;
	case DataType::SingleCellMultiome:
		return static_cast<SingleCellMultiome*>(this->data_)->species_;
		break;
		break;
	default:
		break;
	}

	return {};
};

QStringList FeatureHandler::get_metadata_names() {

	switch (this->type_)
	{
	case DataType::BulkRna:
		return static_cast<BulkRna*>(this->data_)->metadata()->mat_.colnames_;
		break;
	case DataType::SingleCellRna:
		return static_cast<SingleCellRna*>(this->data_)->metadata()->mat_.colnames_;
		break;
	case DataType::SingleCellAtac:
		return static_cast<SingleCellAtac*>(this->data_)->metadata()->mat_.colnames_;
		break;
	case DataType::SingleCellMultiome:
		return static_cast<SingleCellMultiome*>(this->data_)->metadata()->mat_.colnames_;
		break;
		break;
	default:
		break;
	}

	return {};
};

FEATURE_DATA FeatureHandler::get_feature_names() {

	FEATURE_DATA res;

	switch (this->type_)
	{
	case DataType::BulkRna:
		res = get_feature_names_bulk_rna();
		break;
	case DataType::SingleCellRna:
		res = get_feature_names_single_cell_rna();
		break;
	case DataType::SingleCellAtac:
		res = get_feature_names_single_cell_atac();
		break;
	case DataType::SingleCellMultiome:
		res = get_feature_names_single_cell_multiome();
		break; 
	case DataType::DataFrame:
		res = get_feature_names_dataframe();
		break;
	default:
		break;
	}

	res.all_names = _Cs unique(res.all_names);

	res.factor_names = _Cs unique(res.factor_names);

	res.completer_names = _Cs unique(res.completer_names);

	return res;
};


FEATURE_DATA FeatureHandler::get_feature_names_bulk_rna() {

	FEATURE_DATA res;

	auto rna = static_cast<BulkRna*>(this->data_);

	auto metadata = rna->metadata();

	if (metadata != nullptr) {

		this->check_custom_matrix_names(metadata->mat_, res);
	}

	auto counts = rna->counts();

	if (counts != nullptr) {

		res.all_names << counts->rownames_;
		res.completer_names << counts->rownames_;

		return res;
	}

	auto normalized = rna->normalized();

	if (normalized != nullptr) {
		res.all_names << normalized->rownames_;
		res.completer_names << normalized->rownames_;

		return res;
	}

	return res;
};

FEATURE_DATA FeatureHandler::get_feature_names_single_cell_rna() {

	FEATURE_DATA res;

	auto single_cell_rna = static_cast<SingleCellRna*>(this->data_);

	auto metadata = single_cell_rna->metadata();

	if (metadata != nullptr) {

		this->check_custom_matrix_names(metadata->mat_, res);
	}

	auto counts = single_cell_rna->counts();

	if (counts != nullptr) {

		res.all_names << counts->rownames_;
		res.completer_names << counts->rownames_;

		return res;
	}

	auto normalized = single_cell_rna->normalized();

	if (normalized != nullptr) {
		res.all_names << normalized->rownames_;
		res.completer_names << normalized->rownames_;

		return res;
	}

	return res;
};

FEATURE_DATA FeatureHandler::get_feature_names_single_cell_atac() {

	FEATURE_DATA res;

	auto single_cell_atac = static_cast<SingleCellAtac*>(this->data_);

	auto metadata = single_cell_atac->metadata();

	if (metadata != nullptr) {

		this->check_custom_matrix_names(metadata->mat_, res);
	}

	auto counts = single_cell_atac->counts();

	if (counts != nullptr) {

		res.all_names << counts->rownames_;
	}
	else {
		auto normalized = single_cell_atac->normalized();

		if (normalized != nullptr) {

			res.all_names << normalized->rownames_;
		}
	}

	auto gene_activity_counts = single_cell_atac->gene_activity_counts();

	if (gene_activity_counts != nullptr) {

		res.all_names << gene_activity_counts->rownames_;
		res.completer_names << gene_activity_counts->rownames_;
	}
	else {
		auto gene_activity_normalized = single_cell_atac->gene_activity_normalized();

		if (gene_activity_normalized != nullptr) {

			res.all_names << gene_activity_normalized->rownames_;
			res.completer_names << gene_activity_normalized->rownames_;
		}
	}
	return res;
};

FEATURE_DATA FeatureHandler::get_feature_names_dataframe() {
	FEATURE_DATA res;

	auto data_frame = static_cast<DataFrame*>(this->data_);

	this->check_custom_matrix_names(data_frame->mat_, res);

	return res;
};

FEATURE_DATA FeatureHandler::get_feature_names_single_cell_multiome() {
	FEATURE_DATA res;

	auto single_cell_multiome = static_cast<SingleCellMultiome*>(this->data_);

	auto metadata = single_cell_multiome->metadata();

	if (metadata != nullptr) {

		this->check_custom_matrix_names(metadata->mat_, res);
	}

	auto rna_counts = single_cell_multiome->rna_counts();

	if (rna_counts != nullptr) {

		res.all_names << rna_counts->rownames_;
		res.completer_names << rna_counts->rownames_;
	}
	else {
		auto rna_normalized = single_cell_multiome->rna_normalized();

		if (rna_normalized != nullptr) {

			res.all_names << rna_normalized->rownames_;
			res.completer_names << rna_normalized->rownames_;
		}
	}

	auto atac_counts = single_cell_multiome->atac_counts();

	if (atac_counts != nullptr) {

		res.all_names << atac_counts->rownames_;
	}
	else {
		auto atac_normalized = single_cell_multiome->atac_normalized();

		if (atac_normalized != nullptr) {

			res.all_names << atac_normalized->rownames_;
		}
	}

	auto gene_activity_counts = single_cell_multiome->trans_counts();

	if (gene_activity_counts != nullptr) {

		res.all_names << gene_activity_counts->rownames_;
		res.completer_names << gene_activity_counts->rownames_;
	}
	else {
		auto gene_activity_normalized = single_cell_multiome->trans_normalized();

		if (gene_activity_normalized != nullptr) {

			res.all_names << gene_activity_normalized->rownames_;
			res.completer_names << gene_activity_normalized->rownames_;
		}
	}
	return res;
};

void FeatureHandler::check_custom_matrix_names(CustomMatrix& mat, FEATURE_DATA& res) {

	for (auto&& [name, type] : mat.data_type_) {

		switch (type)
		{
		case CustomMatrix::DataType::QStringFactor:
		case CustomMatrix::DataType::IntegerFactor:
			res.factor_names << name;
		case CustomMatrix::DataType::DoubleNumeric:
		case CustomMatrix::DataType::IntegerNumeric:
			res.all_names << name;
			res.completer_names << name;
			break;
		default:
			break;
		}
	}

};

void FeatureHandler::check_custom_matrix(CustomMatrix& mat, QUERY_DATA& res, const QString& name) {

	auto type = mat.data_type_.at(name);

	switch (type)
	{
	case CustomMatrix::DataType::DoubleNumeric:
		res.type = QUERY_DATA::DataType::numeric;
		res.dd = mat.get_const_double_reference(name);
		break;
	case CustomMatrix::DataType::IntegerNumeric:
		res.type = QUERY_DATA::DataType::integer;
		res.di = mat.get_const_integer_reference(name);
		break;
	case CustomMatrix::DataType::QString:
		res.type = QUERY_DATA::DataType::string;
		res.ds = mat.get_const_qstring_reference(name);
		break;
	case CustomMatrix::DataType::QStringFactor:
		res.type = QUERY_DATA::DataType::string;
		res.ds = mat.get_const_qstring_reference(name);
		res.dsl = mat.string_factors_.at(name);
		break;
	case CustomMatrix::DataType::IntegerFactor:
		res.type = QUERY_DATA::DataType::integer;
		res.di = mat.get_const_integer_reference(name);
		res.dil = mat.integer_factors_.at(name);
		break;
	default:
		break;
	}
};

QUERY_DATA FeatureHandler::get_data(QUERY_INFO info) {

	switch (this->type_)
	{
	case DataType::BulkRna:
		return get_data_bulk_rna(info);
		break;
	case DataType::SingleCellRna:
		return get_data_single_cell_rna(info);
		break;
	case DataType::SingleCellAtac:
		return get_data_single_cell_atac(info);
		break;
	case DataType::SingleCellMultiome:
		return get_data_single_cell_multiome(info);
		break;
	case DataType::DataFrame:
		return get_data_single_cell_dataframe(info);
		break;
	default:
		break;
	}

	return {};
};

QUERY_DATA FeatureHandler::get_data_single_cell_atac(QUERY_INFO info) {

	QUERY_DATA res{ info.name };

	auto single_cell_atac = static_cast<SingleCellAtac*>(this->data_);

	auto metadata = single_cell_atac->metadata();

	if (metadata != nullptr) {
		if (metadata->mat_.contains(info.name)) {

			this->check_custom_matrix(metadata->mat_, res, info.name);
			res.info["Source"] = "Metadata";
			return res;
		}
	}

	if (info.gene_activity) {

		res.info["Source"] = "Gene Activity";

		if (info.normalize) {
			auto gene_activity_normalized = single_cell_atac->gene_activity_normalized();

			if (gene_activity_normalized == nullptr) {
				return res;
			}

			auto index = gene_activity_normalized->rownames_.indexOf(info.name);

			if (index == -1) {
				return res;
			}

			res.type = QUERY_DATA::DataType::numeric;
			res.dd = _Cs cast<QVector>(Eigen::ArrayXd(gene_activity_normalized->mat_.row(index)));
			return res;
		}
		else {
			auto gene_activity_counts = single_cell_atac->gene_activity_counts();

			if (gene_activity_counts == nullptr) {
				return res;
			}

			auto index = gene_activity_counts->rownames_.indexOf(info.name);

			if (index == -1) {
				return res;
			}

			res.type = QUERY_DATA::DataType::integer;
			res.di = _Cs cast<QVector>(Eigen::ArrayXi(gene_activity_counts->mat_.row(index)));
			return res;
		}
	}
	else {

		res.info["Source"] = "ATAC";

		if (info.normalize) {
			auto normalized = single_cell_atac->normalized();

			if (normalized == nullptr) {
				return res;
			}

			auto index = normalized->rownames_.indexOf(info.name);

			if (index == -1) {
				return res;
			}

			res.type = QUERY_DATA::DataType::numeric;
			res.dd = _Cs cast<QVector>(Eigen::ArrayXd(normalized->mat_.row(index)));
			return res;
		}
		else {
			auto counts = single_cell_atac->counts();

			if (counts == nullptr) {
				return res;
			}

			auto index = counts->rownames_.indexOf(info.name);

			if (index == -1) {
				return res;
			}

			res.type = QUERY_DATA::DataType::integer;
			res.di = _Cs cast<QVector>(Eigen::ArrayXi(counts->mat_.row(index)));
			return res;
		}
	}
}

QUERY_DATA FeatureHandler::get_data_single_cell_dataframe(QUERY_INFO info) {

	QUERY_DATA res{ info.name };

	auto data_frame = static_cast<DataFrame*>(this->data_);

	if (data_frame->mat_.contains(info.name)) {

		this->check_custom_matrix(data_frame->mat_, res, info.name);
	}

	return res;
};

QUERY_DATA FeatureHandler::get_data_single_cell_multiome(QUERY_INFO info) {
	QUERY_DATA res{ info.name };

	auto single_cell_multiome = static_cast<SingleCellMultiome*>(this->data_);

	auto metadata = single_cell_multiome->metadata();

	if (metadata != nullptr) {
		if (metadata->mat_.contains(info.name)) {

			this->check_custom_matrix(metadata->mat_, res, info.name);
			res.info["Source"] = "Metadata";
			return res;
		}
	}

	if (info.gene_activity) {

		res.info["Source"] = "Gene Activity";

		if (info.normalize) {
			auto trans_normalized = single_cell_multiome->trans_normalized();

			if (trans_normalized == nullptr) {
				return res;
			}

			auto index = trans_normalized->rownames_.indexOf(info.name);

			if (index == -1) {
				return res;
			}

			res.type = QUERY_DATA::DataType::numeric;
			res.dd = _Cs cast<QVector>(Eigen::ArrayXd(trans_normalized->mat_.row(index)));
			return res;
		}
		else {
			auto trans_counts = single_cell_multiome->trans_counts();

			if (trans_counts == nullptr) {
				return res;
			}

			auto index = trans_counts->rownames_.indexOf(info.name);

			if (index == -1) {
				return res;
			}

			res.type = QUERY_DATA::DataType::integer;
			res.di = _Cs cast<QVector>(Eigen::ArrayXi(trans_counts->mat_.row(index)));
			return res;
		}
	}

	for (auto&& [_, df] : SUBMODULES(*single_cell_multiome, DataField)) {

		if (info.normalize) {
			auto normalized = df.normalized();

			if (normalized == nullptr) {
				continue;
			}

			auto index = normalized->rownames_.indexOf(info.name);

			if (index == -1) {
				continue;
			}

			if (df.data_type_ == DataField::DataType::Rna) {
				res.info["Source"] = "RNA";
			}
			if (df.data_type_ == DataField::DataType::Atac) {
				res.info["Source"] = "ATAC";
			}

			res.type = QUERY_DATA::DataType::numeric;
			res.dd = _Cs cast<QVector>(Eigen::ArrayXd(normalized->mat_.row(index)));
			return res;
		}
		else {
			auto counts = df.counts();

			if (counts == nullptr) {
				continue;
			}

			auto index = counts->rownames_.indexOf(info.name);

			if (index == -1) {
				continue;
			}

			if (df.data_type_ == DataField::DataType::Rna) {
				res.info["Source"] = "RNA";
			}
			if (df.data_type_ == DataField::DataType::Atac) {
				res.info["Source"] = "ATAC";
			}

			res.type = QUERY_DATA::DataType::integer;
			res.di = _Cs cast<QVector>(Eigen::ArrayXi(counts->mat_.row(index)));
			return res;
		}
	}

	return res;
}


QUERY_DATA FeatureHandler::get_data_bulk_rna(QUERY_INFO info) {

	QUERY_DATA res{ info.name };

	auto single_cell_rna = static_cast<BulkRna*>(this->data_);

	auto metadata = single_cell_rna->metadata();

	if (metadata != nullptr) {
		if (metadata->mat_.contains(info.name)) {

			this->check_custom_matrix(metadata->mat_, res, info.name);

			res.info["Source"] = "Metadata";

			return res;
		}
	}

	res.info["Source"] = "RNA";

	if (info.normalize) {
		auto normalized = single_cell_rna->normalized();

		if (normalized == nullptr) {
			return res;
		}

		auto index = normalized->rownames_.indexOf(info.name);

		if (index == -1) {
			return res;
		}
		res.type = QUERY_DATA::DataType::numeric;
		res.dd = _Cs cast<QVector>(Eigen::ArrayXd(normalized->mat_.row(index)));
		return res;
	}
	else {
		auto counts = single_cell_rna->counts();

		if (counts == nullptr) {
			return res;
		}

		auto index = counts->rownames_.indexOf(info.name);

		if (index == -1) {
			return res;
		}

		res.type = QUERY_DATA::DataType::integer;
		res.di = _Cs cast<QVector>(Eigen::ArrayXi(counts->mat_.row(index)));
		return res;
	}
};

QUERY_DATA FeatureHandler::get_data_single_cell_rna(QUERY_INFO info) {

	QUERY_DATA res{ info.name };

	auto single_cell_rna = static_cast<SingleCellRna*>(this->data_);

	auto metadata = single_cell_rna->metadata();

	if (metadata != nullptr) {
		if (metadata->mat_.contains(info.name)) {

			this->check_custom_matrix(metadata->mat_, res, info.name);

			res.info["Source"] = "Metadata";

			return res;
		}
	}

	res.info["Source"] = "RNA";

	if (info.normalize) {
		auto normalized = single_cell_rna->normalized();

		if (normalized == nullptr) {
			return res;
		}

		auto index = normalized->rownames_.indexOf(info.name);

		if (index == -1) {
			return res;
		}
		res.type = QUERY_DATA::DataType::numeric;
		res.dd = _Cs cast<QVector>(Eigen::ArrayXd(normalized->mat_.row(index)));
		return res;
	}
	else {
		auto counts = single_cell_rna->counts();

		if (counts == nullptr) {
			return res;
		}

		auto index = counts->rownames_.indexOf(info.name);

		if (index == -1) {
			return res;
		}

		res.type = QUERY_DATA::DataType::integer;
		res.di = _Cs cast<QVector>(Eigen::ArrayXi(counts->mat_.row(index)));
		return res;
	}
};