#pragma once

#include "Identifier.h"

#include "DataFrame.h"
#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"
#include "BulkRna.h"

#include "custom.h"

struct FEATURE_DATA {

	QStringList all_names;

	QStringList factor_names;

	QStringList numeric_names;
};

struct QUERY_DATA {

	QString name;

	enum class DataType : int {notype = 0, string = 1, integer = 2, numeric = 3};

	DataType type{ DataType::notype};

	QVector<int> di;

	QVector<double> dd;

	QStringList ds;

	QStringList dsl;

	QVector<int> dil;

	QMap<QString, QString> info;
	QSet<QString> message;

	qsizetype length() const { return _Cs max(di.size(), dd.size(), ds.size()); };
	bool is_valid() const { return type != DataType::notype; };
	bool is_continuous() const { return type == DataType::numeric || (type == DataType::integer && dil.isEmpty()); }
	QVector<double> get_continuous() const {

		if (type == DataType::numeric) {
			return dd;
		}

		if (type == DataType::integer && dil.isEmpty()) {
			return _Cs cast<double>(di);
		}

		return {};
	}

	bool is_factor() const { return !(dil.isEmpty() && dsl.isEmpty()); }
	QStringList get_factor() const {

		if (type == DataType::string) {
			return ds;
		}

		if (type == DataType::integer) {
			return _Cs cast<QString>(di);
		}

		return {};
	}

	QStringList get_levels() const {

		if (!dsl.isEmpty()) {
			return dsl;
		}

		if (!dil.isEmpty()) {
			return _Cs cast<QString>(dil);
		}

		return {};
	}

	template<typename SliceType>
	requires _Cs is_slice_container<SliceType>
	void slice(const SliceType& slice) {
		if (!dd.isEmpty()) {
			dd = _Cs sliced(dd, slice);
		}

		if (!di.isEmpty()) {
			di = _Cs sliced(di, slice);
			if (!dil.isEmpty()) {
				dil = _Cs unique(dil);
			}
		}

		if (!ds.isEmpty()) {
			ds = _Cs sliced(ds, slice);
			if (!dsl.isEmpty()) {
				dsl = _Cs unique(dsl);
			}
		}
	}
};

struct QUERY_INFO {

	QString name;

	bool normalize{false};

	bool gene_activity{ false };
};

class FeatureHandler
{
public:

	FeatureHandler() = default;
	explicit FeatureHandler(Metadata* d) : type_(DataType::Metadata), data_(d) {};
	explicit FeatureHandler(DataFrame* d) : type_(DataType::DataFrame), data_(d) {};
	explicit FeatureHandler(SingleCellRna* d) : type_(DataType::SingleCellRna), data_(d) {};
	explicit FeatureHandler(SingleCellAtac* d) : type_(DataType::SingleCellAtac), data_(d) {};
	explicit FeatureHandler(SingleCellMultiome* d) : type_(DataType::SingleCellMultiome), data_(d) {};
	explicit FeatureHandler(BulkRna* d) : type_(DataType::BulkRna), data_(d) {};

	void set(Metadata* d) {
		this->type_ = DataType::Metadata;
		this->data_ = d;
	};

	void set(DataFrame* d){
		this->type_ = DataType::DataFrame;
		this->data_ = d;
	};

	void set(SingleCellRna* d){
		this->type_ = DataType::SingleCellRna;
		this->data_ = d;
	};

	void set(BulkRna* d) {
		this->type_ = DataType::BulkRna;
		this->data_ = d;
	};

	void set(SingleCellAtac* d){
		this->type_ = DataType::SingleCellAtac;
		this->data_ = d;
	};

	void set(SingleCellMultiome* d){
		this->type_ = DataType::SingleCellMultiome;
		this->data_ = d;
	};

	QUERY_DATA get_data(QUERY_INFO);

	FEATURE_DATA get_feature_names();

	QStringList get_metadata_names();

	soap::Species get_species();

	enum class DataType : int { NoType, SingleCellRna, SingleCellAtac, SingleCellMultiome, DataFrame,
		BulkRna, Metadata};

	DataType type_{ DataType::NoType };

	void* data_{ nullptr };

private:

	void check_custom_matrix_names(CustomMatrix& mat, FEATURE_DATA& res);
	void check_custom_matrix(CustomMatrix& mat, QUERY_DATA&, const QString& name);

	QUERY_DATA get_data_bulk_rna(QUERY_INFO);
	QUERY_DATA get_data_single_cell_rna(QUERY_INFO);
	QUERY_DATA get_data_single_cell_atac(QUERY_INFO);
	QUERY_DATA get_data_single_cell_multiome(QUERY_INFO);
	QUERY_DATA get_data_dataframe(QUERY_INFO);
	QUERY_DATA get_data_metadata(QUERY_INFO);

	FEATURE_DATA get_feature_names_bulk_rna();
	FEATURE_DATA get_feature_names_single_cell_rna();
	FEATURE_DATA get_feature_names_single_cell_atac();
	FEATURE_DATA get_feature_names_single_cell_multiome();
	FEATURE_DATA get_feature_names_dataframe();
	FEATURE_DATA get_feature_names_metadata();
};

