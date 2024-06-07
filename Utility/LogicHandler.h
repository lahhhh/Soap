#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class LogicHandler
{
public:
	LogicHandler(SingleCellRna*);
	LogicHandler(SingleCellAtac*);
	LogicHandler(SingleCellMultiome*);
	LogicHandler(DataField*);
	LogicHandler(CustomMatrix*);

	const QStringList& data_names();

	enum class SourceType : int {SingleCellRna, SingleCellAtac, SingleCellMultiome, DataField, CustomMatrix};

	enum class DataType : int {NoType, String, StringFactor, IntegerFactor, IntegerNumeric, DoubleNumeric};

	DataType get_type(const QString& name) const;

	std::tuple<bool, DataType, QStringList> get_content(const QString& name) const;

	Eigen::ArrayX<bool> resolve(const QString& s) const;

private:

	Eigen::ArrayX<bool> resolve_triplet(const QString& feature, const QString& comp, const QString& val) const;

	Eigen::ArrayX<bool> resolve_triplet_single_cell_rna(const QString& feature, const QString& comp, const QString& val) const;

	Eigen::ArrayX<bool> resolve_triplet_single_cell_atac(const QString& feature, const QString& comp, const QString& val) const;

	Eigen::ArrayX<bool> resolve_triplet_single_cell_multiome(const QString& feature, const QString& comp, const QString& val) const;

	Eigen::ArrayX<bool> resolve_triplet_cm(const QString& feature, const QString& comp, const QString& val) const;

	Eigen::ArrayX<bool> resolve_triplet_data_field(const QString& feature, const QString& comp, const QString& val) const;

	void process_cm(CustomMatrix* mat);

	void finalize();

	void* source_{nullptr};

	int logic_width_{ 0 };

	SourceType source_type_;

	QStringList data_names_;

	std::unordered_set<QString> string_factor_names_;
	std::unordered_set<QString> integer_factor_names_;
	std::unordered_set<QString> integer_numeric_names_;
	std::unordered_set<QString> double_numeric_names_;

	std::unordered_map<QString, QStringList> string_factor_contents_;
	std::unordered_map<QString, QStringList> integer_factor_contents_;
};

