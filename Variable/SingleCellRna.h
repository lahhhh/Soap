#pragma once

#include "Identifier.h"

#include "SparseInt.h"
#include "SparseDouble.h"
#include "DenseDouble.h"
#include "Embedding.h"
#include "Metadata.h"
#include "Gsea.h"
#include "CellChat.h"
#include "Cnv.h"
#include "VelocytoBase.h"
#include "DifferentialAnalysis.h"
#include "Monocle3.h"
#include "Cicero.h"

class SingleCellRna
{
public:

	enum class DataType : int { Plain = 0, Integrated = 1 };

	G_CLASS_FUNCTION_DEFAULT(SingleCellRna);

	G_SET_IDENTIFIER("SingleCellRna");

	G_QUICK_ACCESS(Metadata, metadata);
	G_QUICK_ACCESS_TYPE(SparseInt, counts, Counts);
	G_QUICK_ACCESS_TYPE(SparseDouble, normalized, Normalized);
	G_QUICK_ACCESS_TYPE(Embedding, pca, Pca);
	G_QUICK_ACCESS_TYPE(Embedding, umap, Umap);
	G_QUICK_ACCESS_TYPE(Embedding, tsne, Tsne);
	G_QUICK_ACCESS_TYPE(Embedding, harmony, Harmony);
	G_QUICK_ACCESS(VelocytoBase, velocyto_base);

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	void row_slice(const SliceType& slice) {

		for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
			data.row_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, SparseDouble)) {
			data.row_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, VelocytoBase)) {
			data.row_slice(slice);
		}

		SUBMODULES(*this, DifferentialAnalysis).clear();
		SUBMODULES(*this, GSEA).clear();
		SUBMODULES(*this, CellChat).clear();
		SUBMODULES(*this, CNV).clear();
		SUBMODULES(*this, Monocle3).clear();
	}

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	[[nodiscard]] SingleCellRna* row_sliced(const SliceType& slice) const {

		SingleCellRna* ret = new SingleCellRna(*this);

		ret->row_slice(slice);

		return ret;
	}

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	void col_slice(const SliceType& slice) {

		for (auto& [name, data] : SUBMODULES(*this, Metadata)) {
			data.row_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
			data.col_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, SparseDouble)) {
			data.col_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, Embedding)) {
			data.slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, VelocytoBase)) {
			data.col_slice(slice);
		}

		SUBMODULES(*this, DifferentialAnalysis).clear();
		SUBMODULES(*this, GSEA).clear();
		SUBMODULES(*this, CellChat).clear();
		SUBMODULES(*this, CNV).clear();
		SUBMODULES(*this, Monocle3).clear();
	}

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	[[nodiscard]] SingleCellRna* col_sliced(const SliceType& slice) const {

		SingleCellRna* ret = new SingleCellRna(*this);

		ret->col_slice(slice);

		return ret;
	}

	template <typename OrderType>
		requires _Cs is_order_container<OrderType>
	void col_reorder(const OrderType& order) {

		for (auto& [name, data] : SUBMODULES(*this, Metadata)) {
			data.row_reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
			data.col_reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, SparseDouble)) {
			data.col_reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, Embedding)) {
			data.reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, VelocytoBase)) {
			data.col_reorder(order);
		}

		SUBMODULES(*this, DifferentialAnalysis).clear();
		SUBMODULES(*this, GSEA).clear();
		SUBMODULES(*this, CellChat).clear();
		SUBMODULES(*this, CNV).clear();
		SUBMODULES(*this, Monocle3).clear();
	}

	template <typename OrderType>
		requires _Cs is_order_container<OrderType>
	[[nodiscard]] SingleCellRna* col_reordered(const OrderType& order) const {

		SingleCellRna* ret = new SingleCellRna(*this);

		ret->col_reorder(order);

		return ret;
	}

	template <typename SliceType, typename SliceType2>
		requires _Cs is_slice_container<SliceType>
	void slice(const SliceType& row_slice, const SliceType2& col_slice) {

		this->row_slice(row_slice);

		this->col_slice(col_slice);
	}

	template <typename SliceType, typename SliceType2>
		requires _Cs is_slice_container<SliceType>
	[[nodiscard]] SingleCellRna* sliced(const SliceType& row_slice, const SliceType2& col_slice) const {

		SingleCellRna* ret = new SingleCellRna(*this);

		ret->slice(row_slice, col_slice);

		return ret;
	}

	/******       data        *****/

	DataType data_type_{ DataType::Plain };

	soap::Species species_ = soap::Species::Undefined;
	int random_state_ = 1997;

	std::map<QString, QString> string_information_;
	std::map<QString, int> integer_information_;
	std::map<QString, double> double_information_;

	std::map<QString, QStringList> string_vectors_;
	std::map<QString, QVector<int> > integer_vectors_;
	std::map<QString, QVector<double> > double_vectors_;	

	SOAP_SUBMODULES(Metadata);
	SOAP_SUBMODULES(SparseInt);
	SOAP_SUBMODULES(SparseDouble);
	SOAP_SUBMODULES(Embedding);
	SOAP_SUBMODULES(VelocytoBase);
	SOAP_SUBMODULES(DifferentialAnalysis);
	SOAP_SUBMODULES(GSEA);
	SOAP_SUBMODULES(CellChat);
	SOAP_SUBMODULES(CNV);	
	SOAP_SUBMODULES(Monocle3);
};

