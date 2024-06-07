#pragma once

#include "Identifier.h"

#include "SparseInt.h"
#include "SparseDouble.h"
#include "DenseDouble.h"
#include "CustomMatrix.h"
#include "Embedding.h"
#include "Metadata.h"
#include "Enrichment.h"
#include "DifferentialAnalysis.h"
#include "GenomicRange.h"
#include "MotifPosition.h"
#include "CoverageTrack.h"
#include "Fragments.h"
#include "Cicero.h"
#include "Monocle3.h"

class SingleCellAtac
{
public:

	enum class DataType : int { Plain = 0, Integrated = 1 };

	G_CLASS_FUNCTION_DEFAULT(SingleCellAtac);

	G_SET_IDENTIFIER("SingleCellAtac");

	G_QUICK_ACCESS(Metadata, metadata);
	G_QUICK_ACCESS_TYPE(SparseInt, counts, Counts);
	G_QUICK_ACCESS_TYPE(SparseDouble, normalized, Normalized);
	G_QUICK_ACCESS_TYPE(Embedding, pca, Pca);
	G_QUICK_ACCESS_TYPE(Embedding, umap, Umap);
	G_QUICK_ACCESS_TYPE(Embedding, tsne, Tsne);
	G_QUICK_ACCESS_TYPE(Embedding, harmony, Harmony);
	G_QUICK_ACCESS_TYPE(SparseInt, gene_activity_counts, GeneActivity);
	G_QUICK_ACCESS_TYPE(SparseDouble, gene_activity_normalized, GeneActivity);
	G_QUICK_ACCESS(Fragments, fragments);
	G_QUICK_ACCESS(MotifPosition, motif_position);
	G_QUICK_ACCESS(Cicero, cicero);

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

		for (auto& [name, data] : SUBMODULES(*this, Fragments)) {
			data.slice(slice);
		}

		SUBMODULES(*this, DifferentialAnalysis).clear();
		SUBMODULES(*this, GenomicRange).clear();
		SUBMODULES(*this, MotifPosition).clear();
		SUBMODULES(*this, CoverageTrack).clear();
		SUBMODULES(*this, Cicero).clear();
		SUBMODULES(*this, Monocle3).clear();
	}

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	[[nodiscard]] SingleCellAtac* col_sliced(const SliceType& slice) const {

		SingleCellAtac* ret = new SingleCellAtac(*this);

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

		for (auto& [name, data] : SUBMODULES(*this, Fragments)) {
			data.reorder(order);
		}

		SUBMODULES(*this, DifferentialAnalysis).clear();
		SUBMODULES(*this, GenomicRange).clear();
		SUBMODULES(*this, MotifPosition).clear();
		SUBMODULES(*this, CoverageTrack).clear();
		SUBMODULES(*this, Cicero).clear();
		SUBMODULES(*this, Monocle3).clear();
	}

	template <typename OrderType>
		requires _Cs is_order_container<OrderType>
	[[nodiscard]] SingleCellAtac* col_reordered(const OrderType& order) const {

		SingleCellAtac* ret = new SingleCellAtac(*this);

		ret->col_reorder(order);

		return ret;
	}

	/******       data        *****/

	DataType data_type_{ DataType::Plain };

	soap::Species species_ = soap::Species::Undefined;
	int random_state_{ 1997 };
	Bias insertion_bias_;

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
	SOAP_SUBMODULES(DifferentialAnalysis);
	SOAP_SUBMODULES(GenomicRange);
	SOAP_SUBMODULES(MotifPosition);
	SOAP_SUBMODULES(CoverageTrack);
	SOAP_SUBMODULES(Fragments);
	SOAP_SUBMODULES(Cicero);
	SOAP_SUBMODULES(Monocle3);
};

