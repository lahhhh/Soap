#pragma once

#include "Identifier.h"

#include "CustomTemplates.h"

#include "DataField.h"
#include "SparseInt.h"
#include "SparseDouble.h"
#include "DenseDouble.h"
#include "CustomMatrix.h"
#include "Embedding.h"
#include "Metadata.h"
#include "Enrichment.h"
#include "Gsea.h"
#include "CellChat.h"
#include "Cnv.h"
#include "GenomicRange.h"
#include "MotifPosition.h"
#include "CoverageTrack.h"
#include "Fragments.h"
#include "VelocytoBase.h"
#include "Pando.h"
#include "Monocle3.h"
#include "Cicero.h"

class SingleCellMultiome
{

public:

	enum class DataType : int { Plain = 0, Integrated = 1 };

	G_CLASS_FUNCTION_DEFAULT(SingleCellMultiome);

	G_SET_IDENTIFIER("SingleCellMultiome");

	G_QUICK_ACCESS(Metadata, metadata);
	G_QUICK_ACCESS(Fragments, fragments);
	G_QUICK_ACCESS(MotifPosition, motif_position);
	G_QUICK_ACCESS(Cicero, cicero);
	G_QUICK_ACCESS(VelocytoBase, velocyto_base);

	DataField* get_field(DataField::DataType field_type);
	DataField* rna_field();
	DataField* atac_field();
	DataField* trans_field();

	DataField& create_field(const QString& field_name, DataField::DataType field_type);

	SparseInt* field_sparseint(DataField::DataType type1, SparseInt::DataType type2);
	SparseDouble* field_sparsedouble(DataField::DataType type1, SparseDouble::DataType type2);
	Embedding* field_embedding(DataField::DataType type1, Embedding::DataType type2);

	SparseInt const * field_sparseint(DataField::DataType type1, SparseInt::DataType type2) const;
	SparseDouble const * field_sparsedouble(DataField::DataType type1, SparseDouble::DataType type2) const;
	Embedding const * field_embedding(DataField::DataType type1, Embedding::DataType type2) const;

	QStringList embedding_names() const;
	Embedding* embedding(const QString& name);

	SparseInt* rna_counts();
	SparseInt* atac_counts();
	SparseInt* trans_counts();

	SparseDouble* rna_normalized();
	SparseDouble* atac_normalized();
	SparseDouble* trans_normalized();

	Embedding* rna_pca();
	Embedding* atac_pca();
	Embedding* trans_pca();

	Embedding* rna_umap();
	Embedding* atac_umap();
	Embedding* trans_umap();

	Embedding* rna_tsne();
	Embedding* atac_tsne();
	Embedding* trans_tsne();

	Embedding* rna_harmony();
	Embedding* atac_harmony();
	Embedding* trans_harmony();

	SparseInt const * rna_counts() const;
	SparseInt const * atac_counts() const;
	SparseInt const * trans_counts() const;

	SparseDouble const * rna_normalized() const;
	SparseDouble const * atac_normalized() const;
	SparseDouble const * trans_normalized() const;

	Embedding const * rna_pca() const;
	Embedding const * atac_pca() const;
	Embedding const * trans_pca() const;

	Embedding const * rna_umap() const;
	Embedding const * atac_umap() const;
	Embedding const * trans_umap() const;

	Embedding const * rna_tsne() const;
	Embedding const * atac_tsne() const;
	Embedding const * trans_tsne() const;

	Embedding const * rna_harmony() const;
	Embedding const * atac_harmony() const;
	Embedding const * trans_harmony() const;

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	void col_slice(const SliceType& slice) {

		for (auto& [name, data] : SUBMODULES(*this, Metadata)) {
			data.row_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, DataField)) {
			data.col_slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, Fragments)) {
			data.slice(slice);
		}

		for (auto& [name, data] : SUBMODULES(*this, VelocytoBase)) {
			data.col_slice(slice);
		}

		SUBMODULES(*this, GSEA).clear();
		SUBMODULES(*this, CellChat).clear();
		SUBMODULES(*this, CNV).clear();
		SUBMODULES(*this, GenomicRange).clear();
		SUBMODULES(*this, MotifPosition).clear();
		SUBMODULES(*this, CoverageTrack).clear();
		SUBMODULES(*this, Pando).clear();
		SUBMODULES(*this, Monocle3).clear();
		SUBMODULES(*this, Cicero).clear();
	}

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	[[nodiscard]] SingleCellMultiome* col_sliced(const SliceType& slice) const {

		SingleCellMultiome* ret = new SingleCellMultiome(*this);

		ret->col_slice(slice);

		return ret;
	}

	template <typename OrderType>
		requires _Cs is_order_container<OrderType>
	void col_reorder(const OrderType& order) {

		for (auto& [name, data] : SUBMODULES(*this, Metadata)) {
			data.row_reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, DataField)) {
			data.col_reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, Fragments)) {
			data.reorder(order);
		}

		for (auto& [name, data] : SUBMODULES(*this, VelocytoBase)) {
			data.col_reorder(order);
		}

		SUBMODULES(*this, GSEA).clear();
		SUBMODULES(*this, CellChat).clear();
		SUBMODULES(*this, CNV).clear();
		SUBMODULES(*this, GenomicRange).clear();
		SUBMODULES(*this, MotifPosition).clear();
		SUBMODULES(*this, CoverageTrack).clear();
		SUBMODULES(*this, Pando).clear();
		SUBMODULES(*this, Monocle3).clear();
		SUBMODULES(*this, Cicero).clear();
	}

	template <typename OrderType>
		requires _Cs is_order_container<OrderType>
	[[nodiscard]] SingleCellMultiome* col_reordered(const OrderType& order) const {

		SingleCellMultiome* ret = new SingleCellMultiome(*this);

		ret->col_reorder(order);

		return ret;
	}

	/*********      data       *********/

	DataType data_type_{ DataType::Plain };

	soap::Species species_ = soap::Species::Undefined;
	int random_state_ = 1997;
	Bias insertion_bias_;

	std::map<QString, QString> string_information_;
	std::map<QString, int> integer_information_;
	std::map<QString, double> double_information_;

	std::map<QString, QStringList> string_vectors_;
	std::map<QString, QVector<int>> integer_vectors_;
	std::map<QString, QVector<double>> double_vectors_;

	SOAP_SUBMODULES(DataField);
	SOAP_SUBMODULES(Metadata);
	SOAP_SUBMODULES(GSEA);
	SOAP_SUBMODULES(CellChat);
	SOAP_SUBMODULES(CNV);
	SOAP_SUBMODULES(GenomicRange);
	SOAP_SUBMODULES(MotifPosition);
	SOAP_SUBMODULES(CoverageTrack);
	SOAP_SUBMODULES(VelocytoBase);
	SOAP_SUBMODULES(Fragments);
	SOAP_SUBMODULES(Pando);
	SOAP_SUBMODULES(Monocle3);
	SOAP_SUBMODULES(Cicero);
};
