#include "SingleCellMultiome.h"

#include <QFileInfo>

DataField* SingleCellMultiome::get_field(DataField::DataType field_type) {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == field_type) {

			return &data;
		}
	}

	return nullptr;
};

DataField* SingleCellMultiome::rna_field() {

	return this->get_field(DataField::DataType::Rna);
};

DataField* SingleCellMultiome::atac_field() {

	return this->get_field(DataField::DataType::Atac);
};

DataField* SingleCellMultiome::trans_field() {

	return this->get_field(DataField::DataType::Trans);
};

DataField& SingleCellMultiome::create_field(const QString& field_name, DataField::DataType field_type) {

	auto& field = SUBMODULES(*this, DataField)[field_name];

	field.data_type_ = field_type;

	return field;
};

SparseInt* SingleCellMultiome::field_sparseint(DataField::DataType type1, SparseInt::DataType type2) {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == type1) {

			for (auto&& [sparseint_name, sparseint_data] : SUBMODULES(data, SparseInt)) {

				if (sparseint_data.data_type_ == type2) {

					return &sparseint_data;
				}
			}
		}
	}

	return nullptr;
};

SparseInt const * SingleCellMultiome::field_sparseint(DataField::DataType type1, SparseInt::DataType type2) const {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == type1) {

			for (auto&& [sparseint_name, sparseint_data] : SUBMODULES(data, SparseInt)) {

				if (sparseint_data.data_type_ == type2) {

					return &sparseint_data;
				}
			}
		}
	}

	return nullptr;
};

SparseInt* SingleCellMultiome::rna_counts() {

	return this->field_sparseint(DataField::DataType::Rna, SparseInt::DataType::Counts);
};

SparseInt* SingleCellMultiome::atac_counts() {

	return this->field_sparseint(DataField::DataType::Atac, SparseInt::DataType::Counts);
};

SparseInt* SingleCellMultiome::trans_counts() {

	return this->field_sparseint(DataField::DataType::Trans, SparseInt::DataType::Counts);
};

SparseInt const * SingleCellMultiome::rna_counts() const {

	return this->field_sparseint(DataField::DataType::Rna, SparseInt::DataType::Counts);
};

SparseInt const * SingleCellMultiome::atac_counts() const {

	return this->field_sparseint(DataField::DataType::Atac, SparseInt::DataType::Counts);
};

SparseInt const * SingleCellMultiome::trans_counts() const {

	return this->field_sparseint(DataField::DataType::Trans, SparseInt::DataType::Counts);
};

SparseDouble* SingleCellMultiome::field_sparsedouble(DataField::DataType type1, SparseDouble::DataType type2) {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == type1) {

			for (auto&& [sparsedouble_name, sparsedouble_data] : SUBMODULES(data, SparseDouble)) {

				if (sparsedouble_data.data_type_ == type2) {

					return &sparsedouble_data;
				}
			}
		}
	}

	return nullptr;
};

SparseDouble const * SingleCellMultiome::field_sparsedouble(DataField::DataType type1, SparseDouble::DataType type2) const {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == type1) {

			for (auto&& [sparsedouble_name, sparsedouble_data] : SUBMODULES(data, SparseDouble)) {

				if (sparsedouble_data.data_type_ == type2) {

					return &sparsedouble_data;
				}
			}
		}
	}

	return nullptr;
};

SparseDouble* SingleCellMultiome::rna_normalized() {

	return this->field_sparsedouble(DataField::DataType::Rna, SparseDouble::DataType::Normalized);
};

SparseDouble* SingleCellMultiome::atac_normalized() {

	return this->field_sparsedouble(DataField::DataType::Atac, SparseDouble::DataType::Normalized);
};

SparseDouble* SingleCellMultiome::trans_normalized() {

	return this->field_sparsedouble(DataField::DataType::Trans, SparseDouble::DataType::Normalized);
};

SparseDouble const * SingleCellMultiome::rna_normalized() const {

	return this->field_sparsedouble(DataField::DataType::Rna, SparseDouble::DataType::Normalized);
};

SparseDouble const * SingleCellMultiome::atac_normalized() const {

	return this->field_sparsedouble(DataField::DataType::Atac, SparseDouble::DataType::Normalized);
};

SparseDouble const * SingleCellMultiome::trans_normalized() const {

	return this->field_sparsedouble(DataField::DataType::Trans, SparseDouble::DataType::Normalized);
};

Embedding* SingleCellMultiome::embedding(const QString& name) {

	for (auto&& [_, data] : SUBMODULES(*this, DataField)) {
		for (auto&& [embedding_name, embedding_data] : SUBMODULES(data, Embedding)) {

			if (embedding_name == name) {
				return &embedding_data;
			}
		}
	}

	return nullptr;
};

QStringList SingleCellMultiome::embedding_names() const{

	QStringList res;

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {
		for (auto&& [embedding_name, _] : SUBMODULES(data, Embedding)) {
			res << embedding_name;
		}
	}

	return res;
};

Embedding* SingleCellMultiome::field_embedding(DataField::DataType type1, Embedding::DataType type2) {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == type1) {

			for (auto&& [embedding_name, embedding_data] : SUBMODULES(data, Embedding)) {

				if (embedding_data.data_type_ == type2) {

					return &embedding_data;
				}
			}
		}
	}

	return nullptr;
};

Embedding const * SingleCellMultiome::field_embedding(DataField::DataType type1, Embedding::DataType type2) const {

	for (auto&& [name, data] : SUBMODULES(*this, DataField)) {

		if (data.data_type_ == type1) {

			for (auto&& [embedding_name, embedding_data] : SUBMODULES(data, Embedding)) {

				if (embedding_data.data_type_ == type2) {

					return &embedding_data;
				}
			}
		}
	}

	return nullptr;
};

Embedding* SingleCellMultiome::rna_pca() {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Pca);
};

Embedding* SingleCellMultiome::atac_pca() {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Pca);
};

Embedding* SingleCellMultiome::trans_pca() {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Pca);
};

Embedding* SingleCellMultiome::rna_umap() {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Umap);
};

Embedding* SingleCellMultiome::atac_umap() {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Umap);
};

Embedding* SingleCellMultiome::trans_umap() {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Umap);
};

Embedding* SingleCellMultiome::rna_tsne() {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Tsne);
};

Embedding* SingleCellMultiome::atac_tsne() {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Tsne);
};

Embedding* SingleCellMultiome::trans_tsne() {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Tsne);
};

Embedding* SingleCellMultiome::rna_harmony() {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Harmony);
};

Embedding* SingleCellMultiome::atac_harmony() {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Harmony);
};

Embedding* SingleCellMultiome::trans_harmony() {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Harmony);
};


Embedding const * SingleCellMultiome::rna_pca() const {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Pca);
};

Embedding const * SingleCellMultiome::atac_pca() const {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Pca);
};

Embedding const * SingleCellMultiome::trans_pca() const {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Pca);
};

Embedding const * SingleCellMultiome::rna_umap() const {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Umap);
};

Embedding const * SingleCellMultiome::atac_umap() const {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Umap);
};

Embedding const * SingleCellMultiome::trans_umap() const {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Umap);
};

Embedding const * SingleCellMultiome::rna_tsne() const {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Tsne);
};

Embedding const * SingleCellMultiome::atac_tsne() const {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Tsne);
};

Embedding const * SingleCellMultiome::trans_tsne() const {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Tsne);
};

Embedding const * SingleCellMultiome::rna_harmony() const {

	return this->field_embedding(DataField::DataType::Rna, Embedding::DataType::Harmony);
};

Embedding const * SingleCellMultiome::atac_harmony() const {

	return this->field_embedding(DataField::DataType::Atac, Embedding::DataType::Harmony);
};

Embedding const * SingleCellMultiome::trans_harmony() const {

	return this->field_embedding(DataField::DataType::Trans, Embedding::DataType::Harmony);
};