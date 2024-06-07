#include "VelocytoBase.h"

const SparseInt* VelocytoBase::get_spliced() const {

	for (const auto& [name, data] : SUBMODULES(*this, SparseInt)) {
		if (data.data_type_ == SparseInt::DataType::Spliced) {
			return &data;
		}
	}

	return nullptr;
};

const SparseInt* VelocytoBase::get_unspliced() const {

	for (const auto& [name, data] : SUBMODULES(*this, SparseInt)) {
		if (data.data_type_ == SparseInt::DataType::Unspliced) {
			return &data;
		}
	}

	return nullptr;
};
