#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"

class Metadata
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(Metadata);

	template<typename Mat> requires (!std::same_as<std::decay_t<Mat>, Metadata>)
		Metadata(Mat&& df) : mat_(std::forward<Mat>(df)) {}

	template <typename SliceType>
		requires custom::is_slice_container<SliceType>
	void row_slice(const SliceType& slice) {
		this->mat_.row_slice(slice);
	}

	template <typename SliceType>
		requires custom::is_slice_container<SliceType>
	[[nodiscard]] Metadata row_sliced(const SliceType& slice) const {
		return Metadata(this->mat_.row_sliced(slice));
	}

	template <typename OrderType>
		requires custom::is_order_container<OrderType>
	void row_reorder(const OrderType& order) {
		this->mat_.row_reorder(order);
	}

	template <typename OrderType>
		requires custom::is_order_container<OrderType>
	[[nodiscard]] Metadata row_reordered(const OrderType& order) const {
		return Metadata(this->mat_.row_reorder(order));
	}

	DataType data_type_{ DataType::Plain };

	CustomMatrix mat_;

	G_SET_IDENTIFIER("Metadata");
};

