#pragma once

#include "Identifier.h"

#include <QStringList>

#include "Custom.h"

class DenseDouble;

class SparseDouble
{
public:

    enum class DataType : int { Plain = 0, Normalized = 1, GeneActivity = 2 };

    G_CLASS_FUNCTION_DEFAULT(SparseDouble);

	template <typename SparseMatrix, typename String1, typename String2>
	SparseDouble(
        DataType data_type,
		SparseMatrix&& mat,
		String1&& rownames,
		String2&& colnames
	) :
        data_type_(data_type),
		mat_(std::forward<SparseMatrix>(mat)),
		rownames_(std::forward<String1>(rownames)),
		colnames_(std::forward<String2>(colnames))
	{ }

	G_SET_IDENTIFIER("SparseDouble");

    DataType data_type_{ DataType::Plain };

	Eigen::SparseMatrix<double> mat_;
	QStringList rownames_;
	QStringList colnames_;

	qsizetype rows() const;
	qsizetype cols() const;

    template <typename SliceType>
    requires _Cs is_slice_container<SliceType>
        void row_slice(const SliceType& slice) {
        this->mat_ = _Cs row_sliced(this->mat_, slice);
        this->rownames_ = _Cs sliced(this->rownames_, slice);
    };

    template <typename SliceType>
    requires _Cs is_slice_container<SliceType>
        void col_slice(const SliceType& slice) {
        this->mat_ = _Cs col_sliced(this->mat_, slice);
        this->colnames_ = _Cs sliced(this->colnames_, slice);
    };

    template <typename SliceType, typename SliceType2>
    requires _Cs is_slice_container<SliceType>&& _Cs is_slice_container<SliceType2>
        void slice(const SliceType& row_slice, const SliceType2& col_slice) {
        this->mat_ = _Cs sliced(this->mat_, row_slice, col_slice);
        this->rownames_ = _Cs sliced(this->rownames_, row_slice);
        this->colnames_ = _Cs sliced(this->colnames_, col_slice);
    }

    template <typename SliceType, typename SliceType2>
    requires _Cs is_slice_container<SliceType>&& _Cs is_slice_container<SliceType2>
        [[nodiscard]] SparseDouble sliced(const SliceType& row_slice, const SliceType2& col_slice) const {
        return SparseDouble(
            this->data_type_, 
            _Cs sliced(this->mat_, row_slice, col_slice), 
            _Cs sliced(this->rownames_, row_slice),
            _Cs sliced(this->colnames_, col_slice)
        );
    }

    template <typename SliceType>
    requires _Cs is_slice_container<SliceType>
        [[nodiscard]] SparseDouble row_sliced(const SliceType& slice) const {
        return SparseDouble(
            this->data_type_, 
            _Cs row_sliced(this->mat_, slice),
            _Cs sliced(this->rownames_, slice), 
            this->colnames_
        );
    };

    template <typename SliceType>
    requires _Cs is_slice_container<SliceType>
        [[nodiscard]] SparseDouble col_sliced(const SliceType& slice) const {
        return SparseDouble(
            this->data_type_, 
            _Cs col_sliced(this->mat_, slice), 
            this->rownames_, 
            _Cs sliced(this->colnames_, slice)
        );
    };

    template <typename OrderType>
    requires _Cs is_order_container<OrderType>
        void row_reorder(const OrderType& order) {
        this->mat_ = _Cs row_reordered(this->mat_, order);
        this->rownames_ = _Cs reordered(this->rownames_, order);
    };

    template <typename QStringContainer>
    requires _Cs is_specific_container<QStringContainer, QString>
        void row_reorder(const QStringContainer& order) {
        this->row_reorder(_Cs index_of(order, this->rownames_));
    }

    template <typename OrderType>
    requires _Cs is_order_container<OrderType>
        void col_reorder(const OrderType& order) {
        this->mat_ = _Cs col_reordered(this->mat_, order);
        this->colnames_ = _Cs reordered(this->colnames_, order);
    };

    template <typename QStringContainer>
    requires _Cs is_specific_container<QStringContainer, QString>
        void col_reorder(const QStringContainer& order) {
        this->col_reorder(_Cs index_of(order, this->colnames_));
    }

    template <typename OrderType, typename OrderType2>
    requires _Cs is_order_container<OrderType>&& _Cs is_order_container<OrderType2>
        void reorder(const OrderType& row_order, const OrderType2& col_order) {
        this->mat_ = _Cs reordered(this->mat_, row_order, col_order);
        this->rownames_ = _Cs reordered(this->rownames_, row_order);
        this->colnames_ = _Cs reordered(this->colnames_, col_order);
    }

    template <typename QStringContainer, typename QStringContainer2>
    requires _Cs is_specific_container<QStringContainer, QString>
        && _Cs is_specific_container<QStringContainer2, QString>
        void reorder(const QStringContainer& row_order, const QStringContainer2& col_order) {
        this->reorder(_Cs index_of(row_order, this->rownames_), _Cs index_of(col_order, this->colnames_));
    }

    template <typename OrderType>
    requires _Cs is_order_container<OrderType>
        [[nodiscard]] SparseDouble reordered(const OrderType& row_order, const OrderType& col_order) const {
        return SparseDouble(
            this->data_type_, 
            _Cs reordered(this->mat_, row_order, col_order), 
            _Cs reordered(this->rownames_, row_order), 
            _Cs reordered(this->colnames_, col_order)
        );
    }

    template <typename QStringContainer, typename QStringContainer2>
    requires _Cs is_specific_container<QStringContainer, QString>
        && _Cs is_specific_container<QStringContainer2, QString>
        [[nodiscard]] SparseDouble reordered(const QStringContainer& row_order, const QStringContainer2& col_order) const {
        return this->reordered(_Cs index_of(row_order, this->rownames_), _Cs index_of(col_order, this->colnames_));
    }

    template <typename OrderType>
    requires _Cs is_order_container<OrderType>
        [[nodiscard]] SparseDouble row_reordered(const OrderType& order) const {
        return SparseDouble(
            this->data_type_, 
            _Cs row_reordered(this->mat_, order), 
            _Cs reordered(this->rownames_, order), 
            this->colnames_
        );
    };

    template <typename QStringContainer>
    requires _Cs is_specific_container<QStringContainer, QString>
        [[nodiscard]] SparseDouble row_reordered(const QStringContainer& order) const {
        return this->row_reordered(_Cs index_of(order, this->rownames_));
    }

    template <typename OrderType>
    requires _Cs is_order_container<OrderType>
        [[nodiscard]] SparseDouble col_reordered(const OrderType& order) const {
        return SparseDouble(
            this->data_type_, 
            _Cs col_reordered(this->mat_, order), 
            this->rownames_, 
            _Cs reordered(this->colnames_, order)
        );
    };

    template <typename QStringContainer>
    requires _Cs is_specific_container<QStringContainer, QString>
        [[nodiscard]] SparseDouble col_reordered(const QStringContainer& order) const {
        return this->col_reordered(_Cs index_of(order, this->colnames_));
    }

	Eigen::ArrayXd get_row(const QString& name) const;
	Eigen::ArrayXd get_column(const QString& name) const;

	Eigen::SparseMatrix<double> get_rows(const QStringList& names) const;

	DenseDouble to_dense() const;

	void clear();
};
