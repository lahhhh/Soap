#pragma once

#include "Identifier.h"

#include <QStringList>

#include "Custom.h"

class DenseDouble;

class DenseInt
{
public:

    enum class DataType : int { Plain = 0, Counts = 1 };

    G_CLASS_FUNCTION_DEFAULT(DenseInt);

    template <typename DenseMatrix, typename String1, typename String2>
    DenseInt(
        DataType data_type,
        DenseMatrix&& mat,
        String1&& rownames,
        String2&& colnames
    ) :
        data_type_(data_type),
        mat_(std::forward<DenseMatrix>(mat)),
        rownames_(std::forward<String1>(rownames)),
        colnames_(std::forward<String2>(colnames))
    { }

    DataType data_type_{ DataType::Plain };

    Eigen::MatrixXi mat_;
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
    [[nodiscard]] DenseInt sliced(const SliceType& row_slice, const SliceType2& col_slice) const {
        return DenseInt(
            this->data_type_,
            _Cs sliced(this->mat_, row_slice, col_slice),
            _Cs sliced(this->rownames_, row_slice),
            _Cs sliced(this->colnames_, col_slice)
        );
    }

    template <typename SliceType>
        requires _Cs is_slice_container<SliceType>
    [[nodiscard]] DenseInt row_sliced(const SliceType& slice) const {
        return DenseInt(
            this->data_type_,
            _Cs row_sliced(this->mat_, slice),
            _Cs sliced(this->rownames_, slice),
            this->colnames_
        );
    };

    template <typename SliceType>
        requires _Cs is_slice_container<SliceType>
    [[nodiscard]] DenseInt col_sliced(const SliceType& slice) const {
        return DenseInt(
            this->data_type_,
            _Cs col_sliced(this->mat_, slice),
            this->rownames_,
            _Cs sliced(this->colnames_, slice)
        );
    };

    template <typename OrderType>
        requires _Cs is_order_container<OrderType>
    void row_reorder(const OrderType& order) {
        this->mat_ = this->mat_(order, Eigen::all).eval();
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
        this->mat_ = this->mat_(Eigen::all, order).eval();
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
        this->mat_ = this->mat_(row_order, col_order).eval();
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
    [[nodiscard]] DenseInt reordered(const OrderType& row_order, const OrderType& col_order) const {
        return DenseInt(
            this->data_type_,
            this->mat_(row_order, col_order),
            _Cs reordered(this->rownames_, row_order),
            _Cs reordered(this->colnames_, col_order)
        );
    }

    template <typename QStringContainer, typename QStringContainer2>
        requires _Cs is_specific_container<QStringContainer, QString>
    && _Cs is_specific_container<QStringContainer2, QString>
        [[nodiscard]] DenseInt reordered(const QStringContainer& row_order, const QStringContainer2& col_order) const {

        return this->reordered(_Cs index_of(row_order, this->rownames_), _Cs index_of(col_order, this->colnames_));
    }

    template <typename OrderType>
        requires _Cs is_order_container<OrderType>
    [[nodiscard]] DenseInt row_reordered(const OrderType& order) const {
        return DenseInt(
            this->data_type_,
            this->mat_(order, Eigen::all),
            _Cs reordered(this->rownames_, order),
            this->colnames_
        );
    };

    template <typename QStringContainer>
        requires _Cs is_specific_container<QStringContainer, QString>
    [[nodiscard]] DenseInt row_reordered(const QStringContainer& order) const {
        return this->row_reordered(_Cs index_of(order, this->rownames_));
    }

    template <typename OrderType>
        requires _Cs is_order_container<OrderType>
    [[nodiscard]] DenseInt col_reordered(const OrderType& order) const {
        return DenseInt(
            this->data_type_,
            this->mat_(Eigen::all, order),
            this->rownames_,
            _Cs reordered(this->colnames_, order)
        );
    };

    template <typename QStringContainer>
        requires _Cs is_specific_container<QStringContainer, QString>
    [[nodiscard]] DenseInt col_reordered(const QStringContainer& order) const {

        return this->col_reordered(_Cs index_of(order, this->colnames_));
    }

    Eigen::ArrayXi get_row(const QString& name) const;
    Eigen::MatrixXi get_rows(const QStringList& names) const;
    Eigen::ArrayXXi get_rows_array(const QStringList& names) const;
    Eigen::ArrayXi get_column(const QString& name) const;

    void clear();

    DenseDouble to_double() const;

    G_SET_IDENTIFIER("DenseInt");
};
