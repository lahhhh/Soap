#pragma once

#include "Identifier.h"

#include "Custom.h"

#include <QStringList>

#include "DenseDouble.h"

class Embedding
{
public:

    enum class DataType : int {
        Plain = 0,
        Pca = 1,
        Umap = 2,
        Tsne = 3,
        Harmony = 4
    };

    G_CLASS_FUNCTION_DEFAULT(Embedding);

    template <
        typename MatrixType,
        typename StringList,
        typename StringList2
    >
    Embedding(
        DataType data_type,
        MatrixType&& mat, 
        StringList&& rownames, 
        StringList2&& embedding_names
    ):
        data_type_(data_type),
        data_(
            DenseDouble::DataType::Plain,
            std::forward<MatrixType>(mat),
            std::forward<StringList>(rownames),
            std::forward<StringList2>(embedding_names)
        )
    { };    

    DataType data_type_{ DataType::Plain };

    DenseDouble data_;

    template <typename SliceType>
        requires _Cs is_slice_container<SliceType>
    void slice(const SliceType& slice) {

        this->data_.row_slice(slice);
    }

    template <typename OrderType>
        requires _Cs is_order_container<OrderType>
    void reorder(const OrderType& order) {

        this->data_.row_reorder(order);
    }


    G_SET_IDENTIFIER("Embedding");
};

