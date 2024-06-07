#pragma once

#include "Identifier.h"

#include "Custom.h"

#include "SparseInt.h"
#include "SparseDouble.h"
#include "Embedding.h"
#include "DifferentialAnalysis.h"

class DataField
{
public:

    enum class DataType : int { Plain = 0, Rna = 1, Atac = 2, Trans = 3 };

    G_CLASS_FUNCTION_DEFAULT(DataField);

    G_SET_IDENTIFIER("DataField");

    DataType data_type_{ DataType::Plain };

    SOAP_SUBMODULES(SparseInt);
    SOAP_SUBMODULES(SparseDouble);
    SOAP_SUBMODULES(Embedding);
    SOAP_SUBMODULES(DifferentialAnalysis);

    G_QUICK_ACCESS_TYPE(SparseInt, counts, Counts);
    G_QUICK_ACCESS_TYPE(SparseDouble, normalized, Normalized);
    G_QUICK_ACCESS_TYPE(Embedding, pca, Pca);
    G_QUICK_ACCESS_TYPE(Embedding, umap, Umap);
    G_QUICK_ACCESS_TYPE(Embedding, tsne, Tsne);
    G_QUICK_ACCESS_TYPE(Embedding, harmony, Harmony);

    void clear() {
        SUBMODULES(*this, SparseInt).clear();
        SUBMODULES(*this, SparseDouble).clear();
        SUBMODULES(*this, Embedding).clear();
        SUBMODULES(*this, DifferentialAnalysis).clear();
    }

    template <typename SliceType>
    requires _Cs is_slice_container<SliceType>
    void col_slice(const SliceType& slice) {

        for (auto&& [name, data] : SUBMODULES(*this, SparseInt)) {
            data.col_slice(slice);
        }

        for (auto&& [name, data] : SUBMODULES(*this, SparseDouble)) {
            data.col_slice(slice);
        }

        for (auto&& [name, data] : SUBMODULES(*this, Embedding)) {
            data.slice(slice);
        }

        SUBMODULES(*this, DifferentialAnalysis).clear();
    }

    template <typename OrderType>
        requires _Cs is_order_container<OrderType>
    void col_reorder(const OrderType& order) {

        for (auto&& [name, data] : SUBMODULES(*this, SparseInt)) {
            data.col_reorder(order);
        }

        for (auto&& [name, data] : SUBMODULES(*this, SparseDouble)) {
            data.col_reorder(order);
        }

        for (auto&& [name, data] : SUBMODULES(*this, Embedding)) {
            data.reorder(order);
        }

        SUBMODULES(*this, DifferentialAnalysis).clear();
    }
};

