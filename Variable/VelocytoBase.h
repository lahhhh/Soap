#pragma once

#include "Identifier.h"

#include "SparseInt.h"

#include "VelocityEstimate.h"
#include "ScveloEstimate.h"

class VelocytoBase
{
public:

    enum class DataType : int { Plain = 0 };

    G_SET_IDENTIFIER("VelocytoBase");

    G_CLASS_FUNCTION_DEFAULT(VelocytoBase);

    DataType data_type_{ DataType::Plain };

    SOAP_SUBMODULES(SparseInt);
    SOAP_SUBMODULES(VelocityEstimate);
    SOAP_SUBMODULES(ScveloEstimate);

    const SparseInt* get_spliced() const;

    const SparseInt* get_unspliced() const;

    template <typename SliceType, typename SliceType2>
        requires _Cs is_slice_container<SliceType>
    && _Cs is_slice_container<SliceType2>
    void slice(const SliceType& row_slice, const SliceType2& col_slice) {

        for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
            data.slice(row_slice, col_slice);
        }

        SUBMODULES(*this, VelocityEstimate).clear();
        SUBMODULES(*this, ScveloEstimate).clear();
    }

    template <typename OrderType>
        requires _Cs is_order_container<OrderType>
    void col_reorder(const OrderType& order) {
    
        for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
            data.col_reorder(order);
        }

        SUBMODULES(*this, VelocityEstimate).clear();
        SUBMODULES(*this, ScveloEstimate).clear();
    }

    template <typename SliceType>
        requires _Cs is_slice_container<SliceType>
    void col_slice(const SliceType& slice) {

        for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
            data.col_slice(slice);
        }

        SUBMODULES(*this, VelocityEstimate).clear();
        SUBMODULES(*this, ScveloEstimate).clear();
    }

    template <typename SliceType>
        requires _Cs is_slice_container<SliceType>
    void row_slice(const SliceType& slice) {

        for (auto& [name, data] : SUBMODULES(*this, SparseInt)) {
            data.row_slice(slice);
        }

        SUBMODULES(*this, VelocityEstimate).clear();
        SUBMODULES(*this, ScveloEstimate).clear();
    }
};

