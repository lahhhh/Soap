#pragma once

#include "Identifier.h"

#include "CustomTemplates.h"

#include <vector>
#include <QMap>
#include <QString>

class Fragments
{
public:

    enum class DataType : int { Plain = 0 };

    G_CLASS_FUNCTION_DEFAULT(Fragments);

    G_SET_IDENTIFIER("Fragments");

    DataType data_type_{ DataType::Plain };

    std::map<QString, std::vector<std::pair<std::vector<int>, std::vector<int>>>> data_;

    QStringList cell_names_;

    auto begin() noexcept {

        return this->data_.begin();
    }

    auto end() noexcept {

        return this->data_.end();
    }
    
    auto begin() const noexcept {

        return this->data_.begin();
    }

    auto end() const noexcept {

        return this->data_.end();
    }

    const std::pair<std::vector<int>, std::vector<int>>& get_fragments(
        const QString& sequence_name,
        const QString& cell_name
    ) const;

    void finalize();

    void adjust_length_by_cell_name_length();

    template <typename SliceType>
        requires custom::is_slice_container<SliceType>
    void slice(const SliceType& slice) {

        this->cell_names_ = custom::sliced(this->cell_names_, slice);

        for (auto& [chr, data] : this->data_) {
            data = custom::sliced(data, slice);
        }
    }

    template <typename SliceType>
        requires custom::is_slice_container<SliceType>
    [[nodiscard]] Fragments sliced(const SliceType& slice) const {

        Fragments ret;
        
        ret.cell_names_ = custom::sliced(this->cell_names_, slice);

        for (const auto& [chr, data] : this->data_) {
            ret.data_[chr] = custom::sliced(data, slice);
        }
        return ret;
    }

    template <typename OrderType>
        requires custom::is_order_container<OrderType>
    void reorder(const OrderType& order) {

        this->cell_names_ = custom::reordered(this->cell_names_, order);

        for (auto& [chr, data] : this->data_) {
            data = custom::reordered(data, order);
        }
    }

    template <typename OrderType>
        requires custom::is_order_container<OrderType>
    [[nodiscard]] Fragments reordered(const OrderType& order) const {

        Fragments ret;
        
        ret.cell_names_ = custom::reordered(this->cell_names_, order);

        for (const auto& [chr, data] : this->data_) {
            ret.data_[chr] = custom::reordered(data, order);
        }
        return ret;
    }
};

