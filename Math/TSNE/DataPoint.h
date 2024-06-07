#pragma once

#include "Identifier.h"
#include <vector>

class DataPoint
{

public:

    DataPoint() = default;
    
    DataPoint(const DataPoint& other) {

        if (this != &other) {

            this->dimension_ = other.dimension_;
            this->index_ = other.index_;
            this->data_ = (double*)malloc(sizeof(double) * this->dimension_);

            std::memcpy(this->data_, other.data_, sizeof(double) * this->dimension_);
        }
    };
    
    DataPoint(DataPoint&& other) noexcept {

        this->dimension_ = other.dimension_;
        this->index_ = other.index_;
        this->data_ = other.data_;

        other.data_ = nullptr;
        other.dimension_ = 0;
        other.index_ = -1;
    }
    
    DataPoint& operator=(const DataPoint& other) {
        if (this != &other) {

            if (this->data_ != nullptr){
                free(this->data_);
            }

            this->dimension_ = other.dimension_;
            this->index_ = other.index_;
            this->data_ = (double*)malloc(sizeof(double) * this->dimension_);

            std::memcpy(this->data_, other.data_, sizeof(double) * this->dimension_);
        }

        return *this;
    }
    
    DataPoint& operator=(DataPoint&& other) noexcept {
    
        if (this->data_ != nullptr) {
            free(this->data_);
        }

        this->dimension_ = other.dimension_;
        this->index_ = other.index_;
        this->data_ = other.data_;

        other.data_ = nullptr;
        other.dimension_ = 0;
        other.index_ = -1;

        return *this;
    };

    ~DataPoint() {
        if (this->data_ != nullptr) {
            free(this->data_);
        }
    }

    double* data_{ nullptr };
    int dimension_{ 0 };
    int index_{ -1 };


    DataPoint(int dimension, int index, const Eigen::MatrixXd& data):
        dimension_(dimension),
        index_(index)
    {
        this->data_ = (double*)malloc(sizeof(double) * dimension);

        std::memcpy(this->data_, data.col(index).data(), sizeof(double) * dimension);
    }

    int index() const { return this->index_; }
    int dimensionality() const { return this->dimension_; }
    double x(int d) const { return this->data_[d]; }
};

