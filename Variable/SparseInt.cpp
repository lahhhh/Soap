#include "SparseInt.h"

#include "DenseInt.h"
#include "SparseDouble.h"

#include "Custom.h"

DenseInt SparseInt::to_dense() const {

    return { DenseInt::DataType::Plain, this->mat_.toDense(), this->rownames_, this->colnames_ };
};

SparseDouble SparseInt::to_double() const {

    return { SparseDouble::DataType::Plain, this->mat_.cast<double>(), this->rownames_, this->colnames_};
};

qsizetype SparseInt::rows() const {
    return this->mat_.rows();
};

qsizetype SparseInt::cols() const {
    return this->mat_.cols();
};

void SparseInt::clear(){
    this->mat_.resize(0, 0);
    this->rownames_.clear();
    this->colnames_.clear();
};

Eigen::ArrayXi SparseInt::get_row(const QString & name) const{
    int index = this->rownames_.indexOf(name);
    return this->mat_.row(index);
};

Eigen::ArrayXi SparseInt::get_column(const QString & name) const{
    int index = this->colnames_.indexOf(name);
    return this->mat_.col(index);
};

bool SparseInt::is_empty() const{
    return this->rows() == 0 || this->cols() == 0;
}
