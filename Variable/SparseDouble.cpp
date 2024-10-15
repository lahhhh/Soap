#include "SparseDouble.h"
#include "Custom.h"
#include "DenseDouble.h"

void SparseDouble::clear(){
    this->mat_.resize(0, 0);
    this->rownames_.clear();
    this->colnames_.clear();
};

qsizetype SparseDouble::rows() const {
    return this->mat_.rows();
};

qsizetype SparseDouble::cols() const {
    return this->mat_.cols();
};

Eigen::ArrayXd SparseDouble::get_row(const QString & name) const{
    int index = this->rownames_.indexOf(name);
    return this->mat_.row(index);
};

Eigen::ArrayXd SparseDouble::get_column(const QString & name) const{
    int index = this->colnames_.indexOf(name);
    return this->mat_.col(index);
};

DenseDouble SparseDouble::to_dense() const{

    return DenseDouble(
        DenseDouble::DataType::Plain,
        this->mat_.toDense(), 
        this->rownames_, 
        this->colnames_);
};

Eigen::SparseMatrix<double> SparseDouble::get_rows(const QStringList & names) const{
    auto index = custom::index_of(names, this->rownames_);
    return custom::row_reordered(this->mat_, index);
};
