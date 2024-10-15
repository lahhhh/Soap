#include "DenseDouble.h"
#include "Custom.h"

#include "SparseDouble.h"

SparseDouble DenseDouble::to_sparse() const {
    
    return SparseDouble{ SparseDouble::DataType::Plain, this->mat_.sparseView(), this->rownames_, this->colnames_ };
};

qsizetype DenseDouble::rows() const {
    return this->mat_.rows();
};

qsizetype DenseDouble::cols() const {
    return this->mat_.cols();
};

void DenseDouble::clear(){
    mat_.resize(0, 0);
    rownames_.clear();
    colnames_.clear();
};

Eigen::ArrayXd DenseDouble::get_row(const QString& name) const {
    qsizetype index = this->rownames_.indexOf(name);
    return this->mat_.row(index);
};

Eigen::MatrixXd DenseDouble::get_rows(const QStringList& names) const {
    return this->mat_(custom::index_of(names, this->rownames_), Eigen::all);
};

Eigen::ArrayXXd DenseDouble::get_rows_array(const QStringList& names) const {
    return this->get_rows(names).array();
};

Eigen::ArrayXd DenseDouble::get_column(const QString& name) const {
    qsizetype index = this->colnames_.indexOf(name);
    return this->mat_.col(index);
};