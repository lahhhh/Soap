#include "DenseInt.h"

#include "Custom.h"
#include "DenseDouble.h"

DenseDouble DenseInt::to_double() const {

    return DenseDouble{ DenseDouble::DataType::Plain, this->mat_.cast<double>(), this->rownames_, this->colnames_ };
};

qsizetype DenseInt::rows() const {
    return this->mat_.rows();
};

qsizetype DenseInt::cols() const {
    return this->mat_.cols();
};

void DenseInt::clear() {
    mat_.resize(0, 0);
    rownames_.clear();
    colnames_.clear();
};

Eigen::ArrayXi DenseInt::get_row(const QString& name) const {
    qsizetype index = rownames_.indexOf(name);
    return this->mat_.row(index);
};

Eigen::MatrixXi DenseInt::get_rows(const QStringList& names) const {
    return this->mat_(_Cs index_of(names, this->rownames_), Eigen::all);
};

Eigen::ArrayXXi DenseInt::get_rows_array(const QStringList& names) const {
    return this->get_rows(names).array();
};

Eigen::ArrayXi DenseInt::get_column(const QString& name) const {
    qsizetype index = this->colnames_.indexOf(name);
    return this->mat_.col(index);
};