#pragma once

#include "Identifier.h"

#include <QStringList>

class VelocityEstimate
{
public:

    enum class DataType : int { Plain = 0 };

    G_SET_IDENTIFIER("VelocityEstimate");

    G_CLASS_FUNCTION_DEFAULT(VelocityEstimate);

    DataType data_type_{ DataType::Plain };

    Eigen::MatrixXd current_;
    Eigen::MatrixXd projected_;
    Eigen::MatrixXd deltaE_;

    QStringList gene_names_;
};

