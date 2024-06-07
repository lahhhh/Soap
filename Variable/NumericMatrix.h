#pragma once

#include "Identifier.h"

#include <QStringList>

#include "Custom.h"

class NumericMatrix
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(NumericMatrix);

	NumericMatrix(const Eigen::MatrixXd& mat) :	data_(mat) {}

	DataType data_type_{ DataType::Plain };

	Eigen::MatrixXd data_;

	G_SET_IDENTIFIER("NumericMatrix");
};



