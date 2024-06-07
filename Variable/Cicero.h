#pragma once

#include "Identifier.h"

#include "DenseInt.h"
#include "DenseDouble.h"
#include "Embedding.h"
#include "DifferentialAnalysis.h"

class Cicero
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(Cicero);

	DataType data_type_{ DataType::Plain };

	G_SET_IDENTIFIER("Cicero");

	Eigen::SparseMatrix<double> connections_;
	QVector<QVector<int>> regulation_groups_;
	QStringList peak_names_;
	DenseInt regulation_group_counts_;
	DenseDouble regulation_group_normalized_;

	SOAP_SUBMODULES(Embedding);
	SOAP_SUBMODULES(DifferentialAnalysis);
};
