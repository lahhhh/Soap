#pragma once

#include "Identifier.h"

#include <QStringList>

#include "DifferentialAnalysis.h"
#include "Embedding.h"

class ChromVAR
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(ChromVAR);

	G_SET_IDENTIFIER("ChromVAR");

	DataType data_type_{ DataType::Plain };
	QStringList motif_names_;
	Eigen::MatrixXd z_;
	Eigen::MatrixXd dev_;

	SOAP_SUBMODULES(DifferentialAnalysis);
	SOAP_SUBMODULES(Embedding);
};

