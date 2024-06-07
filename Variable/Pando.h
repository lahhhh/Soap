#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"
#include "GeneName.h"

class Pando
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(Pando);

	template<typename Mat> requires (!std::same_as<std::decay_t<Mat>, Pando>)
		Pando(Mat&& df) : mat_(std::forward<Mat>(df)) {}

	DataType data_type_{ DataType::Plain };

	CustomMatrix mat_;

	SOAP_SUBMODULES(GeneName);

	G_SET_IDENTIFIER("Pando");
};

