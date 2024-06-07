#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"

class Enrichment
{
public:

	enum class DataType : int { 
		Plain = 0, 
		GoEnrichment = 1, 
		KeggEnrichment = 2, 
		MotifEnrichment = 3 
	};

	G_CLASS_FUNCTION_DEFAULT(Enrichment);

	template<typename Mat> requires (!std::same_as<std::decay_t<Mat>, Enrichment>)
		Enrichment(Mat&& df) : mat_(std::forward<Mat>(df)) {}

	DataType data_type_{ DataType::Plain };

	CustomMatrix mat_;

	G_SET_IDENTIFIER("Enrichment");
};
