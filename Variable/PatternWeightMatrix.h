#pragma once

#include "Identifier.h"

#include "DenseDouble.h"
// modified from TFBSTools : https://github.com/ge11232002/TFBSTools

class PatternWeightMatrix
{

public:
	G_CLASS_FUNCTION_DEFAULT(PatternWeightMatrix);

	enum class DataType { Plain = 0, Frequency = 1, Log2 = 2, Log = 3 };
	
	DataType data_type_{ DataType::Plain };

	DenseDouble weight_;

	QString motif_name_;

	QString transcriptional_factor_name_;

	void from_frequency_to_log2();

	void convert(const Eigen::Array4d& background_frequency);

	Eigen::Index size() const;
};

