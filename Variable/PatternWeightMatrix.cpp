#include "PatternWeightMatrix.h"

#include <QFile>
#include <QTextStream>
#include "Custom.h"
#include <fstream>

void PatternWeightMatrix::convert(const Eigen::Array4d& background_frequency) {
	if (this->data_type_ == DataType::Log2) {
		for (int i = 0; i < 4; ++i) {
			this->weight_.mat_.row(i).array() -= log2(0.25) - log2(background_frequency[i]);
		}
	}
	else if (this->data_type_ == DataType::Log) {
		for (int i = 0; i < 4; ++i) {
			this->weight_.mat_.row(i).array() -= log(0.25) - log(background_frequency[i]);
		}
	}
};

void PatternWeightMatrix::from_frequency_to_log2() {

	if (this->data_type_ == DataType::Log2) {
		return;
	}

	double pseudo_count = 0.8;
	double sequence_count = this->weight_.mat_.col(0).sum() + pseudo_count;
	double sub_pseudo = pseudo_count / 4;
	this->weight_.mat_ = log2((this->weight_.mat_.array() + sub_pseudo) / sequence_count * 4).matrix();
	this->data_type_ = DataType::Log2;
};

Eigen::Index PatternWeightMatrix::size() const {
	return this->weight_.mat_.cols();
};