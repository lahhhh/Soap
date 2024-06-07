#pragma once

#include "Identifier.h"

#include <QStringList>
#include <QMap>

class CNV
{
public:

	enum class DataType : int { Plain = 0, InferCnv = 1, SciCnv = 2 };
		
	G_CLASS_FUNCTION_DEFAULT(CNV);

	template<
		typename MatrixType,
		typename InfoType,
		typename InfoType2
	>
	CNV(
		MatrixType&& mat,
		InfoType&& cluster_info,
		InfoType2&& chromosome_info,
		DataType data_type = DataType::SciCnv
	) :
		mat_(std::forward<MatrixType>(mat)),
		cluster_info_(std::forward<InfoType>(cluster_info)),
		chromosome_info_(std::forward<InfoType2>(chromosome_info)),
		data_type_(data_type)
	{}

	template<
		typename MatrixType,
		typename InfoType
	>
	CNV(
		MatrixType&& mat,
		InfoType&& chromosome_info,
		DataType data_type = DataType::InferCnv
	) :
		mat_(std::forward<MatrixType>(mat)),
		chromosome_info_(std::forward<InfoType>(chromosome_info)),
		data_type_(data_type)
	{}

	DataType data_type_{ DataType::Plain };

	Eigen::MatrixXd mat_;

	std::vector<std::tuple<QString, int, int> > cluster_info_;
	std::vector<std::tuple<QString, int, int> > chromosome_info_;	

	G_SET_IDENTIFIER("CNV");
};
