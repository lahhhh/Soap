#pragma once

#include "Identifier.h"

#include <QString>
#include <QStringList>
#include <QMap>

// deprecated

class CellChat
{
public:
		
	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(CellChat);

	CellChat(
		const QString& identity, 
		const QStringList& ligand_receptor_index, 
		const QStringList& pathway_index, 
		const QStringList& levels, 
		const Eigen::ArrayXXi& counts, 
		const Eigen::ArrayXXd& weights,
		const Eigen::ArrayX<Eigen::ArrayXXd>& ligand_receptor_probability, 
		const Eigen::ArrayX<Eigen::ArrayXXd>& ligand_receptor_p_value, 
		const QMap<QString, Eigen::ArrayXXd>& pathway_probability
	);

	DataType data_type_{DataType::Plain };
	
	QString identity_;

	QStringList ligand_receptor_index_;
	QStringList pathway_index_;
	QStringList levels_;

	Eigen::ArrayXXi counts_;
	Eigen::ArrayXXd weights_;

	Eigen::ArrayX<Eigen::ArrayXXd> ligand_receptor_probability_;
	Eigen::ArrayX<Eigen::ArrayXXd> ligand_receptor_p_value_;

	QMap<QString, Eigen::ArrayXXd> pathway_probability_;

	G_SET_IDENTIFIER("CellChat");
};

