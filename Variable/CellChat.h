#pragma once

#include "Identifier.h"

#include <QString>
#include <QStringList>
#include <QMap>

#include "CustomMatrix.h"

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
		const QMap<QString, Eigen::ArrayXXd>& pathway_probability,
		const CustomMatrix& interaction_summary,
		const CustomMatrix& pathway_summary
	):
		identity_(identity),
		ligand_receptor_index_(ligand_receptor_index),
		pathway_index_(pathway_index),
		levels_(levels),
		counts_(counts),
		weights_(weights),
		ligand_receptor_probability_(ligand_receptor_probability),
		ligand_receptor_p_value_(ligand_receptor_p_value),
		pathway_probability_(pathway_probability),
		interaction_summary_(interaction_summary),
		pathway_summary_(pathway_summary)
	{}

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

	CustomMatrix interaction_summary_;
	CustomMatrix pathway_summary_;

	G_SET_IDENTIFIER("CellChat");
};

