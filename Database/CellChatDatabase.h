#pragma once

#include "Identifier.h"

#include <random>

#include "CustomMatrix.h"
#include "SparseInt.h"
#include "SparseDouble.h"
#include "DenseDouble.h"
#include "CellChat.h"

/*
* Modified from CellChat software package :
* Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan CH, Myung P, Plikus MV, Nie Q.
* Inference and analysis of cell-cell communication using CellChat.
* Nat Commun. 2021 Feb 17;12(1):1088. doi: 10.1038/s41467-021-21246-9. PMID: 33597522; PMCID: PMC7889871.
*/

// deprecated

class CellChatDatabase
{
public:
	CellChatDatabase() = default;

	CellChatDatabase(const CellChatDatabase&) = delete;

	CellChat cellchat(
		const QString& identity,
		SparseDouble* normalized,
		soap::Species species,
		const QStringList& metadata,
		double minimum_percentage,
		const QString& p_adjust_method,
		int random_state,
		int n_boot,
		int minimum_cell_number,
		const QString& annotation_type
	);

private:

	void load_human_database();

	void load_mouse_database();

	void subset_data();

	void subset_database(const QString& type);

	void identify_overexpressed_gene();

	void identify_overexpressed_interactions();

	void project_data();

	void compute_communication_probability();

	void compute_pathway_communication_probability();

	void filter_levels(SparseDouble* normalized, const QStringList& metadata, int minimum_cell_number);

	void subset_communication();

	void clear();

	Eigen::ArrayXXd compute_ligand_receptor_expression(const QStringList& gene_ligand_receptor, const DenseDouble& data_use);

	Eigen::ArrayXXd compute_complex_expression(const DenseDouble& data_use, const QStringList& complex);

	Eigen::ArrayXXd compute_coreceptor_expression(const DenseDouble& data_use, const QString& type);

	Eigen::ArrayXd compute_agonist_group_expression(const DenseDouble& data_use, int index, const QStringList& group);

	Eigen::ArrayXd compute_antagonist_group_expression(const DenseDouble& data_use, int index, const QStringList& group);

	SparseDouble* normalized_{ nullptr };

	DenseDouble projected_;

	QStringList metadata_;
	QStringList levels_;

	QString p_adjust_method_;

	double minimum_percentage_;

	int random_state_;
	int n_boot_;
	int n_levels_;

	QMap<QString, CustomMatrix> interaction_database_;
	QMap<QString, CustomMatrix> gene_information_database_;
	QMap<QString, CustomMatrix> complex_database_;
	QMap<QString, CustomMatrix> cofactor_database_;
	QMap<QString, SparseInt> protein_protein_interaction_database_;

	CustomMatrix* interactions_ = nullptr;
	CustomMatrix* gene_information_ = nullptr;
	CustomMatrix* complex_ = nullptr;
	CustomMatrix* cofactor_ = nullptr;

	CustomMatrix markers_;
	CustomMatrix significant_ligand_receptor_;

	SparseInt* protein_protein_interaction_ = nullptr;

	Eigen::ArrayX<Eigen::ArrayXXd> probability_;
	Eigen::ArrayX<Eigen::ArrayXXd> p_value_;

	QMap<QString, Eigen::ArrayXXd> pathway_probability_;

	Eigen::ArrayXXi counts_;
	Eigen::ArrayXXd weights_;

	QStringList significant_pathway_;

	std::default_random_engine random_engine_;
	std::uniform_int_distribution<unsigned> unsigned_distribution_;
};

