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

class CellchatWorker :
	public QObject
{
	Q_OBJECT
public:

	CellchatWorker(
		const QString& identity,
		const SparseDouble* normalized,
		soap::Species species,
		const QStringList& metadata,
		double minimum_percentage,
		const QString& p_adjust_method,
		int random_state,
		int n_boot,
		int minimum_cell_number,
		const QString& annotation_type)
		:
		identity_(identity),
		origin_(normalized),
		species_(species),
		metadata_(metadata),
		minimum_percentage_(minimum_percentage),
		p_adjust_method_(p_adjust_method),
		random_state_(random_state),
		n_boot_(n_boot),
		minimum_cell_number_(minimum_cell_number),
		annotation_type_(annotation_type)
	{};

private:

	bool load_database();

	void subset_data();

	void subset_database(const QString& type);

	void identify_overexpressed_gene();

	void identify_overexpressed_interactions();

	void project_data();

	void compute_communication_probability();

	void compute_pathway_communication_probability();

	void filter_levels();

	void subset_communication();

	Eigen::ArrayXXd compute_ligand_receptor_expression(const QStringList& gene_ligand_receptor, const DenseDouble& data_use) const;

	Eigen::ArrayXXd compute_complex_expression(const DenseDouble& data_use, const QStringList& complex) const;

	Eigen::ArrayXXd compute_coreceptor_expression(const DenseDouble& data_use, const QString& type) const;

	Eigen::ArrayXd compute_agonist_group_expression(const DenseDouble& data_use, int index, const QStringList& group);

	Eigen::ArrayXd compute_antagonist_group_expression(const DenseDouble& data_use, int index, const QStringList& group);

	QString identity_;

	const SparseDouble* origin_{ nullptr };
	SparseDouble normalized_;

	soap::Species species_;

	DenseDouble projected_;

	QStringList metadata_;
	QStringList levels_;

	QString p_adjust_method_;

	double minimum_percentage_;

	int random_state_;
	int n_boot_;
	int minimum_cell_number_;
	QString annotation_type_;
	int n_levels_;

	std::unique_ptr<CustomMatrix> interactions_;
	std::unique_ptr<CustomMatrix> gene_information_;
	std::unique_ptr<CustomMatrix> complex_;
	std::unique_ptr<CustomMatrix> cofactor_;

	CustomMatrix markers_;
	CustomMatrix significant_ligand_receptor_;
	CustomMatrix interaction_summary_;
	CustomMatrix pathway_summary_;

	Eigen::ArrayX<Eigen::ArrayXXd> probability_;
	Eigen::ArrayX<Eigen::ArrayXXd> p_value_;

	QMap<QString, Eigen::ArrayXXd> pathway_probability_;

	Eigen::ArrayXXi counts_;
	Eigen::ArrayXXd weights_;

	QStringList significant_pathway_;

	std::default_random_engine random_engine_;
	std::uniform_int_distribution<unsigned> unsigned_distribution_;

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_cellchat_ready(CellChat);
};

