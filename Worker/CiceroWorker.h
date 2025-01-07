#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"

/*
* Modified from Cicero software package :
* Pliner, Hannah A et al. 
* ¡°Cicero Predicts cis-Regulatory DNA Interactions from Single-Cell Chromatin Accessibility Data.¡± 
* Molecular cell vol. 71,5 (2018): 858-871.e8. doi:10.1016/j.molcel.2018.06.044
*/

class CiceroWorker :
	public QObject
{
	Q_OBJECT
public:

	CiceroWorker(
		const SparseInt* atac_counts,
		const Embedding* embedding
	) :
		atac_counts_(atac_counts),
		embedding_(embedding)
	{};

	const SparseInt* atac_counts_{ nullptr };
	const Embedding* embedding_{ nullptr };

	Eigen::MatrixXi cicero_counts_;

	int window_{ 500000 };
	int sample_num_{ 100 };
	int max_it_{ 100 };
	int distance_constraint_{ 250000 };
	int max_sample_windows_{ 500 };
	int max_element_{ 200 };
	int k_{ 50 };

	double distance_parameter_convergence_{ 1e-22 };
	double coaccess_cutoff_{ 0.2 };

	QVector<double> distance_parameters_;

	GenomicRange windows_;

	QStringList peak_chr_names_;
	QVector<int> peak_starts_;
	QVector<int> peak_ends_;

	std::vector<std::pair<QVector<int>, Eigen::MatrixXd>> models;
	Eigen::SparseMatrix<double> connections_;
	QVector<QVector<int>> regulation_groups_;
	Eigen::MatrixXi regulation_group_counts_;

	std::unique_ptr<Cicero> res_{ nullptr };

	bool estimate_distance_parameter();

	bool generate_cicero_models();

	bool assemble_connections();

	bool generate_ccans();

	void create_ccan_matrix();

	double find_ccan_cutoff();

	bool aggregate();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_cicero_ready(Cicero*);
};

