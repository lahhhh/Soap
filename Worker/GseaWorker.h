#pragma once

#include "Identifier.h"

#include <mutex>

#include "SparseDouble.h"
#include "Gsea.h"

/*
* Modified from GSEA software package :
* Subramanian, Aravind et al. 
* ¡°Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.¡± 
* Proceedings of the National Academy of Sciences of the United States of America vol. 102,43 (2005): 15545-50. doi:10.1073/pnas.0506580102
*/

class GseaWorker 
	: public QObject
{
	Q_OBJECT
public:
	GseaWorker(
		const SparseDouble& sd,
		const QStringList& metadata,
		const QStringList& comparison,
		soap::Species species,
		const QString& database,
		int minimum_overlap,
		int minimum_size,
		int maximum_size,
		int random_state,
		const QString& permutation_type,
		int n_permutation = 1000
	);

	SparseDouble sd_;

	QStringList metadata_;
	QStringList comparison_;
	QStringList gene_list_;
	QStringList pathway_list_;

	soap::Species species_;

	QString database_;
	QString permutation_type_;

	int minimum_overlap_;
	int minimum_size_;
	int maximum_size_;
	int random_state_;
	int n_permutation_;

	QVector<int> sorted_gene_list_;
	QVector<int> gene_locations_;

	std::mutex mutex_;
	std::mutex mutex2_;

	Eigen::ArrayXd original_correlations_;
	Eigen::ArrayXd correlations_;
	Eigen::ArrayXd weights_;

	QVector<QVector<int>> permuted_gene_location_;
	QVector<Eigen::ArrayXd > weights_list_;

	QMap<QString, QStringList> pathway_to_symbol_;
	QMap<QString, int> pathway_to_size_;
	QMap<QString, QVector<int>> pathway_to_location_;

	std::default_random_engine random_engine_;
	std::uniform_int_distribution<unsigned> unsigned_distribution_;

	QMap<QString, QVector<double>> enrichment_scores_;

	QMap<QString, double> enrichment_score_;
	QMap<QString, double> normalized_enrichment_score_;
	QMap<QString, QVector<int> > gene_location_;
	QMap<QString, QVector<double>> point_x_;
	QMap<QString, QVector<double>> point_y_;
	QMap<QString, double> p_values_;
	QMap<QString, double> fdr_;
	QMap<QString, double> fwer_;

	bool load_database();

	void preprocess();

	void get_order();

	void gsea_phenotype();

	void gsea_gene_set();

	void calculate_enrichment_score();

	void permutation_score_phenotype(const QString& geneSetName);

	void permutate_phenotype();

	void calculate_p_value();

	void calculate_fdr();

	void calculate_fwer();

public slots:
	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_gsea_ready(GSEA);

};

