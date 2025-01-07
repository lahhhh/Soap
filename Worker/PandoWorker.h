#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"
#include "GenomeIndex.h"

/*
* modified from r package Pando
* Fleck, Jonas Simon et al. ¡°Inferring and perturbing cell fate regulomes in human brain organoids.¡± 
* Nature vol. 621,7978 (2023): 365-372. doi:10.1038/s41586-022-05279-8
*/

class PandoWorker
	: public QObject
{
	Q_OBJECT
public:

	PandoWorker(
		const SingleCellMultiome* single_cell_multiome,
		const QString& motif_database,
		bool filter_cell,
		bool use_tss,
		const Eigen::ArrayX<bool> cell_filter,
		int n_feature,
		double peak_correlation_threshold,
		double motif_correlation_threshold,
		const QString& p_adjust_method
	):
		single_cell_multiome_(single_cell_multiome),
		motif_database_(motif_database),
		filter_cell_(filter_cell),
		use_tss_(use_tss),
		cell_filter_(cell_filter),
		n_feature_(n_feature),
		peak_correlation_threshold_(peak_correlation_threshold),
		motif_correlation_threshold_(motif_correlation_threshold),
		p_adjust_method_(p_adjust_method)
	{}

	const SingleCellMultiome* single_cell_multiome_{ nullptr };
	QString motif_database_;

	bool filter_cell_{ true };
	bool use_tss_{ true };

	Eigen::ArrayX<bool> cell_filter_;

	int n_feature_{ 0 };
	double peak_correlation_threshold_{ 0.0 };
	double motif_correlation_threshold_{ 0.0 };

	int upstream_{ 100000 };
	int downstream_{ 0 };

	QString p_adjust_method_;

	GenomicRange annotation_;

	GenomicRange peaks_;
	GenomeIndex<IndexMode::ID> peak_index_;

	MotifPosition motif_position_;

	QStringList gene_use_;
	QList<QList<int>> gene_to_peak_;

	Pando res_;

	bool initiate_grn();

	bool infer_grn();

	QStringList find_variable_features(int n_feature);

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_pando_ready(Pando);
};

