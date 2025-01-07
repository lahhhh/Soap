#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"
#include "TwoBitFileProcessor.h"

/*
* Modified from ChromVAR software package :
* Schep AN, Wu B, Buenrostro JD and Greenleaf WJ (2017). 
* ¡°chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data.¡± 
* Nature Methods. doi: 10.1038/nmeth.4401.
*/

class ChromVARWorker :
	public QObject
{
	Q_OBJECT
public:

	ChromVARWorker(
		soap::Species species,
		const MotifPosition* motif_position,
		const SparseInt* atac_counts
	) :
		species_(species),
		motif_position_(motif_position),
		atac_counts_(atac_counts)
	{};

	soap::Species species_;
	const MotifPosition* motif_position_{ nullptr };
	const SparseInt* atac_counts_{ nullptr };

	TwoBitFileProcessor genome_;

	std::vector<std::string> peak_sequences_;

	Eigen::ArrayXd gc_;
	Eigen::ArrayXd expectations_;

	Eigen::MatrixXi background_peaks_;

	std::unique_ptr<ChromVAR> res_{ nullptr };

	bool validate_input();

	bool add_gc_bias();

	bool get_background_peaks();

	bool compute_deviations();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_chromvar_ready(ChromVAR*);
};

