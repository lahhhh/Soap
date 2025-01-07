#pragma once

#include "Identifier.h"

#include <unordered_map>

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"
#include "GenomeIndex.h"

class TssPlotWorker :
	public QObject
{
	Q_OBJECT
public:

	TssPlotWorker(soap::Species species, const Fragments* fragments);

	const Fragments* fragments_{nullptr};

	soap::Species species_;

	GenomicRange genome_;

	GenomeIndex<IndexMode::IdStrand> tss_index_;

	Eigen::ArrayXXd tss_matrix_;

	QVector<double> tss_vector_;

	void calculate_tss();

	void get_tss_position();

	void find_insertion_location_in_tss(
		const QString& seq_name,
		int start,
		int end
	);

	void get_results();

	void create_index();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_tss_ready(Eigen::ArrayXXd tss_matrix, QVector<double> tss_vector);
	
};

