#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"
#include "GenomeIndex.h"

class CalculateCountsByGenomicRangeWorker
	: public QObject
{
	Q_OBJECT
public:

	CalculateCountsByGenomicRangeWorker(const Fragments* fragments, const GenomicRange& genomic_range);

	const Fragments* fragments_ = nullptr;

	GenomicRange genomic_range_;

	GenomeIndex<IndexMode::ID> genome_;

	Eigen::MatrixXi dense_counts_;

	SparseInt* counts_ = nullptr;

	static SparseInt* calculate_counts(const Fragments* fragments, const GenomicRange& genomic_range);

private:

	void create_index();

	bool calculate_counts();

	void build_matrix();

	void find_row(const QString& seq_name, int cell_loc, int start, int end);

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_peak_counts_ready(SparseInt*);	
};

