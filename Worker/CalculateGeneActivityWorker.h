#pragma once

#include "Identifier.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

#include "GenomeIndex.h"

class CalculateGeneActivityWorker
	: public QObject
{
	Q_OBJECT
public:

	CalculateGeneActivityWorker(
		soap::Species species,
		const Fragments* fragments,
		bool use_tss
	)
		:
		species_(species),
		fragments_(fragments),
		use_tss_(use_tss)
	{};

	soap::Species species_;

	const Fragments* fragments_{ nullptr };

	bool use_tss_{ true };

	int tss_width{ 1000 };

	GenomeIndex<IndexMode::ID> genome_;

	QStringList gene_names_;

	Eigen::MatrixXi dense_counts_;

	SparseInt* counts_ = nullptr;

	bool create_index();

	bool calculate_counts();

	int create_index_human();

	void build_matrix();

	void find_row(const QString& seq_name, int cell_loc, int start, int end);

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_gene_activity_ready(SparseInt*);	

};

