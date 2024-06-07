#pragma once

#include "Identifier.h"

#include <unordered_map>

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"
#include "GenomeIndex.h"

/*
* modified from signac : Stuart T, Srivastava A, Madad S, Lareau CA, Satija R. Single-cell chromatin state analysis with Signac.
	Nat Methods. 2021 Nov;18(11):1333-1341. doi: 10.1038/s41592-021-01282-5. Epub 2021 Nov 1. Erratum in: Nat Methods. 2022 Feb;19(2):257. PMID: 34725479; PMCID: PMC9255697.
*/

class FragmentsQualityViewWorker
	: public QObject
{
	Q_OBJECT
public:

	FragmentsQualityViewWorker(const SingleCellMultiome* single_cell_multiome, const Fragments* fragments);

	FragmentsQualityViewWorker(const SingleCellAtac* single_cell_atac, const Fragments* fragments);

	enum class WorkMode {AtacFragmentsObject, MultiomeFragmentsObject};

	WorkMode mode_ = WorkMode::AtacFragmentsObject;

	const SingleCellMultiome* single_cell_multiome_{ nullptr };
	const SingleCellAtac* single_cell_atac_{ nullptr };

	const Fragments* fragments_ = nullptr;

	soap::Species species_;

	std::unordered_map<QString, int> cell_index_;

	QVector<std::size_t> length_distribution_;

	QVector<int> mono_nucleosome_count_;
	QVector<int> nucleosome_free_count_;

	GenomicRange genome_;

	QVector<int> center_counts_;
	QVector<int> flank_counts_;
	QVector<int> blacklist_counts_;
	QVector<int> not_blacklist_counts_;

	QVector<int> fragments_in_peak_;
	QVector<int> fragments_not_in_peak_;

	GenomeIndex<IndexMode::Position> flank_index_;
	GenomeIndex<IndexMode::Position> blacklist_index_;
	GenomeIndex<IndexMode::Position> peak_index_;
	GenomeIndex<IndexMode::Position> center_index_;

	bool calculate_metric();

	void get_results();

	bool create_blacklist_index(const QString& bed_file);

public slots:

	void run();

private:

	void create_index();

	void get_tss_position();

	bool get_blacklist_region();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_qc_ready(QMap<QString, QList<double>>);
};

