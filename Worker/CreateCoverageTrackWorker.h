#pragma once

#include "Identifier.h"

#include "SignalEmitter.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class CreateCoverageTrackWorker :
	public QObject
{
	Q_OBJECT
public:

	CreateCoverageTrackWorker(
		const Metadata* metadata,
		const Fragments* fragments,
		const QString& group_name
	) :
		metadata_(metadata),
		fragments_(fragments),
		group_name_(group_name) 
	{};

	const Metadata* metadata_{ nullptr };

	const Fragments* fragments_{ nullptr };

	QString group_name_;

	QString coverage_track_name_;

	QMap<QString, std::vector< std::vector< std::pair<int, float> > > > temp_data_;

	std::vector<int> cell_index_;

	std::vector<double> cell_size_;

	QStringList group_factors_;
	QMap<QString, int> group_distribution_;

	std::unique_ptr<CoverageTrack> res_{ nullptr };

	void build_index();

	bool calculate_fragments_size();

	bool calculate_matrix();

	bool load_annotation();

	void create_track();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_coverage_track_ready(CoverageTrack*);

};

