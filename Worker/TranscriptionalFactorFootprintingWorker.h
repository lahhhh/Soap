#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"
#include "TwoBitFileProcessor.h"
#include "Footprint.h"
#include "GenomeIndex.h"

class TranscriptionalFactorFootprintingWorker 
	: public QObject
{
	Q_OBJECT
public:

	TranscriptionalFactorFootprintingWorker(
		const Metadata* metadata,
		std::pair<std::vector<std::string>, std::vector<double>> bias,
		soap::Species species,
		const MotifPosition* motif_position,
		const QStringList& transcriptional_factor_names,
		const Fragments* fragments
	);

	TranscriptionalFactorFootprintingWorker(
		const Metadata* metadata,
		std::pair<std::vector<std::string>, std::vector<double>> bias,
		soap::Species species,
		const GenomicRange& genomic_range,
		const QString& file_name,
		const Fragments* fragments,
		const QString& factor_name,
		const QStringList& levels,
		const QStringList& factors,
		const GraphSettings& graph_settings,
		int height,
		int width
	);

	TranscriptionalFactorFootprintingWorker(
		const Metadata* metadata,
		std::pair<std::vector<std::string>, std::vector<double>> bias,
		soap::Species species,
		const MotifPosition* motif_position,
		const QStringList& transcriptional_factor_names,
		const Fragments* fragments,
		const QString& output_directory,
		const QString& factor_name,
		const QStringList& levels,
		const QStringList& factors,
		const GraphSettings& graph_settings,
		int height,
		int width
	);

	enum class WorkMode{ Object, Batch, GenomicRange};

	WorkMode mode_{ WorkMode::Object };

	const Metadata* metadata_{ nullptr };
	std::pair<std::vector<std::string>, std::vector<double>> bias_;
	soap::Species species_;
	const MotifPosition* motif_position_{ nullptr };

	std::vector<std::string> transcriptional_factor_names_;

	QString output_directory_;
	QString picture_name_;
	QString factor_name_;

	const Fragments* fragments_object_;

	GenomicRange designated_location_;

	QStringList cell_names_;

	QStringList levels_;
	QStringList factors_;
	GraphSettings graph_settings_;
	int height_{ 0 };
	int width_{ 0 };

	TwoBitFileProcessor genome_;

	int upstream_ = 250;

	int downstream_ = 250;

	bool load_genome();

	bool calculate_insertion_bias();

	QVector<double> find_expected_insertions(const std::vector<std::string>& motif_sequence) const;

	std::tuple<bool, std::vector<std::string>, GenomicRange> get_motif_location_sequence(const QString& tf_name);

	std::pair<bool, Eigen::MatrixXi> get_insertion_matrix(
		const GenomeIndex<IndexMode::Strand>& motif_index, int insertion_width) const;

	void create_index();

	void count_position(
		const GenomeIndex<IndexMode::Strand>& motif_index,
		Eigen::MatrixXi& insertion_matrix,
		const QString& sequence_name,
		const int start,
		const int end,
		const int cell_index) const;

	void mode2();

public slots:

	void run();
	 
signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_footprint_ready(Footprint);
	
};

