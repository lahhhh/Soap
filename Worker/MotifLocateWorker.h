#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"
#include "TwoBitFileProcessor.h"
#include "PatternWeightMatrix.h"

// modified from motifmatchr : https://github.com/GreenleafLab/motifmatchr

class MotifLocateWorker : public QObject
{
	Q_OBJECT
public:

	MotifLocateWorker(
		const QStringList& peak_names,
		const QString& database_name,
		soap::Species species) :
		peak_names_(peak_names),
		database_name_(database_name),
		species_(species) 
	{};

	QStringList peak_names_;
	QString database_name_;
	soap::Species species_;

	TwoBitFileProcessor genome_file_;

	GenomicRange peak_;

	QStringList peak_sequences_;

	Eigen::Array4d sequence_background_;

	std::map<QString, PatternWeightMatrix> pattern_database_;

	MotifPosition mp_;

	void get_peak();

	bool get_peak_sequence();

	void get_base_background();

	void get_result();


public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_motif_location_ready(MotifPosition);	
};

