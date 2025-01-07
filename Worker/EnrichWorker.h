#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"
#include "MotifPosition.h"

class EnrichWorker :
	public QObject
{
	Q_OBJECT
public:
	EnrichWorker(
		const QString& database,
		const QString& enrich_type,
		const QStringList& gene_names,
		const QString& ontology,
		const soap::Species& species,
		const QString& p_adjust_method,
		double p_threshold
	);

	EnrichWorker(
		const QStringList& peak_names,
		const MotifPosition* motif_position,
		const QString& enrich_type,
		const QString& p_adjust_method,
		double p_threshold
	);

	enum class WorkMode { EnrichPathway, EnrichMotif };

	WorkMode mode_{ WorkMode::EnrichPathway };

	QStringList gene_names_;
	QStringList peak_names_;

	const MotifPosition* motif_position_{ nullptr };

	QString database_;
	QString enrich_type_;
	QString ontology_;
	QString p_adjust_method_;

	soap::Species species_;

	double p_threshold_;

	CustomMatrix res_;

	QString res_name_;

	void enrich_pathway();

	void enrich_motif();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_enrichment_ready(CustomMatrix, QString);

};

