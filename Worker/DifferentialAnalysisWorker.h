#pragma once

#include "Identifier.h"

#include "BulkRna.h"
#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class DifferentialAnalysisWorker : public QObject
{
	Q_OBJECT
public:
	DifferentialAnalysisWorker(
		const SingleCellRna* single_cell_rna,
		const QString& metadata_name,
		const QStringList& metadata,
		const QStringList& comparison,
		double minimum_percentage,
		const QString& p_adjust_method
	) :
		single_cell_rna_(single_cell_rna),
		metadata_name_(metadata_name),
		metadata_(metadata),
		comparison_(comparison),
		feature_names_(single_cell_rna->counts()->rownames_),
		minimum_percentage_(minimum_percentage),
		p_adjust_method_(p_adjust_method),
		mode_(WorkMode::SingleCellRna)
	{}

	DifferentialAnalysisWorker(
		const BulkRna* bulk_rna,
		const QString& metadata_name,
		const QStringList& metadata,
		const QStringList& comparison,
		double minimum_percentage,
		const QString& p_adjust_method
	) :
		bulk_rna_(bulk_rna),
		metadata_name_(metadata_name),
		metadata_(metadata),
		comparison_(comparison),
		feature_names_(bulk_rna->counts()->rownames_),
		minimum_percentage_(minimum_percentage),
		p_adjust_method_(p_adjust_method),
		mode_(WorkMode::BulkRna)
	{}

	DifferentialAnalysisWorker(
		const SingleCellAtac* single_cell_atac,
		const QString& metadata_name,
		const QStringList& metadata,
		const QStringList& comparison,
		double minimum_percentage,
		const QString& p_adjust_method
	) :
		single_cell_atac_(single_cell_atac),
		metadata_name_(metadata_name),
		metadata_(metadata),
		comparison_(comparison),
		feature_names_(single_cell_atac->counts()->rownames_),
		minimum_percentage_(minimum_percentage),
		p_adjust_method_(p_adjust_method),
		mode_(WorkMode::SingleCellAtac)
	{}

	DifferentialAnalysisWorker(
		const DataField* data_field,
		const QString& metadata_name,
		const QStringList& metadata,
		const QStringList& comparison,
		double minimum_percentage,
		const QString& p_adjust_method
	) :
		data_field_(data_field),
		metadata_name_(metadata_name),
		metadata_(metadata),
		comparison_(comparison),
		feature_names_(data_field->counts()->rownames_),
		minimum_percentage_(minimum_percentage),
		p_adjust_method_(p_adjust_method),
		mode_(WorkMode::SingleCellMultiome)
	{}

	DifferentialAnalysisWorker(
		const ChromVAR* chrom_var,
		const QString& metadata_name,
		const QStringList& metadata,
		const QStringList& comparison,
		const QString& p_adjust_method
	) :
		chrom_var_(chrom_var),
		metadata_name_(metadata_name),
		metadata_(metadata),
		comparison_(comparison),
		feature_names_(chrom_var->motif_names_),
		p_adjust_method_(p_adjust_method),
		mode_(WorkMode::ChromVAR)
	{}

	DifferentialAnalysisWorker(
		const Cicero* cicero,
		const QString& metadata_name,
		const QStringList& metadata,
		const QStringList& comparison,
		const QString& p_adjust_method,
		double minimum_percentage
	) :
		cicero_(cicero),
		metadata_name_(metadata_name),
		metadata_(metadata),
		comparison_(comparison),
		feature_names_(cicero->regulation_group_counts_.rownames_),
		p_adjust_method_(p_adjust_method),
		minimum_percentage_(minimum_percentage),
		mode_(WorkMode::Cicero)
	{}

	const BulkRna* bulk_rna_{ nullptr };
	const SingleCellRna* single_cell_rna_{ nullptr };
	const SingleCellAtac* single_cell_atac_{ nullptr };
	const DataField* data_field_{ nullptr };
	const ChromVAR* chrom_var_{ nullptr };
	const Cicero* cicero_{ nullptr };

	Eigen::SparseMatrix<double> sparse_data_;
	Eigen::MatrixXd dense_data_;

	QString metadata_name_;

	QStringList metadata_;
	QStringList comparison_;
	QStringList feature_names_;

	double minimum_percentage_{ 0.0 };

	QString p_adjust_method_;

	enum class WorkMode{BulkRna, SingleCellRna, SingleCellAtac, SingleCellMultiome, ChromVAR, Cicero};

	WorkMode mode_ = WorkMode::SingleCellRna;

	void sparse(DifferentialAnalysis::DataType);
	void dense(DifferentialAnalysis::DataType);

	void scrna_mode();

	void scatac_mode();

	void scmultiome_mode();

	void chromvar_mode();

	void cicero_mode();

	void bulkrna_mode();


public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_differential_analysis_ready(DifferentialAnalysis, QString);
};
