#include "EnrichWorker.h"

#include "EnrichmentUtility.h"
#include "phyper.h"

EnrichWorker::EnrichWorker(
	const QString& database,
	const QString& enrich_type,
	const QStringList& gene_names,
	const QString& ontology,
	const soap::Species& species,
	const QString& p_adjust_method,
	double p_threshold
) :
	mode_(WorkMode::EnrichPathway),
	database_(database),
	enrich_type_(enrich_type),
	gene_names_(gene_names),
	ontology_(ontology),
	species_(species),
	p_adjust_method_(p_adjust_method),
	p_threshold_(p_threshold)
{}

EnrichWorker::EnrichWorker(
	const QStringList& peak_names,
	const MotifPosition* motif_position,
	const QString& enrich_type,
	const QString& p_adjust_method,
	double p_threshold
) :
	mode_(WorkMode::EnrichMotif),
	peak_names_(peak_names),
	motif_position_(motif_position),
	enrich_type_(enrich_type),
	p_adjust_method_(p_adjust_method),
	p_threshold_(p_threshold)
{};

bool EnrichWorker::enrich_pathway() {

	CustomMatrix enrichment = enrich(
		this->database_,
		this->gene_names_,
		this->ontology_,
		this->species_,
		this->p_adjust_method_,
		this->p_threshold_);

	if (enrichment.is_empty()) {
		G_TASK_WARN("Enrichment in " + this->database_ + " database failed.");
		return false;
	}
	else {

		this->res_ = enrichment;

		this->res_name_ = this->database_ + " " + this->enrich_type_ + " " + this->ontology_;

		return true;
	}
};

// modified from signac FindMotifs

bool EnrichWorker::enrich_motif() {

	auto index = custom::valid_index_of(this->peak_names_, this->motif_position_->peak_names_);

	if (index.isEmpty()) {
		G_TASK_WARN("No peaks detected in motif position");
		return false;
	}

	Eigen::MatrixXi frequency = this->motif_position_->get_motif_matrix().cwiseMin(1);

	Eigen::MatrixXi query = frequency(index, Eigen::all);

	Eigen::ArrayXi observed = query.array().colwise().sum();
	Eigen::ArrayXi background = frequency.array().colwise().sum();

	Eigen::ArrayXd percent_observed = observed.cast<double>() / query.rows();
	Eigen::ArrayXd percent_background = background.cast<double>() / frequency.rows();
	Eigen::ArrayXd fold_enrichment = custom::divide(percent_observed, percent_background);

	const int n_motif = frequency.cols(), n_peak = frequency.rows(), n_query = query.rows();
	Eigen::ArrayXd p_value(n_motif);
	for (int i = 0; i < n_motif; ++i) {
		p_value[i] = phyper(observed[i] - 1, background[i], n_peak - background[i], n_query);
	}

	Eigen::ArrayXd p_adjusted = custom::adjust_p_value(p_value, this->p_adjust_method_);

	CustomMatrix ret(this->motif_position_->motif_names_);

	ret.update(METADATA_ENRICHMENT_P_VALUE, custom::cast<QVector>(p_value));
	ret.update(METADATA_ENRICHMENT_ADJUSTED_P_VALUE, custom::cast<QVector>(p_adjusted));
	ret.update(METADATA_ENRICHMENT_OBSERVED_MOTIF_PERCENTAGE, custom::cast<QVector>(percent_observed));
	ret.update(METADATA_ENRICHMENT_BACKGROUND_MOTIF_PERCENTAGE, custom::cast<QVector>(percent_background));
	ret.update(METADATA_ENRICHMENT_OBSERVED_MOTIF_COUNT, custom::cast<QVector>(observed));
	ret.update(METADATA_ENRICHMENT_BACKGROUND_MOTIF_COUNT, custom::cast<QVector>(background));
	ret.update(METADATA_ENRICHMENT_FOLD_CHANGE, custom::cast<QVector>(fold_enrichment));
	ret.row_slice(custom::less_than(p_adjusted, this->p_threshold_));

	if (ret.is_empty()) {
		G_TASK_WARN("No Results meet requirements.");

		return false;
	}
	else {

		this->res_ = ret;

		this->res_name_ =  "Motif " + this->enrich_type_;

		return true;
	}
};

bool EnrichWorker::work() {

	if (this->mode_ == WorkMode::EnrichPathway) {
		return this->enrich_pathway();
	}
	else {
		return this->enrich_motif();
	}
};

void EnrichWorker::run() {

	G_TASK_LOG("Enrichment start...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_enrichment_ready(this->res_, this->res_name_);

	G_TASK_LOG("Enrichment finished.");
	G_TASK_END;
};
