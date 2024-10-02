#include "GseaWorker.h"

#include "Custom.h"
#include <QtConcurrent>

GseaWorker::GseaWorker(
	const SparseDouble& sd,
	const QStringList& metadata,
	const QStringList& comparison,
	soap::Species species,
	const QString& database,
	int minimum_overlap,
	int minimum_size,
	int maximum_size,
	int random_state,
	const QString& permutation_type,
	int n_permutation
) :
	sd_(sd),
	metadata_(metadata),
	comparison_(comparison),
	species_(species),
	database_(database),
	minimum_overlap_(minimum_overlap),
	minimum_size_(minimum_size),
	maximum_size_(maximum_size),
	random_state_(random_state),
	permutation_type_(permutation_type),
	n_permutation_(n_permutation)
{}

bool GseaWorker::load_database() {

	QString database_file;

	if (this->species_ == soap::Species::Human) {
		if (this->database_ == "Curated") {
			database_file = FILE_HUMAN_GSEA_CURATED;
		}
		else {
			database_file = FILE_HUMAN_GSEA_ONTOLOGY;
		}
	}
	else {
		if (this->database_ == "Curated") {
			database_file = FILE_MOUSE_GSEA_CURATED;
		}
		else {
			database_file = FILE_MOUSE_GSEA_ONTOLOGY;
		}
	}

	QStringList database_symbols;

	QFile file(database_file);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QString line = in.readLine();
	while (!line.isNull()) {
		QStringList tmp = line.split('\t');
		QString path = tmp[0];
		QStringList gene_names = tmp.sliced(2);

		int size = gene_names.size();
		if (size >= this->minimum_size_ && size <= this->maximum_size_) {

			int overlap = _Cs intersect_length(gene_names, this->sd_.rownames_);

			if (overlap >= this->minimum_overlap_) {
				this->pathway_to_symbol_[path] = gene_names;
				this->pathway_to_size_[path] = size;
				this->pathway_list_ << path;

				database_symbols << gene_names;
			}
		}
		line = in.readLine();
	}

	if (this->pathway_to_symbol_.isEmpty()) {
		G_TASK_NOTICE("No Pathway is valid.");
		return false;
	}

	database_symbols = _Cs unique(database_symbols);
	auto filter = _Cs in(this->sd_.rownames_, database_symbols);

	if (filter.count() < 500) {
		return false;
	}

	this->sd_.row_slice(filter);

	auto fun2 = [this](const QString& path) {
		QVector<int> loc = _Cs valid_index_of(this->pathway_to_symbol_[path], this->sd_.rownames_);

		this->mutex_.lock();
		this->pathway_to_location_[path] = loc;
		this->mutex_.unlock();
	};

	QStringList params = this->pathway_to_symbol_.keys();
	auto f = QtConcurrent::map(params, fun2);
	f.waitForFinished();

	this->gene_locations_ = _Cs seq_n(0, this->sd_.rownames_.size());

	return true;
};

void GseaWorker::preprocess() {

	Eigen::ArrayX<bool> group1 = _Cs equal(this->metadata_, this->comparison_[1]), group2;

	if (this->comparison_[2] != "REST") {
		group2 = _Cs equal(this->metadata_, this->comparison_[2]);
	}
	else {
		group2 = !group1;
	}

	Eigen::ArrayX<bool> group12 = group1 + group2;
	Eigen::ArrayX<bool> filter = (_Cs col_sliced_row_sum(this->sd_.mat_, group12) > 0);

	this->sd_.row_slice(filter);

	auto order = _Cs which(group1) << _Cs which(group2);

	this->sd_.col_reorder(order);
	this->metadata_ = _Cs reordered(this->metadata_, order);
};


void GseaWorker::get_order() {
	Eigen::ArrayX<bool> group1 = _Cs equal(this->metadata_, this->comparison_[1]), group2;

	if (this->comparison_[2] != "REST") {
		group2 = _Cs equal(this->metadata_, this->comparison_[2]);
	}
	else {
		group2 = !group1;
	}

	Eigen::ArrayXXd exp1 = _Cs col_sliced(this->sd_.mat_, group1).toDense().array();
	Eigen::ArrayXXd exp2 = _Cs col_sliced(this->sd_.mat_, group2).toDense().array();

	const int gene_number = this->sd_.mat_.rows();

	Eigen::ArrayXd signal_to_noise(gene_number), mean1 = exp1.rowwise().mean(), mean2 = exp2.rowwise().mean();
	Eigen::ArrayXd std1 = ((exp1.colwise() -= mean1).rowwise().squaredNorm().array() / exp1.cols()).sqrt(),
		std2 = ((exp2.colwise() -= mean2).rowwise().squaredNorm().array() / exp2.cols()).sqrt();
	Eigen::ArrayXd mean12 = mean1 - mean2, std12 = std1 + std2;

	for (int i = 0; i < gene_number; ++i) {
		if (std12[i] == 0) {
			if (mean1[i] > mean2[i]) {
				signal_to_noise[i] = 1;
			}
			else if (mean1[i] < mean2[i]) {
				signal_to_noise[i] = -1;
			}
			else {
				signal_to_noise[i] = 0;
			}
		}
		else {
			signal_to_noise[i] = mean12[i] / std12[i];
		}
	}

	this->sorted_gene_list_ = _Cs reordered(this->gene_locations_, _Cs order(signal_to_noise, true));
	this->correlations_ = _Cs sorted(signal_to_noise, true);
	this->weights_ = this->correlations_.abs();
};

void GseaWorker::calculate_enrichment_score() {

	for (auto&& path : this->pathway_list_) {

		auto gene_location = _Cs which(_Cs in(this->sorted_gene_list_, this->pathway_to_location_[path]));
		int n_gene = this->sorted_gene_list_.size(), n_hit = gene_location.size(), n_loss = n_gene - n_hit;

		Eigen::ArrayXd weights = this->weights_(gene_location);
		double weight_sum = weights.sum();

		Eigen::ArrayXd up = weights / weight_sum;
		Eigen::ArrayXd gaps(n_hit);

		gaps[0] = gene_location[0];

		for (int i = 1; i < n_hit; ++i) {
			gaps[i] = gene_location[i] - gene_location[i - 1] - 1;
		}

		Eigen::ArrayXd down = gaps / n_loss;
		Eigen::ArrayXd ridge_enrichment_score = _Cs cumsum(up - down);
		Eigen::ArrayXd valleys = ridge_enrichment_score - up;

		double maximum_enrichment_score = ridge_enrichment_score.maxCoeff(), minimum_enrichment_score = valleys.minCoeff();

		double enrichment_score = maximum_enrichment_score > -minimum_enrichment_score ? maximum_enrichment_score : minimum_enrichment_score;

		this->enrichment_score_[path] = enrichment_score;
		this->gene_location_[path] = gene_location;

		QVector<double> plot_point_x, plot_point_y;
		plot_point_x << -1;
		plot_point_y << 0;
		for (int i = 0; i < n_hit; ++i) {
			plot_point_x << gene_location[i] - 1 << gene_location[i];
			plot_point_y << valleys[i] << ridge_enrichment_score[i];
		}
		plot_point_x << n_gene - 1;
		plot_point_y << 0;
		this->point_x_[path] = plot_point_x;
		this->point_y_[path] = plot_point_y;
	}
};

void GseaWorker::permutation_score_phenotype(const QString& gene_set_name) {

	auto& gene_location = this->pathway_to_location_[gene_set_name];

#pragma omp parallel for
	for (int i = 0; i < this->n_permutation_; ++i) {

		auto& gene_list = this->permuted_gene_location_[i];
		Eigen::ArrayXd& weights = this->weights_list_[i];

		auto new_location = _Cs which(_Cs in(gene_list, gene_location));

		int n_gene = gene_list.size(), n_hit = new_location.size(), n_loss = n_gene - n_hit;
		double down_unit = 1 / (double)n_loss;

		Eigen::ArrayXd score = Eigen::ArrayXd::Constant(n_gene, -down_unit);
		Eigen::ArrayXd sub_weights = weights(new_location);

		sub_weights /= sub_weights.sum();
		score(new_location) = sub_weights;

		double maximum_enrichment_score = 0, minimum_enrichment_score = 0, current_enrichment_score = 0;

		for (int i = 0; i < n_gene; ++i) {
			current_enrichment_score += score[i];
			maximum_enrichment_score = maximum_enrichment_score > current_enrichment_score ?
				maximum_enrichment_score : current_enrichment_score;
			minimum_enrichment_score = minimum_enrichment_score < current_enrichment_score ?
				minimum_enrichment_score : current_enrichment_score;
		}

		double enrichment_score = maximum_enrichment_score > -minimum_enrichment_score ?
			maximum_enrichment_score : minimum_enrichment_score;

	#pragma omp critical
		{
			this->enrichment_scores_[gene_set_name] << enrichment_score;
		}
	}
};

void GseaWorker::gsea_gene_set() {

	int count{ 0 };

	for (const auto& path : this->pathway_list_) {

		auto gene_location = _Cs which(_Cs in(this->sorted_gene_list_, this->pathway_to_location_[path]));

		int n_gene = this->sorted_gene_list_.size(), n_hit = gene_location.size(), n_loss = n_gene - n_hit;

		Eigen::ArrayXd weights = this->weights_(gene_location);

		double weight_sum = weights.sum();

		Eigen::ArrayXd up = weights / weight_sum;
		Eigen::ArrayXd gaps(n_hit);

		gaps[0] = gene_location[0];

		for (int i = 1; i < n_hit; ++i) {
			gaps[i] = gene_location[i] - gene_location[i - 1] - 1;
		}
		Eigen::ArrayXd down = gaps / n_loss;
		Eigen::ArrayXd ridge_enrichment_score = _Cs cumsum(up - down);
		Eigen::ArrayXd valleys = ridge_enrichment_score - up;

		double maximum_enrichment_score = ridge_enrichment_score.maxCoeff(), minimum_enrichment_score = valleys.minCoeff();

		double enrichment_score = maximum_enrichment_score > -minimum_enrichment_score ?
			maximum_enrichment_score : minimum_enrichment_score;

		this->enrichment_score_[path] = enrichment_score;
		this->gene_location_[path] = gene_location;

		QVector<double> plot_point_x, plot_point_y;
		plot_point_x << -1;
		plot_point_y << 0;
		for (int i = 0; i < n_hit; ++i) {
			plot_point_x << gene_location[i] - 1 << gene_location[i];
			plot_point_y << valleys[i] << ridge_enrichment_score[i];
		}
		plot_point_x << n_gene - 1;
		plot_point_y << 0;

		this->point_x_[path] = plot_point_x;
		this->point_y_[path] = plot_point_y;

		auto fun = [this, &path](int) {

			int size = this->sorted_gene_list_.size();
			auto& gene_location = this->pathway_to_location_[path];

			QVector<int> shuffled_location = this->gene_locations_;

			this->mutex_.lock();
			std::default_random_engine local_random_engine(this->unsigned_distribution_(this->random_engine_));
			this->mutex_.unlock();

			std::shuffle(shuffled_location.begin(), shuffled_location.end(), local_random_engine);

			auto new_location = _Cs which(_Cs in(shuffled_location, gene_location));

			int n_gene = size, n_hit = new_location.size(), n_loss = n_gene - n_hit;

			double down_unit = 1 / (double)n_loss;

			Eigen::ArrayXd score = Eigen::ArrayXd::Constant(n_gene, -down_unit);
			Eigen::ArrayXd sub_weights = this->weights_(new_location);

			sub_weights /= sub_weights.sum();
			score(new_location) = sub_weights;

			double maximum_enrichment_score = 0, minimum_enrichment_score = 0, current_enrichment_score = 0;

			for (int i = 0; i < n_gene; ++i) {
				current_enrichment_score += score[i];
				maximum_enrichment_score = maximum_enrichment_score > current_enrichment_score ?
					maximum_enrichment_score : current_enrichment_score;
				minimum_enrichment_score = minimum_enrichment_score < current_enrichment_score ?
					minimum_enrichment_score : current_enrichment_score;
			}
			double enrichment_score = maximum_enrichment_score > -minimum_enrichment_score ?
				maximum_enrichment_score : minimum_enrichment_score;

			this->mutex2_.lock();
			this->enrichment_scores_[path] << enrichment_score;
			this->mutex2_.unlock();
		};
		auto params = _Cs seq_n(0, this->n_permutation_);
		QFuture<void> f = QtConcurrent::map(params, fun);
		f.waitForFinished();
	}
};

void GseaWorker::permutate_phenotype() {

	std::shuffle(this->metadata_.begin(), this->metadata_.end(), this->random_engine_);

	Eigen::ArrayX<bool> group1 = _Cs equal(this->metadata_, this->comparison_[1]), group2;
	if (this->comparison_[2] != "REST") {
		group2 = _Cs equal(this->metadata_, this->comparison_[2]);
	}
	else {
		group2 = !group1;
	}
	Eigen::ArrayXXd exp1, exp2;
	Eigen::ArrayXd mean1, mean2, std1, std2;
	const int gene_number = this->sd_.mat_.rows();
#pragma omp sections
	{
	#pragma omp section
		{
			exp1 = _Cs col_sliced(this->sd_.mat_, group1).toDense().array();
			mean1 = exp1.rowwise().mean();
			std1 = ((exp1.colwise() -= mean1).rowwise().squaredNorm().array() / exp1.cols()).sqrt();
		}
	#pragma omp section
		{
			exp2 = _Cs col_sliced(this->sd_.mat_, group2).toDense().array();
			mean2 = exp2.rowwise().mean();
			std2 = ((exp2.colwise() -= mean2).rowwise().squaredNorm().array() / exp2.cols()).sqrt();
		}
	}

	Eigen::ArrayXd signal_to_noise(gene_number);
	Eigen::ArrayXd mean12 = mean1 - mean2, std12 = std1 + std2;
	for (int i = 0; i < gene_number; ++i) {
		if (std12[i] == 0) {
			if (mean1[i] > mean2[i]) {
				signal_to_noise[i] = 1;
			}
			else if (mean1[i] < mean2[i]) {
				signal_to_noise[i] = -1;
			}
			else {
				signal_to_noise[i] = 0;
			}
		}
		else {
			signal_to_noise[i] = mean12[i] / std12[i];
		}
	}

	this->permuted_gene_location_ << _Cs reordered(this->gene_locations_, _Cs order(signal_to_noise, true));
	Eigen::ArrayXd correlations = _Cs sorted(signal_to_noise, true);
	this->weights_list_ << correlations.abs();
};


void GseaWorker::calculate_p_value() {
	for (auto path : this->pathway_to_symbol_.keys()) {

		double pathway_enrichment_score = this->enrichment_score_[path];
		auto& enrichment_scores = this->enrichment_scores_[path];
		int size = enrichment_scores.size();

		if (pathway_enrichment_score >= 0) {
			QVector<double> pos;
			for (int i = 0; i < size; ++i) {
				if (enrichment_scores[i] >= 0) {
					pos << enrichment_scores[i];
				}
			}
			if (pos.isEmpty()) {
				this->normalized_enrichment_score_[path] = 2 * pathway_enrichment_score;
				enrichment_scores = _Cs multiply(enrichment_scores, 2.0);

				if (pathway_enrichment_score > 0.5) {
					this->p_values_[path] = 1 / (double)this->n_permutation_;
				}
				else {
					this->p_values_[path] = 1;
				}
			}
			else {
				double mean_enrichment_score = _Cs mean(pos);

				this->normalized_enrichment_score_[path] = pathway_enrichment_score / mean_enrichment_score;
				enrichment_scores = _Cs divide(enrichment_scores, mean_enrichment_score);

				this->p_values_[path] = _Cs greater_equal(pos, pathway_enrichment_score).count() / (double)pos.size();
			}
		}
		else {
			QVector<double> neg;
			for (int i = 0; i < size; ++i) {
				if (enrichment_scores[i] < 0) {
					neg << enrichment_scores[i];
				}
			}

			if (neg.isEmpty()) {
				this->normalized_enrichment_score_[path] = 2 * pathway_enrichment_score;
				enrichment_scores = _Cs multiply(enrichment_scores, 2.0);
				if (pathway_enrichment_score < -0.5) {
					this->p_values_[path] = 1 / (double)this->n_permutation_;
				}
				else {
					this->p_values_[path] = 1;
				}
			}
			else {
				double mean_enrichment_score = _Cs mean(neg);

				this->normalized_enrichment_score_[path] = pathway_enrichment_score / (-mean_enrichment_score);
				enrichment_scores = _Cs divide(enrichment_scores, -mean_enrichment_score);

				this->p_values_[path] = _Cs less_equal(neg, pathway_enrichment_score).count() / (double)neg.size();
			}
		}
	}
};

void GseaWorker::calculate_fwer() {

	QVector<double> maximum_enrichment_score_positive(this->n_permutation_, 0), 
		minimum_enrichment_score_negative(this->n_permutation_, 0);

	for (const auto& path : this->pathway_list_) {
		auto& enrichment_scores = this->enrichment_scores_[path];
		int size = enrichment_scores.size();
		for (int i = 0; i < size; ++i) {
			if (maximum_enrichment_score_positive[i] < enrichment_scores[i]) {
				maximum_enrichment_score_positive[i] = enrichment_scores[i];
			}
			if (minimum_enrichment_score_negative[i] > enrichment_scores[i]) {
				minimum_enrichment_score_negative[i] = enrichment_scores[i];
			}
		}
	}

	for (const auto& path : this->pathway_list_) {
		double pathway_normalized_enrichment_score = this->normalized_enrichment_score_[path], fwer;
		if (pathway_normalized_enrichment_score >= 0) {
			fwer = _Cs greater_equal(maximum_enrichment_score_positive, pathway_normalized_enrichment_score).count() /
				(double)this->n_permutation_;
		}
		else {
			fwer = _Cs less_equal(minimum_enrichment_score_negative, pathway_normalized_enrichment_score).count() / 
				(double)this->n_permutation_;
		}
		this->fwer_[path] = fwer;
	}
}

void GseaWorker::calculate_fdr() {

	QVector<double> all_scores;
	for (auto&& s : this->enrichment_scores_) {
		all_scores << s;
	}

	QVector<double> all_nes = this->normalized_enrichment_score_.values();

	_Cs sort(all_scores);
	_Cs sort(all_nes);

	int n_score = all_scores.size();
	int n_nes = all_nes.size();
	int all_pos = std::ranges::count(all_scores, true, [](double score) {return score >= 0.0; });
	int all_neg = n_score - all_pos;

	int nes_pos = std::ranges::count(all_nes, true, [](double score) {return score >= 0.0; });
	int nes_neg = n_nes - nes_pos;

	for (const auto &path : this->pathway_list_) {
		double pathway_normalized_enrichment_score = this->normalized_enrichment_score_[path];

		if (pathway_normalized_enrichment_score >= 0.0) {

			if (all_pos == 0) {
				this->fdr_[path] = 1.0;
				continue;
			}

			double g_score = std::ranges::count(all_scores, true, 
				[pathway_normalized_enrichment_score](double score) {return score >= pathway_normalized_enrichment_score; });
			double g_nes = std::ranges::count(all_nes, true,
				[pathway_normalized_enrichment_score](double score) {return score >= pathway_normalized_enrichment_score; });

			double pi_norm = g_score / all_pos;
			double pi_obs = g_nes / nes_pos;

			if (pi_obs == 0.0) {
				this->fdr_[path] = 1.0;
			}
			else {
				double fwer = pi_norm / pi_obs;
				if (fwer > 1.0) {
					fwer = 1.0;
				}

				this->fdr_[path] = fwer;
			}
		}
		else {
			if (all_neg == 0) {
				this->fdr_[path] = 1.0;
				continue;
			}

			double l_score = std::ranges::count(all_scores, true,
				[pathway_normalized_enrichment_score](double score) {return score < pathway_normalized_enrichment_score; });
			double l_nes = std::ranges::count(all_nes, true,
				[pathway_normalized_enrichment_score](double score) {return score < pathway_normalized_enrichment_score; });

			double pi_norm = l_score / all_neg;
			double pi_obs = l_nes / nes_neg;

			if (pi_obs == 0.0) {
				this->fdr_[path] = 1.0;
			}
			else {
				double fwer = pi_norm / pi_obs;
				if (fwer > 1.0) {
					fwer = 1.0;
				}

				this->fdr_[path] = fwer;
			}
		}
	}

}

void GseaWorker::gsea_phenotype() {
	this->calculate_enrichment_score();

	for (int i = 0; i < this->n_permutation_; ++i) {
		this->permutate_phenotype();
	}

	for (const auto& path : this->pathway_list_) {
		this->permutation_score_phenotype(path);
	}
};

void GseaWorker::run() {

	this->random_engine_.seed(this->random_state_);

	this->preprocess();

	if (!this->load_database()) {
		G_TASK_END;
	}

	this->get_order();

	this->original_correlations_ = this->correlations_;

	if (this->permutation_type_ == "Phenotype") {
		this->gsea_phenotype();
	}
	else {
		this->gsea_gene_set();
	}

	this->calculate_p_value();

	this->calculate_fdr();

	this->calculate_fwer();

	emit x_gsea_ready(GSEA(
		this->database_,
		this->comparison_,
		this->enrichment_score_,
		this->normalized_enrichment_score_,
		this->p_values_,
		this->fdr_,
		this->fwer_,
		this->pathway_to_size_,
		_Cs cast<QVector>(this->original_correlations_),
		this->gene_location_,
		this->point_x_,
		this->point_y_
	));

	G_TASK_END;
}
