#include "PandoWorker.h"

#include "ItemDatabase.h"
#include "Custom.h"
#include "GenomeUtility.h"
#include "MotifLocateWorker.h"

#include "glm.h"

bool PandoWorker::work() {

	if (this->single_cell_multiome_->species_ != soap::Species::Human) {
		G_TASK_WARN("Now only human is supported.");
		return false;
	}

	bool success = ItemDatabase::read_item(FILE_HUMAN_GENOME_GENOMIC_RANGE_SIF, this->annotation_);

	if (!success) {
		G_TASK_WARN("Loading failed.");
		return false;
	}

	if (!this->initiate_grn()) {
		return false;
	}

	if (!this->infer_grn()) {
		return false;
	}

	return true;
};

void PandoWorker::run() {

	G_TASK_LOG("Start pando...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_pando_ready(this->res_);

	G_TASK_LOG("Pando finished.");

	G_TASK_END;
};

bool PandoWorker::initiate_grn() {

	this->peaks_ = custom::stringlist_to_genomic_range(this->single_cell_multiome_->atac_counts()->rownames_);

	if (this->peaks_.size() == 0) {
		G_TASK_WARN("Invalid peak name.");
		return false;
	}

	// get exon
	auto exon_range = this->annotation_.row_sliced(
		custom::equal(this->annotation_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_TYPE), QString("exon")));

	if (exon_range.size() == 0) {
		G_TASK_WARN("Invalid annotation file.");
		return false;
	}

	exon_range.merge_nearby_range();

	this->peaks_.subtract(exon_range);

	if (this->peaks_.size() == 0) {
		G_TASK_WARN("No peak left after remove exon.");
		return false;
	}

	MotifLocateWorker worker(this->peaks_.get_range_names(), this->motif_database_, soap::Species::Human);

	if (!worker.work()) {
		G_TASK_WARN("Motif location Failed.");
		return false;
	}

	this->motif_position_ = worker.res_;

	this->peak_index_.set(this->peaks_);

	return true;
};

QStringList PandoWorker::find_variable_features(int n_feature) {

	Eigen::SparseMatrix<double> mat = custom::normalize(this->single_cell_multiome_->rna_counts()->mat_, 10000.0);

	const int ncol = mat.cols();
	const int nrow = mat.rows();
	double clipmax = std::sqrt(ncol);
	Eigen::ArrayXd row_mean = custom::row_mean(mat);
	Eigen::ArrayXd row_var = custom::row_var(mat);

	Eigen::ArrayX<bool> not_const = row_var > 0;
	const int not_const_row = not_const.count();
	if (not_const_row < n_feature) {
		return QStringList();
	}
	auto not_const_index = custom::which(not_const);

	Eigen::ArrayXd not_const_row_means = custom::sliced(row_mean, not_const);
	Eigen::ArrayXd not_const_row_var = custom::sliced(row_var, not_const);

	Eigen::ArrayXd fitted = custom::loess_mt(log10(not_const_row_var), log10(not_const_row_means), 2, 0.3, 50);

	Eigen::ArrayXd expected_var = Eigen::ArrayXd::Zero(nrow);
	expected_var(not_const_index) = pow(10, fitted);

	Eigen::ArrayXd expected_sd = sqrt(expected_var);
	Eigen::ArrayXd standardized_row_var = Eigen::ArrayXd::Zero(nrow);
	Eigen::ArrayXi n_zero = Eigen::ArrayXi::Constant(nrow, ncol);

	for (int i = 0; i < ncol; ++i) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
		{
			--n_zero[it.row()];
			standardized_row_var[it.row()] += std::pow(std::min(clipmax, std::abs((it.value() - row_mean[it.row()]) / expected_sd[it.row()])), 2);
		}
	}

	for (int i = 0; i < nrow; ++i) {
		if (expected_sd[i] == 0)continue;
		standardized_row_var[i] += std::pow(row_mean[i] / expected_sd[i], 2) * n_zero[i];
	}
	standardized_row_var /= (ncol - 1);
	Eigen::ArrayXi standardized_row_var_sorted_index = custom::order(standardized_row_var, true);
	Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(nrow, false);

	for (int i = 0; i < n_feature; ++i) {
		int index = standardized_row_var_sorted_index[i];
		filter[index] = true;
	}

	return custom::sliced(this->single_cell_multiome_->rna_counts()->rownames_, filter);
};

// bug 
bool PandoWorker::infer_grn() {

	auto gene_location = custom::get_hg38_gene_location();
	auto gene_names = this->find_variable_features(this->n_feature_);

	auto gene_filter = custom::in(gene_location.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME), gene_names);
	
	if (gene_filter.count() == 0) {
		G_TASK_WARN("Error in HG38 gene location File.");
		return false;
	}
	QStringList normalized_genes = this->single_cell_multiome_->rna_normalized()->rownames_;
	QStringList normalized_peaks = this->single_cell_multiome_->atac_normalized()->rownames_;
	gene_location.row_slice(gene_filter);
	gene_names = gene_location.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME);
	auto gene_index = custom::index_of(gene_names, normalized_genes);
	const int n_gene = gene_location.size();

	QVector<QVector<int>> gene_to_peak(n_gene);
#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {

		auto [seq_name, start, width, strand] = gene_location[i];

		int end{ 0 };
		if (this->use_tss_) {
			int tss{ 0 };
			if (strand == '-') {
				tss = start + width;
			}
			else {
				tss = start;
			}

			start = tss - this->upstream_;
			end = tss + this->downstream_;
		}
		else {
			end = start + width + this->downstream_;
			start -= this->upstream_;
		}

		if (start < 0) {
			start = 0;
		}

		gene_to_peak[i] = this->peak_index_.find_overlap_ranges(seq_name, start, end);
	}


	auto& peak_index = this->peaks_.metadata_.get_const_integer_reference(METADATA_GENOMIC_RANGE_PEAK_INDEX);

	auto& peak_data = this->single_cell_multiome_->atac_normalized()->mat_;
	auto& gene_data = this->single_cell_multiome_->rna_normalized()->mat_;

	auto&& peak_names = this->single_cell_multiome_->atac_normalized()->rownames_;

	QVector<QString> gene_name_res;
	QVector<QString> peak_name_res;
	QVector<QString> tf_name_res;
	QVector<double> rsq_res;
	QVector<double> estimate_res;
	QVector<double> std_error_res;
	QVector<double> t_val_res;
	QVector<double> p_val_res;
	QVector<int> n_variable_res;

#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {

		auto& peak_near = gene_to_peak[i];

		if (peak_near.isEmpty()) {
			continue;
		}

		int n_peak = peak_near.size();

		QVector<int> peak_reserve;
		QVector<double> peak_reserve_cor;

		Eigen::ArrayXd _gene = gene_data.row(gene_index[i]);

		if (this->filter_cell_) {
			_gene = custom::sliced(_gene, this->cell_filter_);
		}
		QString gene_name = normalized_genes[gene_index[i]];

		for (int j = 0; j < n_peak; ++j) {

			Eigen::ArrayXd _peak = peak_data.row(peak_index[peak_near[j]]);
			if (this->filter_cell_) {
				_peak = custom::sliced(_peak, this->cell_filter_);
			}

			double correlation = custom::correlation_pearson(_gene, _peak);

			if (std::abs(correlation) > this->peak_correlation_threshold_) {
				peak_reserve << peak_near[j];
				peak_reserve_cor << correlation;
			}
		}

		if (peak_reserve.isEmpty()) {
			continue;
		}

		auto peak_tfs = custom::sapply(peak_reserve,
			[this, &peak_index](auto ind) {return this->motif_position_.peak_to_tf(ind); }
		);

		if (custom::sum(custom::sapply(peak_tfs, [](auto&& motif) {return motif.size(); })) == 0) {
			continue;
		}

		auto tf_names = custom::unique(custom::unroll(peak_tfs));

		tf_names.removeOne(gene_names[i]);

		if (tf_names.isEmpty()) {
			continue;
		} 
		tf_names = custom::intersect(tf_names, normalized_genes);

		if (tf_names.isEmpty()) {
			continue;
		}
		auto tf_index = custom::index_of(tf_names, normalized_genes);
		int n_tf = tf_index.size();

		QStringList tf_use;
		QVector<double> tf_use_correlation;

		for (int j = 0; j < n_tf; ++j) {

			Eigen::ArrayXd _tf = gene_data.row(tf_index[j]);
			if (this->filter_cell_) {
				_tf = custom::sliced(_tf, this->cell_filter_);
			}

			double correlation = custom::correlation_pearson(_gene, _tf);

			if (std::abs(correlation) > this->motif_correlation_threshold_) {
				tf_use << tf_names[j];
				tf_use_correlation << correlation;
			}
		}

		if (tf_use.isEmpty()) {
			continue;
		}
		n_peak = peak_reserve.size();

		const int n_peak_tf_pair = custom::sum(custom::sapply(peak_tfs, [&tf_use](auto&& tfs)
		{return custom::intersect(tfs, tf_use).size(); }
		));

		if (n_peak_tf_pair == 0) {
			continue;
		}

		int peak_tf_pair_count{ 0 };
		const int n_cell = _gene.size();
		if (n_cell < n_peak_tf_pair) {
			continue;
		}

		Eigen::MatrixXd fit_x(n_cell, n_peak_tf_pair);

		QStringList gene_peak_names;
		QStringList gene_tf_names;
		QVector<double> gene_peak_cor;
		QVector<double> gene_tf_cor;
		for (int j = 0; j < n_peak; ++j) {


			auto peak_tf_use = custom::intersect(peak_tfs[j], tf_use);

			if (peak_tf_use.isEmpty()) {
				continue;
			}

			int peak_ind = peak_index[peak_reserve[j]];

			QString peak_name = normalized_peaks[peak_ind];

			Eigen::ArrayXd _peak = peak_data.row(peak_ind);
			if (this->filter_cell_) {
				_peak = custom::sliced(_peak, this->cell_filter_);
			}

			int n_peak_tf_use = peak_tf_use.size();

			gene_peak_names << QStringList(n_peak_tf_use, peak_name);
			gene_peak_cor << QVector<double>(n_peak_tf_use, peak_reserve_cor[j]);

			for (int k = 0; k < n_peak_tf_use; ++k) {

				int tf_ind = tf_use.indexOf(peak_tf_use[k]);

				gene_tf_cor << tf_use_correlation[tf_ind];

				int ind = normalized_genes.indexOf(peak_tf_use[k]);

				gene_tf_names << peak_tf_use[k];

				Eigen::ArrayXd _tf = gene_data.row(ind);
				if (this->filter_cell_) {
					_tf = custom::sliced(_tf, this->cell_filter_);
				}

				fit_x.col(peak_tf_pair_count++) = _tf * _peak;
			}
		}

		auto res = glm_gaussian(_gene, fit_x);

		double rsq = 1.0 - res.deviance / res.null_deviance;

		bool hasnan = res.coefficients.hasNaN();
		if (hasnan) {
			continue;
		}

		auto summary = glm_summary(res);
		if (summary.hasnan) {
			continue;
		}

	#pragma omp critical
		{
			gene_name_res << QStringList(n_peak_tf_pair, gene_name);
			peak_name_res << gene_peak_names;
			tf_name_res << gene_tf_names;
			rsq_res << QVector<double>(n_peak_tf_pair, rsq);
			estimate_res << custom::cast<QVector>(summary.estimate.segment(1, n_peak_tf_pair));
			std_error_res << custom::cast<QVector>(summary.std_err.segment(1, n_peak_tf_pair));
			t_val_res << custom::cast<QVector>(summary.t_value.segment(1, n_peak_tf_pair));
			p_val_res << custom::cast<QVector>(summary.p.segment(1, n_peak_tf_pair));
			n_variable_res << QVector<int>(n_peak_tf_pair, fit_x.cols());
		}
	}

	if (gene_name_res.isEmpty()) {
		G_TASK_WARN("No Results.");
		return false;
	}

	this->res_.mat_.update(METADATA_PANDO_GENE_NAME, gene_name_res);
	this->res_.mat_.update(METADATA_PANDO_PEAK_NAME, peak_name_res);
	this->res_.mat_.update(METADATA_PANDO_TF_NAME, tf_name_res);
	this->res_.mat_.update(METADATA_PANDO_R_SQUARED_NAME, rsq_res);
	this->res_.mat_.update(METADATA_PANDO_ESTIMATE_NAME, estimate_res);
	this->res_.mat_.update(METADATA_PANDO_STD_ERR_NAME, std_error_res);
	this->res_.mat_.update(METADATA_PANDO_T_VAL_NAME, t_val_res);
	this->res_.mat_.update(METADATA_PANDO_P_VAL_NAME, p_val_res);
	this->res_.mat_.update(METADATA_PANDO_N_VARIABLE_NAME, n_variable_res);

	return true;	
};