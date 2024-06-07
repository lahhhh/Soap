#include "InferCnvWorker.h"

#include "FileIO.h"

#include "Custom.h"
#include "nls.h"

#include "PcaWorker.h"
#include "LeidenPartitionWorker.h"

#include "lm.h"
#include "pval.h"

void InferCnvWorker::finalize() {

	QVector<int> cell_ind;
	std::vector<std::tuple<QString, int, int> > cell_location;

	int cell_count{ 0 };

	for (auto&& [type, ind] : this->ref_grouped_cell_indices_) {

		auto subcluster = _Cs reordered(this->subclusters_, ind);

		auto clusters = _Cs unique(subcluster);

		for (auto&& cluster : clusters) {

			auto ind2 = _Cs match(subcluster, cluster);

			cell_ind << _Cs reordered(ind, ind2);
		}

		int n_type_cell = ind.size();

		cell_location.emplace_back(std::tuple{ type, cell_count, n_type_cell });

		cell_count += n_type_cell;
	}

	for (auto&& [type, ind] : this->obs_grouped_cell_indices_) {

		auto subcluster = _Cs reordered(this->subclusters_, ind);

		auto clusters = _Cs unique(subcluster);

		for (auto&& cluster : clusters) {

			auto ind2 = _Cs match(subcluster, cluster);

			cell_ind << _Cs reordered(ind, ind2);
		}

		int n_type_cell = ind.size();

		cell_location.emplace_back(std::tuple{ type, cell_count, n_type_cell });

		cell_count += n_type_cell;
	}

	QString chr = this->chrs_[0];
	int loc{ 0 };
	int length{ 1 };

	std::vector<std::tuple<QString, int, int> > chromosome_location;

	int n_gene = this->chrs_.size();

	QString now = chr;

	for (int i = 1; i < n_gene; ++i) {
		if (this->chrs_[i] == now) {
			++length;
		}
		else {
			chromosome_location.emplace_back(std::tuple{ now, loc, length });
			now = this->chrs_[i];
			loc += length;
			length = 1;
		}
	}

	this->expr_ = this->expr_(Eigen::all, cell_ind).eval();

	chromosome_location.emplace_back(std::tuple{ now, loc, length });

	CNV* cnv = new CNV(this->expr_, cell_location, chromosome_location, CNV::DataType::InferCnv);

	emit x_cnv_ready(cnv);

};

void InferCnvWorker::run() {

	if (this->species_ != soap::Species::Human && this->species_ != soap::Species::Mouse) {
		G_TASK_WARN("Unsupport Species.");
		G_TASK_END;
	}

	if (!this->prepare_0()) {
		G_TASK_END;
	}
	if (!this->remove_insufficiently_expressed_genes_1()) {
		G_TASK_END;
	}
	if (!this->normalization_by_sequencing_depth_2()) {
		G_TASK_END;
	}
	if (!this->log_transformation_3()) {
		G_TASK_END;
	}
	if (!this->subtract_average_reference_5_8()) {
		G_TASK_END;
	}
	if (!this->smooth_6()) {
		G_TASK_END;
	}
	if (!this->center_7()) {
		G_TASK_END;
	}
	if (!this->subtract_average_reference_5_8()) {
		G_TASK_END;
	}
	if (!this->invert_log_transform_10()) {
		G_TASK_END;
	}
	if (!this->compute_tumor_subcluster_11()) {
		G_TASK_END;
	}
	this->finalize();
	G_TASK_END;
};

bool InferCnvWorker::prepare_0() {

	std::unique_ptr<CustomMatrix> gene_location;

	if (this->species_ == soap::Species::Human) {
		gene_location.reset(read_sv(FILE_CNV_GENE_LOCATION_HUMAN, '\t'));
	}
	else if (this->species_ == soap::Species::Mouse) {
		gene_location.reset(read_sv(FILE_CNV_GENE_LOCATION_MOUSE, '\t'));
	}
	else {
		G_TASK_WARN("InferCNV now only support human and mouse data.");
		return false;
	}

	if (gene_location == nullptr) {
		G_TASK_WARN("Resource File Loading Failed.");
		return false;
	}

	QStringList gene_list = gene_location->get_qstring("Column 1");
	QStringList chromosome_list = gene_location->get_qstring("Column 2");
	QVector<int> starts = _Cs cast<int>(gene_location->get_qstring("Column 3"));
	QVector<int> stops = _Cs cast<int>(gene_location->get_qstring("Column 4"));
	std::ranges::for_each(chromosome_list, [](auto&& t) {t = _Cs standardize_chromosome_name(t); });

	QStringList unique_chromosome = _Cs unique(chromosome_list);

	if (this->species_ == soap::Species::Human) {
		unique_chromosome = _Cs sliced(soap::HumanChromosomeOrder, _Cs in(soap::HumanChromosomeOrder, unique_chromosome));
	}
	else if (this->species_ == soap::Species::Mouse) {
		unique_chromosome = _Cs sliced(soap::MouseChromosomeOrder, _Cs in(soap::MouseChromosomeOrder, unique_chromosome));
	}

	QVector<int> gene_index;

	for (const auto& chr : unique_chromosome) {
		auto chr_filter = _Cs equal(chromosome_list, chr);
		QStringList chr_gene_list = _Cs sliced(gene_list, chr_filter);
		auto gene_filter = _Cs in(chr_gene_list, this->counts_.rownames_);
		int n_gene = gene_filter.count();
		if (n_gene > 20) {

			this->chrs_ << QStringList(n_gene, chr);

			gene_index << _Cs valid_index_of(chr_gene_list, this->counts_.rownames_);

			this->starts_ << _Cs sliced(_Cs sliced(starts, chr_filter), gene_filter);

			this->stops_ << _Cs sliced(_Cs sliced(stops, chr_filter), gene_filter);
		}
	}

	if (gene_index.isEmpty()) {
		G_TASK_WARN("No Enough Genes Found.");
		return false;
	}

	this->counts_.row_reorder(gene_index);

	Eigen::ArrayX<bool> filter = _Cs col_sum(this->counts_.mat_) > this->min_counts_per_cell_;

	if (filter.count() < 100) {
		G_TASK_WARN("Two few cells meeting requirements.");
		return false;
	}

	this->counts_.col_slice(filter);

	this->metadata_ = _Cs sliced(this->metadata_, filter);
	this->ref_group_names_ = _Cs intersect(this->ref_group_names_, _Cs unique(this->metadata_));

	if (this->ref_group_names_.isEmpty()) {
		G_TASK_WARN("No Qualified Reference Cells.");
		return false;
	}

	this->obs_group_names_ = _Cs set_difference(_Cs unique(this->metadata_), this->ref_group_names_);
	
	for (auto&& level : this->ref_group_names_) {
		this->ref_grouped_cell_indices_[level] = _Cs match(this->metadata_, level);
	}

	for (auto&& level : this->obs_group_names_) {
		this->obs_grouped_cell_indices_[level] = _Cs match(this->metadata_, level);
	}

	for (auto&& [type, ind] : this->ref_grouped_cell_indices_) {
		if (ind.size() < 30) {
			G_TASK_WARN("Cluster " + type + " has too few cells. Consider remove it from query data.");
			return false;
		}
	}

	for (auto&& [type, ind] : this->obs_grouped_cell_indices_) {
		if (ind.size() < 30) {
			G_TASK_WARN("Cluster " + type + " has too few cells. Consider remove it from query data.");
			return false;
		}
	}

	return true;
};

bool InferCnvWorker::remove_insufficiently_expressed_genes_1() {

	Eigen::ArrayX<bool> filter1 = _Cs row_mean(this->counts_.mat_) >= this->cut_off_;

	Eigen::ArrayX<bool> filter2 = _Cs row_count(this->counts_.mat_, 1) > this->min_cells_per_gene_;

	Eigen::ArrayX<bool> filter = filter1 * filter2;

	if (filter.count() < 100) {
		G_TASK_WARN("Less Than 100 genes passed filter. Task Terminated.");
		return false;
	}

	this->counts_.row_slice(filter);
	this->chrs_ = _Cs sliced(this->chrs_, filter);
	this->starts_ = _Cs sliced(this->starts_, filter);
	this->stops_ = _Cs sliced(this->stops_, filter);
	return true;
};

bool InferCnvWorker::normalization_by_sequencing_depth_2() {

	this->expr_ = this->counts_.mat_.toDense().cast<double>();

	this->counts_.clear();

	_Cs normalize_in_place(this->expr_);

	if (!this->hmm_) {
		return true;
	}

	// add hspike

	//int n_cells{ 100 };
	//int n_genes_per_chr{ 400 };
	//int n_total_gene = this->expr_.rows();

	//int n_remaining = n_total_gene - 10 * n_genes_per_chr;
	//if (n_remaining < n_genes_per_chr) {
	//	n_remaining = n_genes_per_chr;
	//}

	//this->fake_chr_info_ << std::make_tuple(QString{ "chrA" }, 1.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chrA");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chr_0" }, 0.01, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chr_0");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chrB" }, 1.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chrB");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chr_0pt5" }, 0.5, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chr_0pt5");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chrC" }, 1.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chrC");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chr_1pt5" }, 1.5, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chr_1pt5");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chrD" }, 1.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chrD");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chr_2pt0" }, 2.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chr_2pt0");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chrE" }, 1.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chrE");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chr_3pt0" }, 3.0, n_genes_per_chr);
	//this->fake_chrs_ << QStringList(n_genes_per_chr, "chr_3pt0");
	//this->fake_chr_info_ << std::make_tuple(QString{ "chrF" }, 1.0, n_remaining);
	//this->fake_chrs_ << QStringList(n_remaining, "chrF");

	//int n_genes = _Cs sum(_Cs sapply(this->fake_chr_info_, [](auto&& t) {return std::get<2>(t); }));

	//auto genes_means_use_idx = _Cs sample_integer(0, this->expr_.rows(), n_genes);

	//Eigen::MatrixXd sim_counts_matrix;

	//int cell_counter{ 0 };
	//for (auto&& [type, index] : this->ref_grouped_cell_indices_) {
	//	Eigen::MatrixXd normal_cells_expr = this->expr_(Eigen::all, index);

	//	Eigen::ArrayXd gene_means_orig = normal_cells_expr.rowwise().mean();
	//	Eigen::ArrayXd gene_means = gene_means_orig(genes_means_use_idx);

	//	gene_means = (gene_means == 0.0).select(0.001, gene_means);

	//	Eigen::MatrixXd sim_normal_matrix; // get_simulated_cell_matrix_using_meanvar_trend

	//	Eigen::ArrayXd hspike_gene_means = gene_means;

	//	for (auto&& [name, cnv, _] : this->fake_chr_info_) {
	//		if (cnv != 1.0) {
	//			hspike_gene_means(_Cs match(this->fake_chrs_, name)) *= cnv;
	//		}
	//	}

	//	Eigen::MatrixXd sim_spiked_cnv_matrix;

	//	if (sim_counts_matrix.size() == 0) {
	//		sim_counts_matrix = _Cs cbind(sim_normal_matrix, sim_spiked_cnv_matrix);
	//	}
	//	else {
	//		sim_counts_matrix = _Cs cbind(sim_counts_matrix, sim_normal_matrix, sim_spiked_cnv_matrix);
	//	}

	//	this->fake_ref_grouped_cell_indices_[type] = _Cs seq_n(cell_counter, n_cells);

	//	cell_counter += n_cells;

	//	this->fake_obs_grouped_cell_indices_[type] = _Cs seq_n(cell_counter, n_cells);

	//	cell_counter += n_cells;

	//}

	//this->fake_expr_ = sim_counts_matrix;
	//sim_counts_matrix.resize(0, 0);

	//_Cs normalize_in_place(this->fake_expr_); // to do

	return true;

};

bool InferCnvWorker::log_transformation_3() {

	this->expr_ = log2(this->expr_.array() + 1.0).eval();

	if (this->hmm_) {
		this->fake_expr_ = log2(this->fake_expr_.array() + 1.0).eval();
	}

	if (this->scale_data_) {
		
		_Cs scale_in_place(this->expr_);

		if (this->hmm_) {
			_Cs scale_in_place(this->fake_expr_);
		}

	}

	return true;
};

bool InferCnvWorker::subtract_average_reference_5_8() {

	int n_ref_group = this->ref_grouped_cell_indices_.size();
	int n_gene = this->expr_.rows();
	int n_cell = this->expr_.cols();

	Eigen::MatrixXd gene_mean = Eigen::MatrixXd::Zero(n_gene, n_ref_group);

	int count{ 0 };
	for (auto&& [type, ind] : this->ref_grouped_cell_indices_) {
		gene_mean.col(count++) = this->expr_(Eigen::all, ind).rowwise().mean();
	}

	if (this->use_bounds_) {
		for (int i = 0; i < n_gene; ++i) {
			auto [min, max] = std::ranges::minmax(gene_mean.row(i));
			for (int j = 0; j < n_cell; ++j) {
				if (this->expr_(i, j) > max) {
					this->expr_(i, j) -= max;
				}
				else if (this->expr_(i, j) < min) {
					this->expr_(i, j) -= min;
				}
				else {
					this->expr_(i, j) = 0.0;
				}
			}

		}
	}
	else {
		this->expr_.colwise() -= gene_mean.rowwise().mean();
	}

	this->expr_ = this->expr_.unaryExpr([this](auto&& val)->double {
		if (val > this->max_centered_threshold_) {
			return this->max_centered_threshold_;
		}
		else if (val < -this->max_centered_threshold_) {
			return -this->max_centered_threshold_;
		}
		else {
			return val;
		}
	});

	if (!this->hmm_) {
		return true;
	}

	int n_fake_ref_group = this->fake_ref_grouped_cell_indices_.size();
	int n_fake_cell = this->fake_expr_.cols();

	Eigen::MatrixXd fake_gene_mean = Eigen::MatrixXd::Zero(n_gene, n_fake_ref_group);

	int fake_count{ 0 };
	for (auto&& [type, ind] : this->fake_ref_grouped_cell_indices_) {
		fake_gene_mean.col(fake_count++) = this->fake_expr_(Eigen::all, ind).rowwise().mean();
	}

	if (this->use_bounds_) {
		for (int i = 0; i < n_gene; ++i) {
			auto [min, max] = std::ranges::minmax(fake_gene_mean.row(i));
			for (int j = 0; j < n_cell; ++j) {
				if (this->fake_expr_(i, j) > max) {
					this->fake_expr_(i, j) -= max;
				}
				else if (this->fake_expr_(i, j) < min) {
					this->fake_expr_(i, j) -= min;
				}
				else {
					this->fake_expr_(i, j) = 0.0;
				}
			}

		}
	}
	else {
		this->fake_expr_.colwise() -= fake_gene_mean.rowwise().mean();
	}

	this->fake_expr_ = this->fake_expr_.unaryExpr([this](auto&& val)->double {
		if (val > this->max_centered_threshold_) {
			return this->max_centered_threshold_;
		}
		else if (val < -this->max_centered_threshold_) {
			return -this->max_centered_threshold_;
		}
		else {
			return val;
		}
	});

	return true;
};

Eigen::MatrixXd running_mean_mt(
	const Eigen::MatrixXd& mat,
	const int window_size
) {
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	int span = (window_size - 1) / 2;

	Eigen::ArrayXd index = Eigen::ArrayXd::LinSpaced(window_size, 1, window_size);
	Eigen::ArrayXd pyramid = index.cwiseMin(index.reverse());

	Eigen::MatrixXd res(nrow, ncol);

#pragma omp parallel for
	for (int i = 0; i < nrow; ++i) {
		int start = i - span;
		if (start < 0) {
			start = 0;
		}

		int ind_start = start - i + span;

		int end = i + span;
		if (end > nrow - 1) {
			end = nrow - 1;
		}

		int ind_end = end - i + span;

		res.row(i) = (mat.block(start, 0, end - start + 1, ncol).array().colwise() * pyramid.segment(ind_start, ind_end - ind_start + 1) / pyramid.segment(ind_start, ind_end - ind_start + 1).sum()).colwise().sum();
	}

	return res;
}

// use pyramidinal
bool InferCnvWorker::smooth_6() {

	auto chrs = _Cs unique(this->chrs_);

	for (auto&& chr : chrs) {

		auto ind = _Cs match(this->chrs_, chr);

		Eigen::MatrixXd m = this->expr_(ind, Eigen::all);
		this->expr_(ind, Eigen::all) = running_mean_mt(m, this->window_length_);
	}

	if (!this->hmm_) {
		return true;
	}

	auto fake_chrs = _Cs unique(this->fake_chrs_);

	for (auto&& chr : fake_chrs) {

		auto ind = _Cs match(this->fake_chrs_, chr);

		Eigen::MatrixXd m = this->fake_expr_(ind, Eigen::all);
		this->fake_expr_(ind, Eigen::all) = running_mean_mt(m, this->window_length_);
	}

	return true;
};

bool InferCnvWorker::center_7() {

	auto meds = _Cs sapply(_Cs seq_n(0, this->expr_.cols()), [this](auto&& i) {return _Cs median(this->expr_.col(i)); });

	this->expr_.array().rowwise() -= _Cs cast<Eigen::ArrayX>(meds).transpose();

	if (!this->hmm_) {
		return true;
	}

	meds = _Cs sapply(_Cs seq_n(0, this->fake_expr_.cols()), [this](auto&& i) {return _Cs median(this->fake_expr_.col(i)); });

	this->fake_expr_.array().rowwise() -= _Cs cast<Eigen::ArrayX>(meds).transpose();

	return true;
};

bool InferCnvWorker::invert_log_transform_10() {

	this->expr_ = this->expr_.unaryExpr([](auto&& val) {return std::pow(2.0, val); });

	if (this->hmm_) {

		this->fake_expr_ = this->fake_expr_.unaryExpr([](auto&& val) {return std::pow(2.0, val); });
	}

	return true;
};

bool InferCnvWorker::compute_tumor_subcluster_11() {

	Eigen::MatrixXd expr = this->expr_;
	QStringList chrs = this->chrs_;

	if (this->z_score_filter_ > 0.0) {

		auto ind = _Cs unroll(_Cs sapply(this->ref_grouped_cell_indices_, [](auto&& it) {return it.second; }));

		Eigen::MatrixXd ref_expr = expr(Eigen::all, ind);
		Eigen::MatrixXd z_score = (ref_expr.array() - ref_expr.mean()) / _Cs sd(ref_expr);
		Eigen::ArrayX<bool> filter = z_score.array().abs().rowwise().mean() > this->z_score_filter_;

		if (filter.count() > 0) {
			expr = _Cs row_sliced(expr, _Cs flip(filter));
			chrs = _Cs sliced(chrs, _Cs flip(filter));
		}
	}

	int n_cell = this->expr_.cols();

	this->subclusters_.resize(n_cell);

	for (auto&& [type, ind] : this->ref_grouped_cell_indices_) {

		Eigen::MatrixXd group_expr = expr(Eigen::all, ind);

		double leiden_resolution = std::pow(11.98 / group_expr.cols(), 1.0 / 1.165);

		auto pca = pca_infercnv(group_expr, _Cs min(2000, group_expr.rows() / 2), 10);
		if (pca.size() == 0) {
			G_TASK_WARN("PCA failed in subcluster.");
			return false;
		}

		auto partition = leiden_cluster(
			pca, "CPM", "Euclidean", 30, 50, leiden_resolution
		);
		_Cs assign(this->subclusters_, partition, ind);
	}

	for (auto&& [type, ind] : this->obs_grouped_cell_indices_) {

		Eigen::MatrixXd group_expr = expr(Eigen::all, ind);

		double leiden_resolution = std::pow(11.98 / group_expr.cols(), 1.0 / 1.165);

		auto pca = pca_infercnv(group_expr, _Cs min(2000, group_expr.rows() / 2), 10);
		if (pca.size() == 0) {
			G_TASK_WARN("PCA failed in subcluster.");
			return false;
		}
		auto partition = leiden_cluster(
			pca, "CPM", "Euclidean", 30, 50, leiden_resolution
		);
		_Cs assign(this->subclusters_, partition, ind);
	}

	return true;
};

bool InferCnvWorker::hmm_13() {

	auto fake_obs_ind = _Cs unroll(_Cs sapply(this->fake_obs_grouped_cell_indices_, [](auto&& it) {return it.second; }));

	Eigen::MatrixXd spike_expr = this->fake_expr_(Eigen::all, fake_obs_ind);

	QVector<double> spike_cnv;

	for (auto&& [chr, cnv, n_cell] : this->fake_chr_info_) {
		spike_cnv << QVector<double>(n_cell, cnv);
	}

	auto cnv_levels = _Cs unique(spike_cnv);

	std::map<double, std::pair<double, double>> cnv_mean_sd;

	int nrounds{ 100 };

	std::map<double, std::pair<double, double>> cnv_level_to_mean_sd_fit;

	for (auto&& level : cnv_levels) {
		Eigen::MatrixXd level_expr = this->fake_expr_(_Cs match(spike_cnv, level), Eigen::all);

		double mean = level_expr.mean();

		double sd = _Cs sd(level_expr);

		cnv_mean_sd[level] = { mean, sd };

		std::size_t s = level_expr.size() - 1;
		std::mt19937 re(this->random_state_);
		std::uniform_int_distribution<std::size_t> dis(0, s);

		const double* d1 = level_expr.data();
		
		Eigen::ArrayXd sds(100);

		for (int i = 2; i < 102; ++i) {
			Eigen::MatrixXd sample(i, nrounds);

			std::size_t n_sample = std::size_t(i) * nrounds;

			double* d2 = sample.data();

			for (std::size_t j = 0; j < n_sample; ++j) {
				d2[j] = d1[dis(re)];
			}

			Eigen::ArrayXd means = sample.rowwise().mean();

			sds[i - 1] = _Cs sd(means);
		}

		Eigen::ArrayX<bool> filter = sds.unaryExpr([](double d) {return !std::isnan(d) && d > 0.0; });

		if (filter.count() < 10) {
			G_TASK_WARN("meeting error in hmm cnv fitting.");
			return false;
		}

		sds = _Cs sliced(sds, filter);

		Eigen::ArrayXd n_cells = Eigen::ArrayXd::LinSpaced(100, 2, 101);

		n_cells = _Cs sliced(n_cells, filter);

		auto fit = lm(log(sds), log(n_cells));

		double slope = fit.coefficients[1];
		double intercept = fit.coefficients[0];
		if (std::isnan(slope) && std::isnan(intercept)) {
			G_TASK_WARN("meeting error in hmm cnv fitting.");
			return false;
		}

		cnv_level_to_mean_sd_fit[level] = { slope, intercept };
	}

	auto chrs = _Cs unique(this->chrs_);

	this->hmm_data_.resize(this->expr_.rows(), this->expr_.cols());

	auto tumor_subclusters = _Cs unique(this->subclusters_);

	Eigen::ArrayXd delta = Eigen::ArrayXd::Constant(6, 1e-6);
	delta(2) = 1.0 - 5e-6;

	Eigen::MatrixXd state_transitions(6, 6);

	double t = 1e-6;

	state_transitions <<
		1 - 5 * t, t, t, t, t, t,
		t, 1 - 5 * t, t, t, t, t,
		t, t, 1 - 5 * t, t, t, t,
		t, t, t, 1 - 5 * t, t, t,
		t, t, t, t, 1 - 5 * t, t,
		t, t, t, t, t, 1 - 5 * t;

	Eigen::MatrixXd logpi = log(state_transitions.array());

	Eigen::ArrayXd mean(6);
	mean[0] = cnv_mean_sd[0.01].first;
	mean[1] = cnv_mean_sd[0.5].first;
	mean[2] = cnv_mean_sd[1].first;
	mean[3] = cnv_mean_sd[1.5].first;
	mean[4] = cnv_mean_sd[2].first;
	mean[5] = cnv_mean_sd[3].first;

	for (auto&& chr : chrs) {

		auto chr_ind = _Cs match(this->chrs_, chr);

		if (chr_ind.size() == 1) {
			this->hmm_data_.row(chr_ind[0]).array() = 3;
			continue;
		}

		for (auto&& cluster : tumor_subclusters) {

			auto subcluster_ind = _Cs match(this->subclusters_, cluster);

			Eigen::ArrayXd gene_expr_vals = this->expr_(chr_ind, subcluster_ind).rowwise().mean();

			int num_cells = subcluster_ind.size();


			std::map<double, std::pair<double, double>> cluster_cnv_mean_sd = cnv_mean_sd;
			for (auto&& [cnv, ms] : cluster_cnv_mean_sd) {
				auto [slope, intercept] = cnv_level_to_mean_sd_fit[cnv];
				ms.second = slope * ms.first + intercept;
			}

			int n = gene_expr_vals.size();
			int m{ 6 };
			Eigen::MatrixXd nu(n, 6);
			Eigen::ArrayXi y(n);
			double pseudo_count{ 1e-20 };
			Eigen::MatrixXd emissions(n, 6);

			double sd = _Cs median(_Cs sapply(cluster_cnv_mean_sd, [](auto&& it) {return it.second.second; }));

			Eigen::ArrayXd emission = log(p_normal(Eigen::ArrayXd((mean - gene_expr_vals[0]).abs() / sd), 0.0, 1.0, false));
			emission = 1.0 / (-emission);
			emission /= emission.sum();
			emissions.row(0) = log(emission);

			nu.row(0) = log(delta) + emissions.row(0).array();

			for (int i = 1; i < n; ++i) {
				Eigen::MatrixXd matrixnu(6, 6);

				matrixnu.colwise() = nu.row(i - 1).transpose();

				emission = log(p_normal(Eigen::ArrayXd((mean - gene_expr_vals[i]).abs() / sd), 0.0, 1.0, false));
				emission = 1.0 / (-emission);
				emission /= emission.sum();
				emissions.row(i) = log(emission);

				nu.row(i) = (matrixnu.array() + logpi.array()).colwise().maxCoeff() + emissions.row(i).array();
			}

			if (_Cs any(_Cs sapply(nu.row(n - 1), [](auto&& d) {return std::isnan(d) || std::isinf(d); }), true)) {
				G_TASK_WARN("Underflow");
				return false;
			}

			y(n - 1) = _Cs argmax(nu.row(n - 1));

			for (int i = n - 2; i > -1; --i) {
				y(i) = _Cs argmax(logpi.col(y(i + 1)) + nu.row(i).transpose());
			}

			this->hmm_data_(chr_ind, subcluster_ind).array().colwise() = y;
		}
	}

	return true;
};

bool InferCnvWorker::generate_reports_14() {

	/*auto subclusters = _Cs unique(this->subclusters_);
	auto chrs = _Cs unique(this->chrs_);

	int n_gene = this->hmm_data_.rows();

	QVector<int> cnv_id;
	QVector<int> state;
	QStringList gene;
	QStringList chromosome;
	QVector<int> start;
	QVector<int> end;

	int id{ 0 };
	int sstate{ -1 };
	for (auto&& cluster : subclusters) {
		Eigen::MatrixXi cell_group_mtx = this->hmm_data_(Eigen::all, _Cs match(this->subclusters_, cluster));
		auto state_consensus = _Cs sapply(_Cs seq_n(0, n_gene),
			[this](int i) {return _Cs find_most_frequent_element(this->hmm_data_.row(i)); });

		for (auto&& chr : chrs) {

			auto gene_idx = _Cs match(this->chrs_, chr);
			auto chr_state = _Cs reordered(state_consensus, gene_idx);


		}
	}*/

	return true;
};
//
//double logistic_midpt_slope(double input, const Eigen::ArrayXd& params)
//{
//
//	return 1.0 / (1.0 + std::exp(-params[1] * (input - params[0])));
//};
//
//Eigen::MatrixXd get_simulated_cell_matrix_using_meanvar_trend(
//	const SparseDouble& expr,
//	const Eigen::ArrayXd& gene_means,
//	int n_cell,
//	const std::map<QString, QVector<int>>& ref_grouped_cell_indices,
//	const std::map<QString, QVector<int>>& obs_grouped_cell_indices) 
//{
//
//	QVector<double> mean, p0, var;
//
//	for (auto&& [group, ind] : ref_grouped_cell_indices) {
//
//		auto expr_sub = expr.col_reordered(ind);
//
//		mean << _Cs cast<QVector>(_Cs row_mean(expr_sub.mat_));
//		p0 << _Cs cast<QVector>(1.0 - _Cs row_count<double, true, false>(expr_sub.mat_, 0).cast<double>() / ind.size());
//		var << _Cs cast<QVector>(_Cs row_var(expr_sub.mat_));
//	}
//
//	for (auto&& [group, ind] : obs_grouped_cell_indices) {
//
//		auto expr_sub = expr.col_reordered(ind);
//
//		mean << _Cs cast<QVector>(_Cs row_mean(expr_sub.mat_));
//		p0 << _Cs cast<QVector>(1.0 - _Cs row_count<double, true, false>(expr_sub.mat_, 0).cast<double>() / ind.size());
//		var << _Cs cast<QVector>(_Cs row_var(expr_sub.mat_));
//	}
//
//	auto filter = _Cs greater_than(mean, 0.0);
//	Eigen::ArrayXd x = log(_Cs cast<Eigen::ArrayX>(_Cs sliced(mean, filter)));
//	Eigen::ArrayXd y = _Cs cast<Eigen::ArrayX>(_Cs sliced(p0, filter));
//
//	Eigen::ArrayXd params = Eigen::ArrayXd::Zero(2);
//	params[0] = x.mean();
//	params[1] = -1.0;
//	
//	nls_LevenbergMarquardt(logistic_midpt_slope, x, y, params);
//
//	double midpt = params[0];
//	double slope = params[1];
//
//	const int ngenes = gene_means.size();
//
//	Eigen::ArrayXd logm = log(_Cs cast<Eigen::ArrayX>(mean) + 1.0);
//	Eigen::ArrayXd logv = log(_Cs cast<Eigen::ArrayX>(var) + 1.0);
//
//	Eigen::MatrixXd sim_cell_matrix = Eigen::MatrixXd::Zero(ngenes, n_cell);
//
//
//};