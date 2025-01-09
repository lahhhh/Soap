#include "CiceroWorker.h"

#include <QFile>
#include <QTextStream>
#include "GenomeUtility.h"

#include "Custom.h"
#include "glasso.h"

bool CiceroWorker::work() {

	this->aggregate();

	if (!this->estimate_distance_parameter()) {
		G_TASK_WARN("No distance parameters calculated.");
		return false;
	}

	try {
		if (!this->generate_cicero_models()) {
			G_TASK_WARN("Models generation failed.");
			return false;
		}
	}
	catch (...) {
		G_TASK_WARN("Meeting problems in generating models");
		return false;
	}

	if (!this->assemble_connections()) {
		return false;
	}

	if (!this->generate_ccans()) {
		G_TASK_WARN("Found no ccan.");
		return false;
	}

	this->create_ccan_matrix();

	this->res_.reset(new Cicero());

	int n_group = this->regulation_group_counts_.rows();
	int n_cell = this->atac_counts_->cols();

	this->res_->connections_ = this->connections_;
	this->res_->regulation_groups_ = this->regulation_groups_;

	this->res_->regulation_group_counts_.colnames_ = this->atac_counts_->colnames_;
	this->res_->regulation_group_counts_.mat_ = this->regulation_group_counts_;
	this->res_->regulation_group_counts_.rownames_ = custom::cast<QString>(custom::seq_n(1, n_group));

	this->res_->regulation_group_normalized_.colnames_ = this->res_->regulation_group_counts_.colnames_;
	this->res_->regulation_group_normalized_.rownames_ = this->res_->regulation_group_counts_.rownames_;
	this->res_->regulation_group_normalized_.mat_ = this->res_->regulation_group_counts_.mat_.cast<double>();

	double mean = this->res_->regulation_group_normalized_.mat_.sum() / n_cell;
	Eigen::ArrayXd depth = this->res_->regulation_group_normalized_.mat_.colwise().sum();
	for (int i = 0; i < n_cell; ++i) {
		if (depth[i] != 0.0) {
			this->res_->regulation_group_normalized_.mat_.col(i) *= mean / depth[i];
		}
	}

	this->res_->regulation_group_normalized_.mat_ = log(this->res_->regulation_group_normalized_.mat_.array() + 1.0).eval();

	this->res_->peak_names_ = this->atac_counts_->rownames_;

	return true;
};

void CiceroWorker::run() {

	G_TASK_LOG("Start cicero...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_cicero_ready(this->res_.release());

	G_TASK_LOG("Cicero finished.");

	G_TASK_END;
};

bool CiceroWorker::aggregate() {

	auto knn_ind = custom::get_knn_mt(this->embedding_->data_.mat_, this->k_);

	int n_cell = knn_ind.rows();

	auto good_choices = custom::seq_n(0, n_cell);

	std::default_random_engine e(1997);

	std::shuffle(good_choices.begin(), good_choices.end(), e);

	int it{ 0 };
	int choose_ind = 0;

	QVector<int> chosen, new_chosen;
	chosen << good_choices[choose_ind];

	++choose_ind;

	while (choose_ind != n_cell && it < 5000) {

		++it;
		int choice = good_choices[choose_ind];
		new_chosen = chosen;
		new_chosen << choice;

		++choose_ind;

		Eigen::MatrixXi cell_sample = knn_ind(new_chosen, Eigen::all);

		int n_chosen = new_chosen.size();
		QVector<int> shared(n_chosen - 1);
		for (int i = 0; i < n_chosen - 1; ++i) {
			Eigen::ArrayXi tmp(2 * this->k_);
			tmp.segment(0, this->k_) = cell_sample.row(i);
			tmp.segment(this->k_, this->k_) = cell_sample.row(n_chosen - 1);
			shared[i] = 2 * this->k_ - (custom::unique(tmp)).size();
		}

		if (std::ranges::max(shared) < 0.9 * this->k_) {
			chosen = new_chosen;
		}
	}

	int n_chosen = new_chosen.size();
	Eigen::MatrixXi cell_sample = knn_ind(new_chosen, Eigen::all);
	int n_feature = this->atac_counts_->rows();

	this->cicero_counts_.resize(n_feature, n_chosen);
	for (int i = 0; i < n_chosen; ++i) {
		this->cicero_counts_.col(i) = custom::column_reorder_and_row_sum(this->atac_counts_->mat_, cell_sample.row(i));
	}

	return true;
};

void CiceroWorker::create_ccan_matrix() {

	int nrow = this->regulation_groups_.size();
	int n_cell = this->atac_counts_->cols();

	this->regulation_group_counts_.resize(nrow, n_cell);
	this->regulation_group_counts_.setZero();

#pragma omp parallel for
	for (int i = 0; i < nrow; ++i) {
		for (auto ind : this->regulation_groups_[i]) {
			this->regulation_group_counts_.row(i) += this->atac_counts_->mat_.row(ind);
		}
	}
};

double CiceroWorker::find_ccan_cutoff() {
	constexpr int tolerance_digits = 2;
	constexpr double tolerance = 0.01;
	double bottom{ 0.0 };
	double top{ 1.0 };

	while ((top - bottom) > tolerance) {
		double test_val = (top + bottom) / 2;
		auto membership = custom::cluster_louvain_igraph(this->connections_, test_val);
		auto unique_cluster = custom::unique(membership);

		int n_cluster{ 0 }, n_cluster2{ 0 };
		int n_vertice = membership.size();

		for (auto c : unique_cluster) {

			int cluster_size{ 0 };
			for (int i = 0; i < n_vertice; ++i) {
				if (membership[i] == c) {
					++cluster_size;
				}
			}

			if (cluster_size > 2) {
				++n_cluster;
			}
		}

		double next_step = test_val;
		while (true) {
			next_step += (top - bottom) / 10;
			membership = custom::cluster_louvain_igraph(this->connections_, next_step);
			unique_cluster = custom::unique(membership);

			n_cluster2 = 0;
			for (auto c : unique_cluster) {

				int cluster_size{ 0 };
				for (int i = 0; i < n_vertice; ++i) {
					if (membership[i] == c) {
						++cluster_size;
					}
				}

				if (cluster_size > 2) {
					++n_cluster2;
				}
			}

			if (n_cluster2 != n_cluster) {
				break;
			}
		}

		if (n_cluster > n_cluster2) {
			top = test_val;
		}
		else {
			bottom = test_val;
		}
	}

	return (top + bottom) / 2;
};

bool CiceroWorker::generate_ccans() {

	this->coaccess_cutoff_ = this->find_ccan_cutoff();

	Eigen::ArrayXi cluster = custom::cluster_louvain_igraph(this->connections_, this->coaccess_cutoff_);

	int n_vertice = cluster.size();

	auto unique_cluster = custom::unique(cluster);

	QVector<int> peak_index;

	for (auto c : unique_cluster) {
		peak_index.clear();

		for (int i = 0; i < n_vertice; ++i) {
			if (cluster[i] == c) {
				peak_index << i;
			}
		}

		if (peak_index.size() > 2) {
			this->regulation_groups_ << peak_index;
		}
	}

	if (this->regulation_groups_.isEmpty()) {
		return false;
	}

	return true;
};

bool CiceroWorker::assemble_connections() {

	std::vector<Eigen::Triplet<double>> triplets;

	for (auto&& [ind, w] : this->models) {

		auto cor = custom::cov2cor(w);

		int n_p = ind.size();

		for (int i = 0; i < n_p; ++i) {
			for (int j = 0; j < n_p; ++j) {

				if (i == j || cor(i, j) == 0.0 ) {
					continue;
				}

				triplets.emplace_back(ind[i], ind[j], cor(i, j));
			}
		}
	}

	int n_peak = this->peak_chr_names_.size();
	this->connections_ = custom::set_from_triplets_mean(triplets, n_peak, n_peak);

	return true;
};

bool CiceroWorker::generate_cicero_models() {

	int n_window = this->windows_.size();
	int n_peak = this->peak_chr_names_.size();
	int n_agg = this->cicero_counts_.cols();

	double distance_parameter = custom::mean(this->distance_parameters_);

#pragma omp parallel for
	for (int i = 0; i < n_window; ++i) {

		auto [chr_name, start, end, strand] = this->windows_.at(i);

		QVector<int> peak_index;
		QVector<double> mean_bp;

		for (int i = 0; i < n_peak; ++i) {

			if (this->peak_chr_names_[i] != chr_name) {
				continue;
			}

			int s = this->peak_starts_[i], e = this->peak_ends_[i];

			if ((s > start && s < end) || (e > start && e < end)) {
				peak_index << i;
				mean_bp << (s + e) / 2;
			}
		}



		int n_peak_in = peak_index.size();

		if (n_peak_in <= 1 || n_peak_in > this->max_element_) {
			continue;
		}

		auto dist_matrix = custom::distance(mean_bp);
		Eigen::MatrixXd vals(n_peak_in, n_agg);
		for (int i = 0; i < n_peak_in; ++i) {
			vals.row(i) = this->cicero_counts_.row(peak_index[i]).cast<double>();
		}

		auto cov_mat = custom::cov(vals.transpose());
		cov_mat.diagonal().array() += 1e-4;

		constexpr double s{ 0.75 };

		Eigen::MatrixXd rho = (1 - pow((1000.0 / dist_matrix.array()), s)) * distance_parameter;
		for (int i = 0; i < n_peak_in; ++i) {
			for (int j = 0; j < n_peak_in; ++j) {
				double val = rho(i, j);
				if (std::isinf(val) || std::isnan(val) || val < 0) {
					rho(i, j) = 0.0;
				}
			}
		}

		auto [w, wi] = glasso(cov_mat, rho);

		bool stop{ false };
		for (int i = 0; i < n_peak_in; ++i) {
			for (int j = 0; j < n_peak_in; ++j) {
				double val = w(i, j);
				if (std::isnan(val) || std::isinf(val)) {
					stop = true;
					break;
				}
			}

			if (stop) {
				break;
			}
		}

		if (stop) {
			continue;
		}

	#pragma omp critical
		{
			this->models.emplace_back(peak_index, w);
		}
	}

	if (this->models.empty()) {
		return false;
	}

	return true;
};

bool CiceroWorker::estimate_distance_parameter() {

	QFile file(FILE_CICERO_HG38_CHROMSIZE);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		G_TASK_WARN("Can not open resource");
		return false;
	}

	QTextStream in(&file);

	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList content = line.split('\t');
		QString chr_name = custom::standardize_chromosome_name(content[0]);
		int chr_size = content[1].toInt();

		if (chr_size <= 1) {
			G_TASK_WARN("Broken Resource File.");
			return false;
		}

		int start = 1;

		while (start < chr_size) {
			this->windows_.append(chr_name, start, this->window_, '*');

			start += this->window_ / 2;
		}
	}

	this->windows_.finalize();
	int n_window = this->windows_.size();

	int distance_parameter_calculated{ 0 };
	int it{ 0 };

	int n_peak = this->atac_counts_->rows();
	int n_agg = this->cicero_counts_.cols();

	for (int i = 0; i < n_peak; ++i) {
		auto [chr_name, start, end, success] = custom::string_to_peak(this->atac_counts_->rownames_[i]);

		if (!success) {
			G_TASK_WARN("Illegal Peak Name : " + this->atac_counts_->rownames_[i]);
			return false;
		}

		this->peak_chr_names_ << chr_name;
		this->peak_starts_ << start;
		this->peak_ends_ << end;
	}

	std::default_random_engine e;
	std::uniform_int_distribution u(0, n_window - 1);

	while (this->sample_num_ > distance_parameter_calculated && it < this->max_sample_windows_) {

		QVector<int> samples;
		for (int ii = 0; ii < 20; ++ii) {
			samples << u(e);
		}

		it += 20;

	#pragma omp parallel for
		for (int ii = 0; ii < 20; ++ii) {

			int sampled = samples[ii];

			auto [chr_name, start, end, strand] = this->windows_.at(sampled);

			QVector<int> peak_index;
			QVector<double> mean_bp;

			for (int i = 0; i < n_peak; ++i) {

				if (this->peak_chr_names_[i] != chr_name) {
					continue;
				}

				int s = this->peak_starts_[i], e = this->peak_ends_[i];

				if ((s > start && s < end) || (e > start && e < end)) {
					peak_index << i;
					mean_bp << (s + e) / 2;
				}
			}

			int n_peak_in = peak_index.size();

			if (n_peak_in <= 1 || n_peak_in > this->max_element_) {
				continue;
			}

			auto dist_matrix = custom::distance(mean_bp);

			if ((dist_matrix.array() > (double)distance_constraint_).count() < 2) {
				continue;
			}

			double starting_max{ 2.0 };
			double distance_parameter{ 2.0 };
			double distance_parameter_max{ 2.0 };
			double distance_parameter_min{ 0.0 };
			int it2{ 0 };
			int maxit{ 100 };

			bool found{ false };
			while (!found && it2 < maxit) {

				Eigen::MatrixXd vals(n_peak_in, n_agg);
				for (int i = 0; i < n_peak_in; ++i) {
					vals.row(i) = this->cicero_counts_.row(peak_index[i]).cast<double>();
				}

				auto cov_mat = custom::cov(vals.transpose());
				cov_mat.diagonal().array() += 1e-4;

				constexpr double s{ 0.75 };

				Eigen::MatrixXd rho = (1 - pow((1000.0 / dist_matrix.array()), s)) * distance_parameter;
				for (int i = 0; i < n_peak_in; ++i) {
					for (int j = 0; j < n_peak_in; ++j) {
						double val = rho(i, j);
						if (std::isinf(val) || std::isnan(val) || val < 0) {
							rho(i, j) = 0.0;
						}
					}
				}

				auto [w, wi] = glasso(cov_mat, rho);


				int big_entries = (dist_matrix.array() > this->distance_constraint_).count();

				double e1{ 0.0 }, e2{ 0.0 };
				for (int i = 0; i < n_peak_in; ++i) {
					for (int j = 0; j < n_peak_in; ++j) {

						if (wi(i, j) == 0.0) {
							e2 += 1.0;
						}
						else if (dist_matrix(i, j) > this->distance_constraint_) {
							e1 += 1.0;
						}
					}
				}

				bool longs_zero{ true };
				if ((e1 / big_entries > 0.05) || (e2 / wi.size() < 0.2)) {
					longs_zero = false;
				}

				if (!longs_zero || (distance_parameter == 0.0)) {
					distance_parameter_min = distance_parameter;
				}
				else {
					distance_parameter_max = distance_parameter;
				}

				double new_distance_parameter = (distance_parameter_max + distance_parameter_min) / 2;
				if (new_distance_parameter == starting_max) {
					new_distance_parameter = 2 * starting_max;
					starting_max = new_distance_parameter;
				}

				if (distance_parameter_convergence_ > std::abs(distance_parameter - new_distance_parameter)) {
					found = true;
				}
				else {
					distance_parameter = new_distance_parameter;
				}

				++it2;
			}

		#pragma omp critical
			{
				distance_parameters_ << distance_parameter;
				++distance_parameter_calculated;
			}
		}
	}

	if (this->distance_parameters_.isEmpty()) {
		return false;
	}

	return true;
};