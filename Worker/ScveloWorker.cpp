#include "ScveloWorker.h"

#include "annoylib.h"
#include "kissrandom.h"
#include "mman.h"
#include "Custom.h"

#include "lm.h"

#define SMOOTH_K_TOLERANCE 1e-5
#define MIN_K_DIST_SCALE 1e-3

bool ScveloWorker::work() {

	if (!this->get_velocity()) {
		return false;
	}

	this->velocity_graph();

	return true;
};

void ScveloWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_scvelo_ready(this->res_.release());

	G_TASK_END;
}


void ScveloWorker::fuzzy_simplicial_set() {

	constexpr int n_trees = 50;

	if (this->metric_ == "Euclidean") {
		std::tie(this->knn_indices_, this->knn_distance_) = custom::get_knn_mt<Euclidean, true>(this->mat_, this->n_neighbors_, n_trees);
	}
	else if (this->metric_ == "Angular") {
		std::tie(this->knn_indices_, this->knn_distance_) = custom::get_knn_mt<Angular, true>(this->mat_, this->n_neighbors_, n_trees);
	}
	else {
		std::tie(this->knn_indices_, this->knn_distance_) = custom::get_knn_mt<Manhattan, true>(this->mat_, this->n_neighbors_, n_trees);
	}

	auto [sigmas, rho] = this->smooth_knn_distance();

	this->compute_membership_strengths(sigmas, rho);
};

std::pair<Eigen::ArrayXd, Eigen::ArrayXd> ScveloWorker::smooth_knn_distance(double local_connectivity, int n_iter) { // bandwidth --> 1.0
	int sample_number = this->mat_.rows();
	double target = log2((double)this->n_neighbors_);
	Eigen::ArrayXd rho = Eigen::ArrayXd::Zero(sample_number), result = rho;
	double mean_distances = this->knn_distance_.mean();
	int index = (int)floor(local_connectivity);
	double interpolation = local_connectivity - index;

	for (int i = 0; i < sample_number; ++i) {
		double lo = 0.0, mid = 1.0, hi = std::numeric_limits<double>::max();
		double cnt = this->knn_distance_.row(i).count();

		Eigen::ArrayXd ith_distances = this->knn_distance_.row(i);
		Eigen::ArrayXd non_zero_dists = custom::sliced(ith_distances, ith_distances > 0);

		if (cnt > local_connectivity) {
			if (index > 0) {
				rho[i] = non_zero_dists[index - 1];
				if (interpolation > SMOOTH_K_TOLERANCE) {
					rho[i] += interpolation * (non_zero_dists[index] - non_zero_dists[index - 1]);
				}
			}
			else {
				rho[i] = interpolation * non_zero_dists[0];
			}
		}
		else if (cnt > 0) {
			rho[i] = non_zero_dists.maxCoeff();
		}

		for (int j = 0; j < n_iter; ++j) {
			double psum = 0.0;
			for (int n = 1; n < this->n_neighbors_; n++) {
				double d = this->knn_distance_(i, n) - rho[i];
				if (d > 0)
					psum += exp(-(d / mid));
				else
					psum += 1.0;
			}
			if (std::abs(psum - target) < SMOOTH_K_TOLERANCE)
				break;

			if (psum > target) {
				hi = mid;
				mid = (lo + hi) / 2.0;
			}
			else {
				lo = mid;
				if (hi == std::numeric_limits<double>::max())
					mid *= 2;
				else
					mid = (lo + hi) / 2.0;
			}
		}
		result[i] = mid;
		if (rho[i] > 0.0) {
			double mean_ith_distances = ith_distances.mean();
			if (result[i] < MIN_K_DIST_SCALE * mean_ith_distances)
				result[i] = MIN_K_DIST_SCALE * mean_ith_distances;
		}
		else {
			if (result[i] < MIN_K_DIST_SCALE * mean_distances)
				result[i] = MIN_K_DIST_SCALE * mean_distances;
		}
	}

	return std::make_pair(result, rho);
};

bool ScveloWorker::normalize_counts() {
	auto spliced_counts = this->velocyto_base_->get_spliced();
	auto unspliced_counts = this->velocyto_base_->get_unspliced();

	if (spliced_counts == nullptr || unspliced_counts == nullptr) {
		G_TASK_WARN("Incomplete Velocyto Base item.");
		return false;
	}

	if (spliced_counts->rows() != unspliced_counts->rows() ||
		spliced_counts->cols() != unspliced_counts->cols() ||
		spliced_counts->is_empty() ||
		unspliced_counts->is_empty()
		)
	{
		G_TASK_WARN("Incomplete Velocyto Base item.");
		return false;
	}

	this->spliced_normalized_ = spliced_counts->mat_.cast<double>();
	this->unspliced_normalized_ = unspliced_counts->mat_.cast<double>();

	Eigen::ArrayXd cell_size_s = custom::col_sum(this->spliced_normalized_).cwiseMax(1.0) / 10000;
	Eigen::ArrayXd cell_size_u = custom::col_sum(this->unspliced_normalized_).cwiseMax(1.0) / 10000;

	custom::rowwise_divide_in_place(this->spliced_normalized_, cell_size_s);
	custom::rowwise_divide_in_place(this->unspliced_normalized_, cell_size_u);

	return true;

};

bool ScveloWorker::get_velocity() {
	if (!this->moments()) {
		return false;
	}

	if(!this->compute_deterministic()){
		return false;
	}

	this->compute_stochastic();

	return true;
};

void ScveloWorker::compute_stochastic() {

	const int n_gene = this->spliced_normalized_.rows(), n_cell = this->spliced_normalized_.cols();

	Eigen::MatrixXd var_ss(n_cell, n_gene), cov_us(n_cell, n_gene);

#pragma omp parallel for
	for(int i = 0; i < n_cell; ++i){
		Eigen::ArrayXd row = Eigen::ArrayXd::Zero(n_gene);

		for (auto ind : this->connectivities_[i]) {
			Eigen::ArrayXd col = this->spliced_normalized_.col(ind);
			row += col.square();
		}

		row /= this->connectivities_[i].size();

		var_ss.row(i) = 2 * row - this->ms_.row(i).array();
	};

#pragma omp parallel for
	for (int i = 0; i < n_cell; ++i) {
		Eigen::ArrayXd row = Eigen::ArrayXd::Zero(n_gene);

		for (auto ind : this->connectivities_[i]) {
			Eigen::ArrayXd col1 = this->unspliced_normalized_.col(ind);
			Eigen::ArrayXd col2 = this->spliced_normalized_.col(ind);
			row += (col1 * col2);
		}

		row /= this->connectivities_[i].size();

		cov_us.row(i) = 2 * row + this->mu_.row(i).array();
	};

	this->spliced_normalized_.resize(0, 0);
	this->unspliced_normalized_.resize(0, 0);

	Eigen::ArrayXd gamma2 = Eigen::ArrayXd::Zero(n_gene);

#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {
		Eigen::ArrayXd s_array = var_ss.col(i);
		Eigen::ArrayXd u_array = cov_us.col(i);
		double s_max = s_array.maxCoeff();
		double u_max = u_array.maxCoeff();

		if (s_max == 0 || u_max == 0) {
			continue;
		}

		auto fit_res = lm(u_array, s_array, false);

		gamma2(i) = fit_res.coefficients[0];
	};

	Eigen::ArrayXd res_std(n_gene), res_std2(n_gene);

#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {
		res_std(i) = custom::sd(this->residual_.col(i));

		Eigen::ArrayXd col = cov_us.col(i).array() - gamma2(i) * var_ss.col(i).array();
		res_std2(i) = custom::sd(col);
	};

#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {

		Eigen::ArrayXd x(2 * n_cell), y(2 * n_cell);
		x.segment(0, n_cell) = this->ms_.col(i).array() / res_std(i);
		x.segment(n_cell, n_cell) = var_ss.col(i).array() / res_std2(i);
		y.segment(0, n_cell) = this->mu_.col(i).array() / res_std(i);
		y.segment(n_cell, n_cell) = cov_us.col(i).array() / res_std2(i);
		
		auto fit_res = lm(y, x, false);

		double coef = fit_res.coefficients[0];

		gamma2(i) = std::isnan(coef) ? 0 : coef;
	};

	this->residual_ = this->mu_.array() - this->ms_.array().rowwise() * gamma2.transpose();

	this->mu_.resize(0, 0);

	//this->residual2_ = (cov_us.array() - 2 * this->ms_.array() * this->mu_.array()) - (var_ss.array() - 2 * this->ms_.array().square()).rowwise() * gamma2.transpose();
};

bool ScveloWorker::compute_deterministic() {
	
	const int n_gene = this->ms_.cols(), n_cell = this->ms_.rows();

	Eigen::ArrayXd gamma = Eigen::ArrayXd::Zero(n_gene);

#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {
		Eigen::ArrayXd s_array = this->ms_.col(i);
		Eigen::ArrayXd u_array = this->mu_.col(i);
		double s_max = s_array.maxCoeff();
		double u_max = u_array.maxCoeff();

		if (s_max == 0 || u_max == 0) {
			continue;
		}

		Eigen::ArrayXd normalized = s_array / s_max + u_array / u_max;
		auto sorted = custom::sorted(normalized);

		double p5 = custom::linear_percentile_no_sort(sorted, 5);
		double p95 = custom::linear_percentile_no_sort(sorted, 95);

		auto filter = (normalized <= p5) + (normalized >= p95);
		s_array = custom::sliced(s_array, filter);
		u_array = custom::sliced(u_array, filter);

		auto fit_res = lm(u_array, s_array, false);

		gamma(i) = fit_res.coefficients[0];
	};

	Eigen::ArrayXd r2 = Eigen::ArrayXd::Zero(n_gene);

	this->residual_ = this->mu_.array() - this->ms_.array().rowwise() * gamma.transpose();

	for (int i = 0; i < n_gene; ++i) {
		double col_mean = this->mu_.col(i).mean();
		double total = (this->mu_.col(i).array() - col_mean).square().sum();

		if (total == 0) {
			continue;
		}

		double residual = this->residual_.col(i).array().square().sum();

		r2(i) = 1 - residual / total;

	}

	constexpr double min_r2 = 0.01, min_ratio = 0.01;

	Eigen::ArrayX<bool> ms_filter = this->ms_.array().colwise().maxCoeff() > 0;
	Eigen::ArrayX<bool> mu_filter = this->mu_.array().colwise().maxCoeff() > 0;
	Eigen::ArrayX<bool> r2_filter = r2 > min_r2;
	Eigen::ArrayX<bool> ratio_filter = gamma > min_ratio;

	Eigen::ArrayX<bool> filter = ms_filter * mu_filter * r2_filter * ratio_filter;

	if (filter.count() < 100) {
		return false;
	}

	this->ms_ = custom::col_sliced(this->ms_, filter);
	this->mu_ = custom::col_sliced(this->mu_, filter);
	this->residual_ = custom::col_sliced(this->residual_, filter);

	this->spliced_normalized_ = custom::row_sliced(this->spliced_normalized_, filter);
	this->unspliced_normalized_ = custom::row_sliced(this->unspliced_normalized_, filter);

	return true;
};

bool ScveloWorker::moments() {

	if (!this->normalize_counts()) {
		return false;
	}

	this->fuzzy_simplicial_set();
	this->get_connectivities();

	// cell pooling
	const int n_cell = this->spliced_normalized_.cols(), n_gene = this->spliced_normalized_.rows();

	this->ms_.resize(n_cell, n_gene);

#pragma omp parallel for
	for (int i = 0; i < n_cell; ++i) {
		Eigen::ArrayXd row = Eigen::ArrayXd::Zero(n_gene);

		for (auto ind : this->connectivities_[i]) {
			row += this->spliced_normalized_.col(ind);
		}

		this->ms_.row(i) = row / this->connectivities_[i].size();
	}

	this->mu_.resize(n_cell, n_gene);

#pragma omp parallel for
	for (int i = 0; i < n_cell; ++i) {
		Eigen::ArrayXd row = Eigen::ArrayXd::Zero(n_gene);

		for (auto ind : this->connectivities_[i]) {
			row += this->unspliced_normalized_.col(ind);
		}

		this->mu_.row(i) = row / this->connectivities_[i].size();
	}

	return true;
};

void ScveloWorker::get_connectivities() {

	const auto size = this->strengths_.nonZeros();
	const auto n_cell = this->strengths_.rows();

	this->connectivities_.resize(n_cell);

#pragma omp parallel for
	for (int i = 0; i < n_cell; ++i) {
		Eigen::ArrayXd c = this->strengths_.col(i);

		auto ind = custom::which(c > 0.0);

		if (ind.size() < this->n_neighbors_) {
			this->connectivities_[i].append_range(ind);
		}
		else {
			Eigen::ArrayXd d = c(ind);
			auto order = custom::order(d, true);
			ind = custom::reordered(ind, order.segment(0, this->n_neighbors_));
			this->connectivities_[i].append_range(ind);
		}

	}

	this->strengths_.resize(0, 0);
};

void ScveloWorker::compute_membership_strengths(const Eigen::ArrayXd& sigmas, const Eigen::ArrayXd& rho) {

	const int n_samples = this->knn_distance_.rows();
	std::vector<Eigen::Triplet<double> > strength_triplets;

	for (int i = 0; i < n_samples; ++i) {
		for (int j = 0; j < this->n_neighbors_; ++j) {
			double val;
			if (this->knn_indices_(i, j) == i) {
				continue;
			}
			else if (this->knn_distance_(i, j) - rho[i] <= 0.0 || sigmas[i] == 0.0) {
				val = 1.0;
			}
			else {
				val = exp(-((this->knn_distance_(i, j) - rho[i]) / sigmas[i]));
			}
			strength_triplets.emplace_back(i, this->knn_indices_(i, j), val);
		}
	}

	this->strengths_.resize(n_samples, n_samples);
	this->strengths_.setFromTriplets(strength_triplets.cbegin(), strength_triplets.cend());

	Eigen::SparseMatrix<double> transpose = this->strengths_.transpose();
	Eigen::SparseMatrix<double> prod_matrix = this->strengths_.cwiseProduct(transpose);

	constexpr double set_op_mix_ratio = 1.0;
	this->strengths_ = /*set_op_mix_ratio * */ (this->strengths_ + transpose - prod_matrix).eval() /* + (1.0 - set_op_mix_ratio) * prod_matrix */;
	this->strengths_ = custom::eliminate_zero(this->strengths_);
};

void ScveloWorker::velocity_graph() {

	Eigen::MatrixXd& V = this->residual_;
	const int n_gene = V.cols(), n_cell = V.rows();
	for (int j = 0; j < n_gene; ++j) {
		for (int i = 0; i < n_cell; ++i) {
			if (V(i, j) > 0) {
				V(i, j) = std::sqrt(V(i, j));
			}
			else {
				V(i, j) = -std::sqrt(-V(i, j));
			}
		}
	}

	V.colwise() -= V.rowwise().mean();

	constexpr int n_recurse_neighbors = 2;

	std::vector<Eigen::Triplet<double>> graph_triplets, graph_neg_triplets;

	for (int i = 0; i < n_cell; ++i) {

		Eigen::ArrayXd r = V.row(i);

		auto [min, max] = std::ranges::minmax(V.row(i));

		if (min == 0 && max == 0) {
			continue;
		}

		QVector<int> neighbor_index;
		neighbor_index.reserve(this->n_neighbors_ * this->n_neighbors_);
		for (auto neighbor : this->knn_indices_.row(i)) {
			neighbor_index << custom::cast<QVector>(this->knn_indices_.row(neighbor));
		}

		neighbor_index = custom::unique(neighbor_index);

		Eigen::MatrixXd dX = this->ms_(neighbor_index, Eigen::all).array().rowwise() - this->ms_.row(i).array();

		int n_neighbor = dX.rows();
		for (int k = 0; k < n_gene; ++k) {
			for (int j = 0; j < n_neighbor; ++j) {
				if (dX(j, k) > 0) {
					dX(j, k) = std::sqrt(dX(j, k));
				}
				else {
					dX(j, k) = -std::sqrt(-dX(j, k));
				}
			}
		}

		Eigen::ArrayXd vals = ScveloWorker::cosine_correlation(dX, V.row(i));

		for (int j = 0; j < n_neighbor; ++j) {
			double val = vals[j];

			if (val > 0) {
				graph_triplets.emplace_back(i, neighbor_index[j], std::min(val, 1.0));
			}
			else if (val < 0) {
				graph_neg_triplets.emplace_back(i, neighbor_index[j], std::max(val, -1.0));
			}
			
		}
	}

	this->res_.reset(new ScveloEstimate());

	this->res_->graph_.resize(n_cell, n_cell);
	this->res_->graph_neg_.resize(n_cell, n_cell);
	this->res_->graph_.setFromTriplets(graph_triplets.cbegin(), graph_triplets.cend());
	this->res_->graph_neg_.setFromTriplets(graph_neg_triplets.cbegin(), graph_neg_triplets.cend());

	Eigen::ArrayXd confidence = custom::row_max(this->res_->graph_);

	this->res_->self_probability_ = (custom::linear_percentile(confidence, 98) - confidence).cwiseMin(1.0).cwiseMax(0.0);
};

Eigen::ArrayXd ScveloWorker::cosine_correlation(const Eigen::MatrixXd& mat, const Eigen::ArrayXd& arr) {
	Eigen::MatrixXd dx = mat.array().colwise() - mat.array().rowwise().mean();
	double norm = arr.matrix().norm();

	if (norm == 0) {
		return Eigen::ArrayXd::Zero(mat.rows());
	}
	else {
		return (dx.array().rowwise() * arr.transpose()).rowwise().sum() / (dx.rowwise().norm().array() * norm);
	}
};