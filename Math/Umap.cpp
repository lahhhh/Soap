#include "Umap.h"

#include "annoylib.h"
#include "kissrandom.h"
#include "mman.h"
#include "Nls.h"
#include "Custom.h"

#define SMOOTH_K_TOLERANCE 1e-5
#define MIN_K_DIST_SCALE 1e-3

UMAP::UMAP(
	Eigen::MatrixXd* mat,
	int n_neighbors,
	const QString& metric,
	double learning_rate,
	const QString& init,
	double minimum_distance,
	double spread,
	double set_op_mix_ratio,
	double repulsion_strength,
	int negative_sample_rate,
	int random_state,
	int n_trees
) :
	mat_(mat),
	n_neighbors_(n_neighbors),
	metric_(metric),
	learning_rate_(learning_rate),
	init_(init),
	minimum_distance_(minimum_distance),
	spread_(spread),
	set_op_mix_ratio_(set_op_mix_ratio),
	repulsion_strength_(repulsion_strength),
	negative_sample_rate_(negative_sample_rate),
	random_state_(random_state),
	n_trees_(n_trees),
	initial_alpha_(learning_rate)
{}

double func(double input, const Eigen::ArrayXd& params)
{

	double A = params[0];
	double B = params[1];
	return 1 / (1 + A * std::pow(input, 2 * B));
}

void UMAP::find_a_b_param(double spread, double minimum_distance) {

	Eigen::ArrayXd x = Eigen::ArrayXd::LinSpaced(300, 0, spread * 3), y = Eigen::ArrayXd::Zero(300), param(2);

	for (int i = 0; i < 300; ++i) {

		if (x[i] < minimum_distance) {
			y[i] = 1;
		}
		else {
			y[i] = exp(-(x[i] - minimum_distance) / spread);
		}
	}
	
	param[0] = param[1] = 1;
	
	nls_LevenbergMarquardt(func, x, y, param);

	this->a_ = param[0];
	this->b_ = param[1];
};

void UMAP::set_disconnection_distance() {

	if (this->metric_ == "correlation" || this->metric_ == "cosine") {
	
		this->disconnection_distance_ = 2.0;
	}
	else if (this->metric_ == "hellinger" || this->metric_ == "jaccard" || this->metric_ == "dice") {
		
		this->disconnection_distance_ = 1.0;
	}
	else {
		
		this->disconnection_distance_ = std::numeric_limits<double>::max();
	}
};

void UMAP::smooth_knn_distance(int n_iter) { // bandwidth --> 1.0

	int sample_number = this->mat_->rows();

	double target = std::log2((double)this->n_neighbors_);

	Eigen::ArrayXd rho = Eigen::ArrayXd::Zero(sample_number);
	Eigen::ArrayXd result = Eigen::ArrayXd::Zero(sample_number);

	double mean_distances = this->knn_distance_.mean();
	
	int index = (int)floor(this->local_connectivity_);
	
	double interpolation = this->local_connectivity_ - index;

#pragma omp parallel for
	for (int i = 0; i < sample_number; ++i) {

		double lo = 0.0, mid = 1.0, hi = std::numeric_limits<double>::max();
		double cnt = this->knn_distance_.row(i).count();

		Eigen::ArrayXd ith_distances = this->knn_distance_.row(i);
		Eigen::ArrayXd non_zero_dists = _Cs sliced(ith_distances, ith_distances > 0);

		if (cnt > this->local_connectivity_) {

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
			
			for (int n = 1; n < this->n_neighbors_; ++n) {

				double d = this->knn_distance_(i, n) - rho[i];
				
				if (d > 0) {
					psum += exp(-(d / mid));
				}
				else {
					psum += 1.0;
				}
			}

			if (std::abs(psum - target) < SMOOTH_K_TOLERANCE) {
				break;
			}

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

			if (result[i] < MIN_K_DIST_SCALE * mean_ith_distances) {
				result[i] = MIN_K_DIST_SCALE * mean_ith_distances;
			}
		}
		else {
			if (result[i] < MIN_K_DIST_SCALE * mean_distances) {
				result[i] = MIN_K_DIST_SCALE * mean_distances;
			}
		}
	}

	this->smooth_sigmas_ = result;
	this->smooth_rho_ = rho;
};

void UMAP::compute_membership_strengths() {

	const int n_samples = this->knn_distance_.rows();
	
	std::vector<Eigen::Triplet<double> > triplets;

	for (int i = 0; i < n_samples; ++i) {
		for (int j = 0; j < this->n_neighbors_; ++j) {

			double val;
			
			if (this->knn_indices_(i, j) == i) {
				continue;
			}
			else if (this->knn_distance_(i, j) - this->smooth_rho_[i] <= 0.0 || this->smooth_sigmas_[i] == 0.0) {
				val = 1.0;
			}
			else {
				val = exp(-((this->knn_distance_(i, j) - this->smooth_rho_[i]) / (this->smooth_sigmas_[i])));
			}

			triplets.emplace_back(i, this->knn_indices_(i, j), val);
		}
	}

	this->knn_distance_.resize(0, 0);
	this->knn_indices_.resize(0, 0);
	this->strengths_.resize(n_samples, n_samples);
	this->strengths_.setFromTriplets(triplets.cbegin(), triplets.cend());

	Eigen::SparseMatrix<double> transpose = this->strengths_.transpose();
	Eigen::SparseMatrix<double> prod_matrix = this->strengths_.cwiseProduct(transpose);

	this->strengths_ = this->set_op_mix_ratio_ * (this->strengths_ + transpose - prod_matrix) +
		(1.0 - this->set_op_mix_ratio_) * prod_matrix;
		
};

void UMAP::fuzzy_simplicial_set() {
	this->n_vertices_ = this->mat_->rows();

	if (this->metric_ == "Euclidean") {

		std::tie(this->knn_indices_, this->knn_distance_) = _Cs get_knn_mt<Euclidean, true>(*this->mat_, this->n_neighbors_, this->n_trees_);
	}
	else if (this->metric_ == "Angular") {

		std::tie(this->knn_indices_, this->knn_distance_) = _Cs get_knn_mt<Angular, true>(*this->mat_, this->n_neighbors_, this->n_trees_);
	}
	else {
		std::tie(this->knn_indices_, this->knn_distance_) = _Cs get_knn_mt<Manhattan, true>(*this->mat_, this->n_neighbors_, this->n_trees_);
	}

	this->smooth_knn_distance();

	this->compute_membership_strengths();
};

void UMAP::spectral_layout() { 
	// TODO : rspectra 
};

void UMAP::make_epoches_per_sample() {

	int length = this->strengths_.nonZeros();
	
	this->epochs_per_sample_ = Eigen::ArrayXd::Constant(length, -1.0);
	
	double max_value = _Cs max(this->strengths_);
	
	double factor = this->n_epoches_ / max_value;

	const double* data = this->strengths_.valuePtr();

	for (int i = 0; i < length; ++i) {
		
		double n_sample = data[i] * factor;

		if (n_sample > 0) {

			this->epochs_per_sample_[i] = this->n_epoches_ / n_sample;
		}
	}
};

SOAP_INLINE double clip(double val) {

	return val > 4.0 ? 4.0 : (val < -4.0 ? -4.0 : val);
};

void UMAP::optimize_layout_euclidean() {

	const std::size_t size = this->embedding_.size();
	std::vector<double> head_embedding(size), tail_embedding(size);
	std::memcpy(head_embedding.data(), this->embedding_.data(), sizeof(double) * size);
	std::memcpy(tail_embedding.data(), this->embedding_.data(), sizeof(double) * size);

	const std::size_t edge_number = this->strengths_.nonZeros();
	std::vector<int> head(edge_number), tail(edge_number);

	int index = 0;
	for (int k = 0; k < this->strengths_.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(this->strengths_, k); it; ++it) {
			head[index] = it.row();
			tail[index] = k;
			++index;
		}
	}

	int n_vertices = this->embedding_.rows();
	double alpha = this->initial_alpha_;
	Eigen::ArrayXd epochs_per_negative_sample = this->epochs_per_sample_ / this->negative_sample_rate_;
	Eigen::ArrayXd epoch_of_next_negative_sample = epochs_per_negative_sample;
	Eigen::ArrayXd epoch_of_next_sample = this->epochs_per_sample_;

	int iteration_length = this->epochs_per_sample_.size();

	std::default_random_engine e;
	e.seed(this->random_state_);
	std::uniform_int_distribution<unsigned> u(0, n_vertices - 1);

	double dys1 = 0.0, dys2 = 0.0;

	for (int n = 0; n < this->n_epoches_; ++n) {

		for (int i = 0; i < iteration_length; ++i) {

			if (epoch_of_next_sample[i] <= n) {

				int dj = head[i] * /* this->n_components_ */ 2;
				int dk = tail[i] * /* this->n_components_ */ 2;

				double dist_squared = 0.0;

				double head_emb1 = head_embedding[dj], head_emb2 = head_embedding[dj + 1];
				double tail_emb1 = tail_embedding[dk], tail_emb2 = tail_embedding[dk + 1];

				double diff1 = head_emb1 - tail_emb1;
				dys1 = diff1;
				dist_squared += diff1 * diff1;

				double diff2 = head_emb2 - tail_emb2;
				dys2 = diff2;
				dist_squared += diff2 * diff2;

				double grad_coeff;
				if (dist_squared > 0) {
					grad_coeff = -2.0 * this->a_ * this->b_ * std::pow(dist_squared, this->b_ - 1.0);


					grad_coeff /= this->a_ * std::pow(dist_squared, this->b_) + 1.0;
				}
				else {
					grad_coeff = 0.0;
				}

				double grad_d1 = clip(grad_coeff * dys1);
				double grad_d2 = clip(grad_coeff * dys2);

				head_embedding[dj] += grad_d1 * alpha;
				head_embedding[dj + 1] += grad_d2 * alpha;

				tail_embedding[dk] += -grad_d1 * alpha;
				tail_embedding[dk + 1] += -grad_d2 * alpha;

				epoch_of_next_sample[i] += this->epochs_per_sample_[i];
				int n_neg_samples = (int)((n - epoch_of_next_negative_sample[i]) / epochs_per_negative_sample[i]);

				for (int p = 0; p < n_neg_samples; ++p) {

					int dk = u(e) * /* this->n_components_ */ 2;
					dist_squared = 0;

					head_emb1 = head_embedding[dj];
					head_emb2 = head_embedding[dj + 1];

					tail_emb1 = tail_embedding[dk];
					tail_emb2 = tail_embedding[dk + 1];

					diff1 = head_emb1 - tail_emb1;
					dys1 = diff1;
					dist_squared += diff1 * diff1;

					diff2 = head_emb2 - tail_emb2;
					dys2 = diff2;
					dist_squared += diff2 * diff2;

					if (dist_squared > 0.0) {
						grad_coeff = 2.0 * this->repulsion_strength_ * this->b_;
						grad_coeff /= (0.001 + dist_squared) * (this->a_ * std::pow(dist_squared, this->b_) + 1);
					}
					else if (dj == dk) {
						continue;
					}
					else {
						grad_coeff = 0.0;
					}

					if (grad_coeff > 0.0) {
						head_embedding[dj] += clip(grad_coeff * dys1) * alpha;
						head_embedding[dj + 1] += clip(grad_coeff * dys2) * alpha;
					}
					else {
						head_embedding[dj] += 4.0 * alpha;
						head_embedding[dj + 1] += 4.0 * alpha;
					}
				}

				epoch_of_next_negative_sample[i] += (n_neg_samples * epochs_per_negative_sample[i]);
			}
		}
		alpha = this->initial_alpha_ * (1.0 - ((double)n / this->n_epoches_));
	}

	this->embedding_.resize(this->n_vertices_, this->n_components_);

	for (int i = 0; i < this->n_components_; ++i) {
		for (int j = 0; j < this->n_vertices_; ++j) {
			this->embedding_(j, i) = head_embedding[this->n_components_ * j + i];
		}
	}
};

void UMAP::fit_embed_data() {
	this->n_epoches_ = this->strengths_.rows() > 10000 ? 200 : 500;
	this->strengths_ = _Cs eliminate_less_than(this->strengths_, _Cs max(this->strengths_) / this->n_epoches_);

	srand(this->random_state_);
	this->embedding_ = Eigen::MatrixXd::Random(this->strengths_.rows(), 2) * 10; // --> this->init_ == random ... to do : spectral_layout

	this->make_epoches_per_sample();

	Eigen::ArrayXd embedding_min = this->embedding_.colwise().minCoeff();
	Eigen::ArrayXd embedding_max = this->embedding_.colwise().maxCoeff();

	for (int i = 0; i < embedding_min.size(); ++i) {
		this->embedding_.col(i) = (10 * (this->embedding_.col(i).array() - embedding_min[i]) / (embedding_max[i] - embedding_min[i])).eval();
	}

	optimize_layout_euclidean();
};

void UMAP::fit() {

	find_a_b_param(this->spread_, this->minimum_distance_);

	fuzzy_simplicial_set();

	fit_embed_data();
};
