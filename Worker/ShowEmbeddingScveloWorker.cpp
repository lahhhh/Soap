#include "ShowEmbeddingScveloWorker.h"

#include <boost\math\distributions\normal.hpp>

void ShowEmbeddingScveloWorker::velocity_embedding() {

	this->velocity_embedding_ = this->embedding_;

	Eigen::SparseMatrix<double> transition_matrix = this->transition_matrix();

	custom::set_diagonal(transition_matrix, 0.0);

	transition_matrix = custom::eliminate_zero(transition_matrix);

	const int n_cell = this->estimate_->graph_.rows();

	for (int i = 0; i < n_cell; ++i) {

		Eigen::MatrixXd dX = this->embedding_(custom::get_row_index(transition_matrix, i), Eigen::all).array().rowwise() - this->embedding_.row(i).array();

		int n_row = dX.rows();

		for (int j = 0; j < n_row; ++j) {
			double norm = dX.row(j).norm();
			if (norm > 0) {
				dX.row(j).array() /= norm;
			}
			else {
				dX.row(j).setZero();
			}
		}

		Eigen::ArrayXd probs = custom::get_row_data(transition_matrix, i);


		this->velocity_embedding_.row(i) = (probs.matrix().transpose() * dX).array() - probs.mean() * dX.array().colwise().sum();

	}

	custom::quiver_scale(this->embedding_, this->velocity_embedding_);
};

void ShowEmbeddingScveloWorker::compute_velocity_on_grid() {

	auto [x_min, x_max] = std::ranges::minmax(this->embedding_.col(0));
	auto [y_min, y_max] = std::ranges::minmax(this->embedding_.col(1));
	double x_span = x_max - x_min, y_span = y_max - y_min;

	x_min -= 0.01 * x_span;
	x_max += 0.01 * x_span;
	y_min -= 0.01 * y_span;
	y_max += 0.01 * y_span;

	constexpr int n_grid = 40;
	Eigen::ArrayXd x_val = Eigen::ArrayXd::LinSpaced(n_grid, x_min, x_max);
	Eigen::ArrayXd y_val = Eigen::ArrayXd::LinSpaced(n_grid, y_min, y_max);
	Eigen::MatrixXd x_grid(n_grid * n_grid, 2), v_grid(n_grid * n_grid, 2);
	for (int i = 0; i < n_grid; ++i) {
		x_grid.col(0).segment(i * n_grid, n_grid) = x_val;
		x_grid.col(1).segment(i * n_grid, n_grid).setConstant(y_val[i]);
	}

	double scale = (x_val[1] + y_val[1] - x_val[0] - y_val[0]) / 2;
	boost::math::normal_distribution<double> norm(0, scale);

	constexpr int n_neighbor = 30;
	auto [neighbors, weight] = custom::get_knn_mt<Euclidean, true>(this->embedding_, x_grid, n_neighbor);
	for (int j = 0; j < n_neighbor; ++j) {
		for (int i = 0; i < n_grid * n_grid; ++i) {
			weight(i, j) = boost::math::pdf(norm, weight(i, j));
		}
	}
	
	Eigen::ArrayXd p_mass = weight.rowwise().sum(), p_mass_adj = p_mass.cwiseMax(1.0);

	for (int i = 0; i < n_grid * n_grid; ++i) {
		v_grid.row(i) = (this->velocity_embedding_(neighbors.row(i), Eigen::all).array().colwise() * weight.row(i).transpose().array()).colwise().sum() / p_mass_adj[i];
	}

	if (this->graph_mode_ == 0) {

		const double min_mass = custom::linear_percentile(p_mass, 99) / 100;

		x_grid = custom::row_sliced(x_grid, p_mass > min_mass);
		v_grid = custom::row_sliced(v_grid, p_mass > min_mass);

		custom::quiver_scale(x_grid, v_grid);

		VELO_GRID_PLOT_ELEMENTS res;
		res.embedding = this->embedding_;
		res.arrows_start = x_grid;
		res.direction = v_grid;
		res.embedding_names = this->embedding_names_;
		res.graph_settings = this->graph_settings_;
		emit x_grid_graph_ready(res);
	}
	else /* if(this->graph_mode == 2) */ {
		

		Eigen::ArrayXd mass = v_grid.col(0).array().square() + v_grid.col(1).array().square();
		double min_mass = 1e-5, mass90 = mass.maxCoeff() * 0.9;

		if (min_mass > mass90) {
			min_mass = mass90;
		}

		Eigen::ArrayX<bool> filter = mass > min_mass;

		Eigen::ArrayXd length(n_grid * n_grid);

		for (int i = 0; i < n_grid * n_grid; ++i) {
			length[i] = this->velocity_embedding_(neighbors.row(i), Eigen::all).array().abs().rowwise().mean().sum();
		}

		filter *= (length > custom::linear_percentile(length, 5));

		Eigen::MatrixX<bool> mask(n_grid, n_grid);
		Eigen::MatrixXd u(n_grid, n_grid), v(n_grid, n_grid);

		for (int i = 0; i < n_grid; ++i) {
			for (int j = 0; j < n_grid; ++j) {
				u(i, j) = v_grid(i * n_grid + j, 0);
				v(i, j) = v_grid(i * n_grid + j, 1);
				mask(i, j) = filter(i * n_grid + j);
			}
		}

		STREAM_PLOT_ELEMENTS res;
		res.embedding = this->embedding_;
		res.x = x_val;
		res.y = y_val;
		res.u = u;
		res.v = v;
		res.mask = mask;
		res.embedding_names = this->embedding_names_;
		res.graph_settings = this->graph_settings_;
		emit x_stream_graph_ready(res);
	}

	
	G_TASK_END;
};

Eigen::SparseMatrix<double> ShowEmbeddingScveloWorker::transition_matrix() {
	Eigen::SparseMatrix<double> T = custom::set_diagonaled(this->estimate_->graph_, this->estimate_->self_probability_);


	constexpr double scale = 10.0;
	custom::elementwise_in_place(T, [](auto& it) {return std::exp(it.value() * scale) - 1; });


	T -= custom::elementwise(this->estimate_->graph_neg_, [](const auto& it) {
		return std::exp(-it.value() * scale) - 1;
	});


	return custom::row_normalize2(T, 1.0);

};

void ShowEmbeddingScveloWorker::run() {

	this->velocity_embedding();

	this->compute_velocity_on_grid();

}