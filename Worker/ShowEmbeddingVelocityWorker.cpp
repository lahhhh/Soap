#include "ShowEmbeddingVelocityWorker.h"


#include <boost\math\distributions\normal.hpp>

#include "Custom.h"

void ShowEmbeddingVelocityWorker::run() {

	this->show_velocity();
};


void ShowEmbeddingVelocityWorker::show_velocity() {

	const Eigen::MatrixXd& em = this->estimate_->current_;
	const Eigen::MatrixXd& nd = this->estimate_->deltaE_;

	Eigen::ArrayXd col_sum1 = em.colwise().sum(), col_sum2 = nd.colwise().sum();

	const int n_col = em.cols(), n_row = em.rows();

	Eigen::MatrixXd cc = this->col_delta_cor_log10(em, log10(nd.array().abs() + 1.0) * custom::sign(nd).array());

	const int n_cell = em.cols();
	if (this->n_cell_neighbor_ > n_cell) {
		this->n_cell_neighbor_ = n_cell;
	}

	for (int i = 0; i < n_cell; ++i) {
		for (int j = 0; j < n_cell; ++j) {
			if (std::isnan(cc(i, j))) {
				cc(i, j) = 0.0;
			}
		}
	}

	for (int i = 0; i < n_cell; ++i) {
		cc(i, i) = 0.0;
	}

	auto knn = custom::balanced_knn_mt(this->embedding_.transpose(), this->n_cell_neighbor_, n_cell, true, "euclidean");
	Eigen::MatrixXd emb_knn = custom::create_matrix_from_knn_index(knn).toDense().cast<double>();
	for (int i = 0; i < n_cell; ++i) {
		emb_knn(i, i) = 1.0;
	}
	constexpr double corr_sigma = 0.05;
	Eigen::MatrixXd tp = exp(cc.array() / corr_sigma) * emb_knn.array();
	cc.resize(0, 0);
	emb_knn.resize(0, 0);

	for (int i = 0; i < n_cell; ++i) {
		double col_sum = tp.col(i).sum();

		if (col_sum != 0.0) {
			tp.col(i) /= col_sum;
		}
	}

	Eigen::MatrixXd arsd = this->emb_arrows(this->embedding_, tp.sparseView()).transpose();
	tp.resize(0, 0);

	constexpr int n_emb = 2, n_grid = 30;
	constexpr double min_grid_cell_mass = 1.0;

	Eigen::MatrixXd arrows_start(n_grid * n_grid, n_emb), direction(n_grid * n_grid, n_emb);

	Eigen::ArrayX<bool> point_filter(n_grid * n_grid);
	auto range_x = custom_plot::utility::get_range(this->embedding_.col(0), 0.0);
	auto range_y = custom_plot::utility::get_range(this->embedding_.col(1), 0.0);
	Eigen::ArrayXd gx = Eigen::ArrayXd::LinSpaced(n_grid, range_x.lower, range_x.upper);
	Eigen::ArrayXd gy = Eigen::ArrayXd::LinSpaced(n_grid, range_y.lower, range_y.upper);

	double grid_sd = std::sqrt(((gy(1) - gy(0)) * (gy(1) - gy(0)) + (gx(1) - gx(0)) * (gx(1) - gx(0))) / 2);
	double min_arrow_size = 0.01 * std::sqrt((gy(1) - gy(0)) * (gy(1) - gy(0)) + (gx(1) - gx(0)) * (gx(1) - gx(0)));

	boost::math::normal_distribution<double> norm(0, grid_sd);

	for (int i = 0; i < n_grid; ++i) {
		Eigen::MatrixXd cell_distance = Eigen::MatrixXd::Zero(n_cell, n_grid), cell_distance2 = cell_distance;
		cell_distance.array().colwise() = this->embedding_.col(1).array();
		cell_distance.array().rowwise() -= gy.transpose();
		cell_distance = cell_distance.array().square();

		cell_distance2.array().colwise() = this->embedding_.col(0).array();
		cell_distance2.array() -= gx(i);
		cell_distance2 = cell_distance2.array().square();

		cell_distance.array() += cell_distance2.array();
		cell_distance = cell_distance.array().sqrt();

		cell_distance2.resize(0, 0);
		for (int j = 0; j < n_grid; ++j) {
			for (int k = 0; k < n_cell; ++k) {
				cell_distance(k, j) = boost::math::pdf(norm, cell_distance(k, j));
			}
		}

		Eigen::ArrayXd gw = cell_distance.array().colwise().sum();
		Eigen::ArrayXd cws = gw.cwiseMax(1.0);
		Eigen::ArrayXd gxd = (cell_distance.array().colwise() * arsd.col(0).array()).colwise().sum() / cws.transpose();
		Eigen::ArrayXd gyd = (cell_distance.array().colwise() * arsd.col(1).array()).colwise().sum() / cws.transpose();


		Eigen::ArrayXd al = (gxd.square() + gyd.square()).sqrt();

		Eigen::ArrayX<bool> filter = (gw > min_grid_cell_mass) * (al > min_arrow_size);

		point_filter.segment(i * n_grid, n_grid) = filter;
		arrows_start.col(0).segment(i * n_grid, n_grid).array() = gx(i);
		arrows_start.col(1).segment(i * n_grid, n_grid).array() = gy;
		direction.col(0).segment(i * n_grid, n_grid).array() = gxd;
		direction.col(1).segment(i * n_grid, n_grid).array() = gyd;
	}
	custom::quiver_scale(arrows_start, direction);

	if (this->graph_mode_ == 0) {
		arrows_start = custom::row_sliced(arrows_start, point_filter);
		direction = custom::row_sliced(direction, point_filter);

		VELO_GRID_PLOT_ELEMENTS res;
		res.embedding = this->embedding_;
		res.arrows_start = arrows_start;
		res.direction = direction;
		res.embedding_names = this->embedding_names_;
		res.graph_settings = this->graph_settings_;
		emit x_grid_graph_ready(res);
	}
	else {
		Eigen::MatrixX<bool> mask(n_grid, n_grid);
		Eigen::MatrixXd u(n_grid, n_grid), v(n_grid, n_grid);

		for (int i = 0; i < n_grid; ++i) {
			for (int j = 0; j < n_grid; ++j) {
				u(i, j) = direction(j * n_grid + i, 0);
				v(i, j) = direction(j * n_grid + i, 1);
				mask(i, j) = point_filter(j * n_grid + i);
			}
		}

		STREAM_PLOT_ELEMENTS res;
		res.embedding = this->embedding_;
		res.x = gx;
		res.y = gy;
		res.u = u;
		res.v = v;
		res.mask = mask;
		res.embedding_names = this->embedding_names_;
		res.graph_settings = this->graph_settings_;
		emit x_stream_graph_ready(res);
	}
	G_TASK_END;

};

Eigen::MatrixXd ShowEmbeddingVelocityWorker::emb_arrows(
	const Eigen::MatrixXd& emb,
	const Eigen::SparseMatrix<double>& tp,
	double arrow_scale
) {
	const int n_cell = emb.rows(), n_emb = 2;
	Eigen::MatrixXd dm(n_emb, n_cell);

	Eigen::SparseMatrix<double> tpb = tp;
	for (int k = 0; k < tpb.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(tpb, k); it; ++it) {
			it.valueRef() = 1.0;
		}
	}
	tpb = custom::normalize(tpb, 1.0);

	Eigen::ArrayXd zv = Eigen::ArrayXd::Zero(n_emb);
	Eigen::MatrixXd temb = emb.transpose();

	for (int i = 0; i < n_cell; ++i) {
		Eigen::MatrixXd di = temb;
		di.array().colwise() -= di.col(i).array();
		di = custom::norm2(di).array() * arrow_scale;
		di.col(i) = zv;
		dm.col(i) = di * tp.col(i) - di * tpb.col(i);
	}

	return dm;
};

Eigen::MatrixXd ShowEmbeddingVelocityWorker::col_delta_cor_log10(
	const Eigen::MatrixXd& e,
	const Eigen::MatrixXd& d,
	const double pseudo_count
) {
	const int ncol = e.cols();
	Eigen::MatrixXd ret(ncol, ncol);
	Eigen::MatrixXd t;

	for (int i = 0; i < ncol; ++i) {
		t = e;
		t.array().colwise() -= e.col(i).array();
		t = log10(t.array().abs() + pseudo_count) * custom::sign(t).array();
		ret.col(i) = custom::cor_mt(t, d.col(i));
	}

	return ret;
};