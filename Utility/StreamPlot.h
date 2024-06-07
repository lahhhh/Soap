#pragma once

#include "Identifier.h"

#include "Custom.h"
#include "qcustomplot.h"

#define Normal 0
#define OutOfBounds 1
#define TerminateTrajectory 2
#define InvalidIndexError 3

class Grid {
public:

	Grid() = default;
	Grid(const Grid&) = default;
	Grid(Grid&&) = default;
	Grid& operator=(const Grid&) = default;
	Grid& operator=(Grid&&) = default;
	~Grid() = default;

	Grid(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y) :
		nx_(x.size()),
		ny_(y.size()),
		x_origin_(x[0]),
		y_origin_(y[0]),
		dx_(x[1] - x[0]),
		dy_(y[1] - y[0]),
		width_(x[x.size() - 1] - x[0]),
		height_(y[y.size() - 1] - y[0])
	{}

	int nx_;
	int ny_;
	double x_origin_;
	double y_origin_;
	double dx_;
	double dy_;
	double width_;
	double height_;

	std::pair<int, int> shape() const{
		return std::make_pair(this->ny_, this->nx_);
	}

	bool within_grid(double xi, double yi) const {
		return (xi >= 0) && (xi <= (this->nx_ - 1)) && (yi >= 0) && (yi <= (this->ny_ - 1));
	}

};

class StreamMask {
public:

	StreamMask() = default;
	StreamMask(const StreamMask&) = default;
	StreamMask(StreamMask&&) = default;
	StreamMask& operator=(const StreamMask&) = default;
	StreamMask& operator=(StreamMask&&) = default;
	~StreamMask() = default;

	explicit StreamMask(int density) :
		nx_(30 * density),
		ny_(30 * density),
		mask_(Eigen::MatrixXi::Zero(30 * density, 30 * density))
	{}

	std::vector<std::pair<int, int>> trajectory_;

	int nx_;
	int ny_;
	int current_x_{ -1 };
	int current_y_{ -1 };

	Eigen::MatrixXi mask_;

	std::pair<int, int> shape() const{
		return std::make_pair(this->ny_, this->nx_);
	}

	int operator()(int row, int col) const{
		return this->mask_(row, col);
	}

	bool update_trajectory(int xm, int ym) {
		if (this->current_x_ != xm || this->current_y_ != ym) {

			if (this->mask_(ym, xm) == 0) {

				this->trajectory_.emplace_back(ym, xm);
				this->mask_(ym, xm) = 1;
				this->current_x_ = xm;
				this->current_y_ = ym;
			}
			else {
				return false;
			}
		}

		return true;
	}

	void undo_trajectory() {
		for (const auto& t : this->trajectory_) {
			this->mask_(t.first, t.second) = 0;
		}
	}

	bool start_trajectory(int xm, int ym) {
		this->trajectory_.clear();

		return this->update_trajectory(xm, ym);
	}
};

class DomainMap {
public:
	DomainMap(
		const Grid& grid,
		const StreamMask& mask
	) :
		grid_(grid),
		mask_(mask)
	{
		this->x_grid2mask_ = (mask.nx_ - 1.0) / (grid.nx_ - 1.0);
		this->y_grid2mask_ = (mask.ny_ - 1.0) / (grid.ny_ - 1.0);

		this->x_mask2grid_ = 1.0 / this->x_grid2mask_;
		this->y_mask2grid_ = 1.0 / this->y_grid2mask_;

		this->x_data2grid_ = 1.0 / grid.dx_;
		this->y_data2grid_ = 1.0 / grid.dy_;
	}

	Grid grid_;
	StreamMask mask_;

	double x_grid2mask_{ 0.0 };
	double y_grid2mask_{ 0.0 };
	double x_mask2grid_{ 0.0 };
	double y_mask2grid_{ 0.0 };
	double x_data2grid_{ 0.0 };
	double y_data2grid_{ 0.0 };

	std::pair<int, int> grid2mask(int xi, int yi) const {
		return std::make_pair(int(xi * this->x_grid2mask_ + 0.5), int(yi * this->y_grid2mask_ + 0.5));
	}

	std::pair<double, double> mask2grid(double xm, double ym) const {
		return std::make_pair(xm * this->x_mask2grid_, ym * this->y_mask2grid_);
	}

	std::pair<double, double> data2grid(double xd, double yd) const {
		return std::make_pair(xd * this->x_data2grid_, yd * this->y_data2grid_);
	}

	std::pair<double, double> grid2data(double xg, double yg) const {
		return std::make_pair(xg / this->x_data2grid_, yg / this->y_data2grid_);
	}

	bool start_trajectory(double xg, double yg) {
		auto [xm, ym] = this->grid2mask(xg, yg);
		return this->mask_.start_trajectory(xm, ym);
	}

	void reset_start_point(double xg, double yg) {
		std::tie(this->mask_.current_x_, this->mask_.current_y_) = this->grid2mask(xg, yg);
	}

	bool update_trajectory(double xg, double yg) {

		if (!this->grid_.within_grid(xg, yg)) {
			return false;
		}

		auto [xm, ym] = this->grid2mask(xg, yg);

		return this->mask_.update_trajectory(xm, ym);
	}

	void undo_trajectory() {
		this->mask_.undo_trajectory();
	}

};


std::pair<int, double> interpgrid(const Eigen::MatrixXd& a, const Eigen::MatrixX<bool>& mask, double xi, double yi);

template <typename Fun>
std::tuple<double, std::vector<double>, std::vector<double>> integrate_rk12(
	double x0,
	double y0, 
	DomainMap& dmap, 
	double max_length,
	const Fun& f) 
{
	constexpr double max_error = 0.003;

	double max_ds = _Cs min(1.0 / dmap.mask_.nx_, 1.0 / dmap.mask_.ny_, 0.1);

	double ds{ max_ds };

	double stotal{ 0.0 }, xi{ x0 }, yi{ y0 };
	std::vector<double> xf_trajectory, yf_trajectory;

	auto [nx, ny] = dmap.grid_.shape();

	while (true) {
		if (dmap.grid_.within_grid(xi, yi)) {
			xf_trajectory.push_back(xi);
			yf_trajectory.push_back(yi);
		}
		else {

outofbounds:if (!xf_trajectory.empty()) {
				auto [ny, nx] = dmap.grid_.shape();

				double _xi = xf_trajectory.back();
				double _yi = yf_trajectory.back();
				auto [_, cx, cy] = f(_xi, _yi);

				double dsx{ 0.0 }, dsy{ 0.0 };

				if (cx == 0.0) {
					dsx = std::numeric_limits<double>::infinity();
				}
				else if (cx < 0.0) {
					dsx = _xi / (-cx);
				}
				else {
					dsx = (nx - 1 - _xi) / cx;
				}

				if (cy == 0.0) {
					dsy = std::numeric_limits<double>::infinity();
				}
				else if (cy < 0.0) {
					dsy = _yi / (-cy);
				}
				else {
					dsy = (ny - 1 - _yi) / cy;
				}

				double _ds = std::min(dsx, dsy);

				xf_trajectory.push_back(_xi + cx * _ds);
				yf_trajectory.push_back(_yi + cy * _ds);

				stotal += _ds;
			}

			break;
		}

		auto [err, k1x, k1y] = f(xi, yi);
		if (err == TerminateTrajectory) {
			break;
		}
		else if (err == OutOfBounds) {
			goto outofbounds;
		}

		auto [err2, k2x, k2y] = f(xi + ds * k1x, yi + ds * k1y);
		if (err2 == TerminateTrajectory) {
			break;
		}
		else if (err2 == OutOfBounds) {
			goto outofbounds;
		}

		double dx1{ ds * k1x };
		double dy1{ ds * k1y };
		double dx2{ ds * 0.5 * (k1x + k2x) };
		double dy2{ ds * 0.5 * (k1y + k2y) };		

		double error = std::sqrt(((dx2 - dx1) / (nx - 1)) * ((dx2 - dx1) / (nx - 1)) + ((dy2 - dy1) / (ny - 1)) * ((dy2 - dy1) / (ny - 1)));

		if (error < max_error) {
			xi += dx2;
			yi += dy2;

			if (!dmap.update_trajectory(xi, yi)) {
				break;
			}

			if (stotal + ds > max_length) {
				break;
			}

			stotal += ds;
		}

		if (error == 0.0) {
			ds = max_ds;
		}
		else {
			ds = std::min(max_ds, 0.85 * ds * std::sqrt(max_error / error));
		}
	}

	return std::make_tuple(stotal, xf_trajectory, yf_trajectory);
};

void stream_plot(
	const Eigen::ArrayXd& x,
	const Eigen::ArrayXd& y,
	const Eigen::MatrixXd& u,
	const Eigen::MatrixXd& v,
	const Eigen::MatrixX<bool>& mask,
	QCustomPlot* draw_area,
	QCPAxisRect* axis_rect
);

