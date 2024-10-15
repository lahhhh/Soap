#include "StreamPlot.h"

#include "CustomPlot.h"

void stream_plot(
	const Eigen::ArrayXd& x,
	const Eigen::ArrayXd& y,
	const Eigen::MatrixXd& u,
	const Eigen::MatrixXd& v,
	const Eigen::MatrixX<bool>& mask,
	QCustomPlot* draw_area,
	QCPAxisRect* axis_rect)
{
	constexpr int density = 1;
	constexpr double min_length = 0.1, max_length = 4.0;
	DomainMap dmap(Grid(x, y), StreamMask(density));

	Eigen::MatrixXd ug = dmap.x_data2grid_ * u.array();
	Eigen::MatrixXd vg = dmap.y_data2grid_ * v.array();

	Eigen::MatrixXd u_ax = ug.array() / (dmap.grid_.nx_ - 1);
	Eigen::MatrixXd v_ax = vg.array() / (dmap.grid_.ny_ - 1);

	Eigen::MatrixXd speed = (v_ax.array().square() + u_ax.array().square()).sqrt();

	std::vector<std::pair<std::vector<double>, std::vector<double>>> trajectories;

	auto forward_time = [&ug, &vg, &speed, &dmap, &mask](double xi, double yi)->std::tuple<int, double, double> {
		if (!dmap.grid_.within_grid(xi, yi)) {
			return std::make_tuple(OutOfBounds, 0.0, 0.0);
		}

		auto [err, ds_dt] = interpgrid(speed, mask, xi, yi);
		if (err == TerminateTrajectory || ds_dt == 0.0) {
			return std::make_tuple(TerminateTrajectory, 0.0, 0.0);
		}

		double dt_ds{ 1.0 / ds_dt };

		auto [err2, ui] = interpgrid(ug, mask, xi, yi);
		if (err2 == TerminateTrajectory) {
			return std::make_tuple(TerminateTrajectory, 0.0, 0.0);
		}

		auto [err3, vi] = interpgrid(vg, mask, xi, yi);
		if (err3 == TerminateTrajectory) {
			return std::make_tuple(TerminateTrajectory, 0.0, 0.0);
		}

		return std::make_tuple(Normal, ui * dt_ds, vi * dt_ds);
	};

	auto backward_time = [&forward_time](double xi, double yi) {
		auto [err, r1, r2] = forward_time(xi, yi);
		return std::make_tuple(err, -r1, -r2);
	};

	auto integrate = [&dmap, min_length, max_length, &forward_time, &backward_time](double x0, double y0)-> std::pair<std::vector<double>, std::vector<double>> {


		if (!dmap.start_trajectory(x0, y0)) {
			return std::pair<std::vector<double>, std::vector<double>>();
		}
		double stotal{ 0.0 };
		std::vector<double> x_traj, y_traj;

		auto [s, xt, yt] = integrate_rk12(x0, y0, dmap, max_length, backward_time);


		stotal += s;
		x_traj.insert(x_traj.end(), xt.rbegin(), xt.rend());
		y_traj.insert(y_traj.end(), yt.rbegin(), yt.rend());

		dmap.reset_start_point(x0, y0);

		std::tie(s, xt, yt) = integrate_rk12(x0, y0, dmap, max_length, forward_time);

		stotal += s;

		if (stotal > min_length) {
			if (xt.size() > 1) {

				x_traj.insert(x_traj.end(), xt.begin() + 1, xt.end());
				y_traj.insert(y_traj.end(), yt.begin() + 1, yt.end());
			}

			return std::make_pair(x_traj, y_traj);
		}
		else {
			dmap.undo_trajectory();
			return std::pair<std::vector<double>, std::vector<double>>();
		}
	};

	auto [nx, ny] = dmap.mask_.shape();
	int xfirst{ 0 }, yfirst{ 1 }, xlast{ nx - 1 }, ylast{ ny - 1 };
	int xm{ 0 }, ym{ 0 };

	// 0 - right ; 1 - up ; 2 - left ; 3 - down
	int direction{ 0 };


	for (int i = 0; i < nx * ny; ++i) {

		if (dmap.mask_(xm, ym) == 0) {
			auto [xg, yg] = dmap.mask2grid(xm, ym);

			auto t = integrate(xg, yg);

			if (!t.first.empty()) {
				trajectories.push_back(t);
			}
		}

		if (direction == 0) {
			++xm;
			if (xm >= xlast) {
				--xlast;
				direction = 1;
			}
		}
		else if (direction == 1) {
			++ym;
			if (ym >= ylast) {
				--ylast;
				direction = 2;
			}
		}
		else if (direction == 2) {
			--xm;
			if (xm <= xfirst) {
				++xfirst;
				direction = 3;
			}
		}
		else if (direction == 3) {
			--ym;
			if (ym <= yfirst) {
				++yfirst;
				direction = 0;
			}
		}
	}

	double max_x = axis_rect->axis(QCPAxis::atBottom)->range().upper, min_x = axis_rect->axis(QCPAxis::atBottom)->range().lower;
	double max_y = axis_rect->axis(QCPAxis::atLeft)->range().upper, min_y = axis_rect->axis(QCPAxis::atLeft)->range().lower;

	double graph_y_x_ratio = (max_y - min_y) / (max_x - min_x);
	double sin = 1 / std::sqrt(1 + graph_y_x_ratio * graph_y_x_ratio), cos = graph_y_x_ratio / std::sqrt(1 + graph_y_x_ratio * graph_y_x_ratio);

	QPen arrow_pen;

	arrow_pen.setColor(Qt::black);
	arrow_pen.setWidth(1);

	for (const auto& t : trajectories) {
		std::vector<double> tx, ty;

		auto&& [xx, yy] = t;
		if (xx.empty() || yy.empty()) {
			continue;
		}

		tx.push_back(xx[0]);
		ty.push_back(yy[0]);

		int np = xx.size();
		double lastx = xx[0];
		double lasty = yy[0];

		for (int i = 1; i < np; ++i) {
			if (xx[i] != lastx || yy[i] != lasty) {
				lastx = xx[i];
				lasty = yy[i];
				tx.push_back(lastx);
				ty.push_back(lasty);
			}
		}

		std::ranges::for_each(tx, [&dmap](auto&& d) {d = d / dmap.x_data2grid_ + dmap.grid_.x_origin_; });
		std::ranges::for_each(ty, [&dmap](auto&& d) {d = d / dmap.y_data2grid_ + dmap.grid_.y_origin_; });

		custom_plot::patch::curve(draw_area, axis_rect, custom::cast<QVector>(tx), custom::cast<QVector>(ty), Qt::black, 2, Qt::SolidLine);

		const int n_point = tx.size();
		std::vector<double> s(n_point, 0.0);    
		for (int i = 1; i < n_point; ++i) {
			double x_start = tx[i - 1], x_end = tx[i], y_start = ty[i - 1], y_end = ty[i];
			double distance = std::sqrt((x_end - x_start) * (x_end - x_start) + (y_end - y_start) * (y_end - y_start));
			s[i] = distance + s[i - 1];
		}

		double arrow_dis = s[n_point - 1] / 2;
		int arrow_loc{ 0 };
		for (int i = 1; i < n_point; ++i) {
			if (s[i] > arrow_dis) {
				if (arrow_dis - s[i - 1] > s[i] - arrow_dis) {
					arrow_loc = i;
				}
				else {
					arrow_loc = i - 1;
				}

				break;
			}
		}

		if (arrow_loc <= n_point - 2) {
			double arrow_start_x = tx[arrow_loc], arrow_start_y = ty[arrow_loc];
			double arrow_end_x = (arrow_start_x + tx[arrow_loc + 1]) / 2, arrow_end_y = (arrow_start_y + ty[arrow_loc + 1]) / 2;

			double x_dis = std::abs(arrow_end_x - arrow_start_x);
			double y_dis = std::abs(arrow_end_y - arrow_start_y);

			if (x_dis < 0.2 && y_dis < 0.2) {
				double times = 0.2 / std::max(x_dis, y_dis);
				double arr_mid_x = (arrow_start_x + arrow_end_x) / 2;
				double arr_mid_y = (arrow_start_y + arrow_end_y) / 2;
				arrow_start_x = arr_mid_x + (arrow_start_x - arr_mid_x) * times;
				arrow_end_x = arr_mid_x + (arrow_end_x - arr_mid_x) * times;
				arrow_start_y = arr_mid_y + (arrow_start_y - arr_mid_y) * times;
				arrow_end_y = arr_mid_y + (arrow_end_y - arr_mid_y) * times;
			}

			double arrow_vertex_x_1 = 1.5 * arrow_start_x - 0.5 * arrow_end_x + (arrow_end_y - arrow_start_y) * sin;
			double arrow_vertex_y_1 = 1.5 * arrow_start_y - 0.5 * arrow_end_y - (arrow_end_x - arrow_start_x) * cos;
			double arrow_vertex_x_2 = 1.5 * arrow_start_x - 0.5 * arrow_end_x - (arrow_end_y - arrow_start_y) * sin;
			double arrow_vertex_y_2 = 1.5 * arrow_start_y - 0.5 * arrow_end_y + (arrow_end_x - arrow_start_x) * cos;

			draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			QVector<QCPCurveData> data(5);

			data[0] = QCPCurveData(0, arrow_start_x, arrow_start_y);
			data[1] = QCPCurveData(1, arrow_vertex_x_1, arrow_vertex_y_1);
			data[2] = QCPCurveData(2, arrow_end_x, arrow_end_y);
			data[3] = QCPCurveData(3, arrow_vertex_x_2, arrow_vertex_y_2);
			data[4] = QCPCurveData(4, arrow_start_x, arrow_start_y);

			shape->setPen(arrow_pen);
			shape->setBrush(QBrush(Qt::black));
			shape->data()->set(data);
		}
	}
};

std::pair<int, double> interpgrid(const Eigen::MatrixXd& a, const Eigen::MatrixX<bool>& mask, double xi, double yi) {

	int nx = a.cols(), ny = a.rows();

	int x{ int(xi) }, y{ int(yi) };
	int xn{ x == (nx - 1) ? x : x + 1 }, yn{ y == (ny - 1) ? y : y + 1 };

	if (!(mask(y, x) || mask(y, xn) || mask(yn, x) || mask(yn, xn))) {
		return std::make_pair(TerminateTrajectory, 0.0);
	}

	double a00{ a(y, x) }, a01{ a(y, xn) }, a10{ a(yn, x) }, a11{ a(yn, xn) };

	double xt{ xi - x }, yt{ yi - y };

	double a0{ a00 * (1 - xt) + a01 * xt };
	double a1{ a10 * (1 - xt) + a11 * xt };
	double ai{ a0 * (1 - yt) + a1 * yt };

	return std::make_pair(Normal, ai);
}

