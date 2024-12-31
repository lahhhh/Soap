#include "CircosPlot.h"

#include "CustomPlot.h"

// row : source; column : target
void circos_plot(
	QCustomPlot* draw_area,
	QCPAxisRect* axis_rect,
	const QStringList& levels,
	const QList<QColor>& colors,
	const Eigen::MatrixXd& data,
	bool color_by_source
) {
	int n_level = levels.size();

	if (n_level < 2) {
		return;
	}

	if (data.rows() != data.cols() || data.rows() != n_level) {
		return;
	}

	double gap_degree = M_PI / 20;

	double degree_free = 2 * M_PI - gap_degree;

	if ((data.array() < 0.0).any()) {
		return;
	}

	double all_len = data.sum() * 2;

	double per_gap_degree = gap_degree / n_level;

	QVector<double> start_degree(n_level), source_degree(n_level), target_degree(n_level);

	start_degree[0] = 0.0;


	for (int i = 0; i < n_level; ++i) {

		double d1 = data.row(i).sum() / all_len * degree_free;
		double d2 = data.col(i).sum() / all_len * degree_free;

		source_degree[i] = start_degree[i] + d1;
		target_degree[i] = start_degree[i] + d1 + d2;

		if (i < n_level - 1) {
			start_degree[i + 1] = start_degree[i] + d1 + d2 + per_gap_degree;
		}

	}

	for (int i = 0; i < n_level; ++i) {

		QColor color = colors[i];

		double degree = (data.row(i).sum() + data.col(i).sum()) / all_len * degree_free;

		if (degree <= 0.0) {
			continue;
		}

		double start_d = start_degree[i];
		double end_degree = start_d + degree;

		auto [x1, y1] = custom_plot::utility::arc(
			0.0,
			0.0,
			1.0,
			start_d,
			end_degree
		);

		auto [x2, y2] = custom_plot::utility::arc(
			0.0,
			0.0,
			1.05,
			end_degree,
			start_d
		);

		x1 << x2;
		y1 << y2;

		custom_plot::patch::shape_borderless(
			draw_area,
			axis_rect,
			x1,
			y1,
			color
		);

	}

	for (int i = 0; i < n_level; ++i) {

		for (int j = 0; j < n_level; ++j) {

			int k = i + j;
			if (k >= n_level) {
				k -= n_level;
			}

			double d = data(i, k);

			if (d <= 0.0) {
				continue;
			}

			QColor source_color = colors[i];
			QColor target_color = colors[k];

			double degree = d / all_len * degree_free;

			double source_end_degree = source_degree[i];
			double source_start_degree = source_end_degree - degree;

			double target_end_degree = target_degree[k];
			double target_start_degree = target_end_degree - degree;

			auto [x1, y1] = custom_plot::utility::arc(
				0.0,
				0.0,
				0.93,
				source_start_degree,
				source_end_degree
			);

			auto [x2, y2] = custom_plot::utility::arc(
				0.0,
				0.0,
				0.96,
				source_end_degree,
				source_start_degree
			);

			x1 << x2;
			y1 << y2;

			custom_plot::patch::shape_borderless(
				draw_area,
				axis_rect,
				x1,
				y1,
				target_color
			);

			auto [x3, y3] = custom_plot::utility::arc(
				0.0,
				0.0,
				0.90,
				source_start_degree,
				source_end_degree
			);

			auto [x4, y4] = custom_plot::utility::circos_curve(
				0.0,
				0.0,
				0.90 * std::cos(source_end_degree),
				0.90 * std::sin(source_end_degree),
				0.90 * std::cos(target_start_degree),
				0.90 * std::sin(target_start_degree)
			);

			auto [x5, y5] = custom_plot::utility::circos_curve(
				0.0,
				0.0,
				0.90 * std::cos(target_end_degree),
				0.90 * std::sin(target_end_degree),
				0.90 * std::cos(source_start_degree),
				0.90 * std::sin(source_start_degree)
			);

			x3 << x4 << 0.96 * std::cos((target_start_degree + target_end_degree) / 2) << x5;
			y3 << y4 << 0.96 * std::sin((target_start_degree + target_end_degree) / 2) << y5;

			QColor c = color_by_source ? source_color : target_color;

			c.setAlpha(64);

			custom_plot::patch::shape_borderless(
				draw_area,
				axis_rect,
				x3,
				y3,
				c
			);

			source_degree[i] -= degree;
			target_degree[k] -= degree;
		}
	}

	custom_plot::patch::remove_left_bottom_axis(axis_rect);
	custom_plot::patch::set_range(axis_rect, { -1.2, 1.2 }, { -1.2, 1.2 });
};