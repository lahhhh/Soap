#pragma once

#include "Identifier.h"
#include <qCustomPlot.h>

#define _CpPatch ::custom_plot::patch::

namespace custom_plot {

	namespace patch {

		QCPGraph* shape(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor border_color,
			QColor fill_color,
			int width
		);

		QCPGraph* shape(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			QColor border_color,
			QColor fill_color,
			int width
		);

		QCPGraph* shape_borderless(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor fill_color
		);

		QCPGraph* shape_borderless(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			QColor fill_color
		);

		QCPGraph* rectangle(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			double start_left_x,
			double start_bottom_y,
			double width,
			double height,
			QColor color,
			int border_width,
			QColor border_color
		);

		QCPGraph* rectangle_borderless(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			double start_left_x,
			double start_bottom_y,
			double width,
			double height,
			QColor color
		);

		QCPGraph* scatter(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			const QColor& color,
			int scatter_size
		);

		QCPGraph* scatter(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			const QColor& color,
			int scatter_size
		);

		void scatter_category(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			const QStringList& values,
			const QStringList& levels,
			const QList<QColor>& colors,
			int scatter_size
		);

		void scatter_category(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			const QStringList& values,
			const QStringList& levels,
			const QList<QColor>& colors,
			int scatter_size
		);

		void scatter_gradient(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			const Eigen::ArrayXd& value,			
			double min,
			double max,
			QColor low_color,
			QColor middle_color,
			QColor high_color,
			int scatter_size
		);

		// note : x can only be ordered, otherwise use curve();
		QCPGraph* line(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor color,
			int width = 1,
			Qt::PenStyle style = Qt::SolidLine
		);

		// note : x can only be ordered, otherwise use curve();
		QCPGraph* line(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			QColor color,
			int width = 1,
			Qt::PenStyle style = Qt::SolidLine
		);

		QCPGraph* curve(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor color,
			int width = 1,
			Qt::PenStyle style = Qt::SolidLine
		);

		QCPGraph* bar(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const Eigen::ArrayXd& bar_location,
			const Eigen::ArrayXd& bar_length,
			int bar_width,
			bool vertical = true
		);

		QCPGraph* bar(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const QVector<double>& bar_location,
			const QVector<double>& bar_length,
			int bar_width,
			bool vertical = true
		);

		void bar_polychrome(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& bar_location,
			const Eigen::ArrayXd& bar_length,
			const QList<QColor>& colors,
			int bar_width,
			bool vertical = true);

		void bar_gradient(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& bar_location,
			const Eigen::ArrayXd& bar_length,
			const Eigen::ArrayXd& value,
			QColor low_color,
			QColor middle_color,
			QColor high_color,
			int bar_width,
			bool vertical = true);

		void bar_single_stack(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QMap<QString, int>& distribution,
			const QStringList& levels,
			const QVector<QColor>& colors,
			double start_point_x,
			double width
		);

		void bar_stack(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QStringList& data,
			const QStringList& data_levels,
			const QVector<QColor>& colors,
			double start_x,
			double width,
			double margin
		);

		QCPGraph* sector(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			double x,
			double y,
			double start_angle,
			double stop_angle,
			double radius,
			int segments
		);

		void pie(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QMap<QString, int>& distribution,
			const QStringList& levels,
			const QVector<QColor>& colors,
			double start_point_x,
			double start_point_y,
			double radius,
			bool white_border = true
		);

		std::pair<double, double> violin(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const QVector<double>& value,
			double center,
			int unit = 16,
			double zero_space = 0.01
		);

		void box(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const QVector<double>& value,
			double center,
			bool draw_outlier,
			int outlier_scatter_size
		);

		void box_batch(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QVector<QColor>& colors,
			const Eigen::ArrayXd& data,
			double start,
			double margin,
			bool draw_outlier,
			int outlier_scatter_size
		);

		std::pair<double, double> violin_facet(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QVector<QColor>& colors,
			const Eigen::ArrayXd& data,
			double center,
			int unit = 16
		);

		std::pair<double, double> violin_batch(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QVector<QColor>& colors,
			const Eigen::ArrayXd& data,
			double start,
			double margin,
			int unit = 16
		);

		void add_title(QCustomPlot* draw_area, QCPLayoutGrid* layout, const QString& title, QFont font);

		void add_label(
			QCustomPlot* draw_area,
			QCPAxisRect* rect,
			const QString& label,
			double x,
			double y,
			QFont font,
			Qt::Alignment alignment = Qt::AlignVCenter | Qt::AlignHCenter,
			double degree = 0);

		void dot(QCustomPlot* draw_area, QCPAxisRect* axis_rect, double x, double y, int size, QColor color);

		void set_proportion_left_axis(QCPAxisRect* axis_rect, QFont tick_label_font);

		void remove_axis(QCPAxisRect* axis_rect, QCPAxis::AxisType type);

		void clear_axis(QCPAxisRect* axis_rect, QCPAxis::AxisType type);

		void remove_left_axis(QCPAxisRect* axis_rect);

		void remove_bottom_axis(QCPAxisRect* axis_rect);

		void clear_left_axis(QCPAxisRect* axis_rect);

		void clear_bottom_axis(QCPAxisRect* axis_rect);

		void remove_all_axis(QCPAxisRect* axis_rect);

		void remove_left_bottom_axis(QCPAxisRect* axis_rect);

		void remove_left_ticks(QCPAxisRect* axis_rect);

		void remove_bottom_ticks(QCPAxisRect* axis_rect);

		void remove_left_subticks(QCPAxisRect* axis_rect);

		void remove_bottom_subticks(QCPAxisRect* axis_rect);

		void remove_left_grid(QCPAxisRect* axis_rect);

		void remove_bottom_grid(QCPAxisRect* axis_rect);

		void remove_grid(QCPAxisRect* axis_rect);

		void set_range(QCPAxisRect* axis_rect, QCPRange bottom_range, QCPRange left_range);

		void set_range(QCPAxisRect* axis_rect, std::pair<QCPRange, QCPRange> range);

		void set_fixed_size(QCPAxisRect* axis_rect, int width, int height);

		void set_single_square_legend(
			QCustomPlot* draw_area,
			QCPAxisRect* legend_rect,
			const QString& label,
			QColor color,
			QFont font,
			double start_point_x,
			double start_point_y,
			int legend_size = 10
		);

		void set_single_round_legend(
			QCustomPlot* draw_area,
			QCPAxisRect* legend_rect,
			const QString& label,
			QColor color,
			QFont font,
			double start_point_x,
			double start_point_y,
			int legend_size = 10
		);

		void set_border_only(
			QCPAxisRect* axis_rect,
			QColor color,
			int width
		);

		QCPAxisRect* new_axis_rect(QCustomPlot* draw_area);

		QCPLayoutGrid* set_legend_layout(
			QCustomPlot* draw_area,
			QCPLayoutGrid* legend_layout);

		void add_gradient_legend(
			QCustomPlot* draw_area,
			QCPLayoutGrid* legend_layout,
			double minimum_value,
			double maximum_value,
			const QString& title,
			const QString& lower_label,
			const QString& upper_label,
			const QFont& legend_title_font,
			const QFont& legend_label_font,
			QColor low_color,
			QColor middle_color,
			QColor high_color
		);

		void set_fixed_width(QCPAxisRect* axis_rect, int width);

		void set_fixed_height(QCPAxisRect* axis_rect, int height);

		std::pair<int, int> find_next_empty_position(QCPLayoutGrid* layout);

		void set_axis_label(
			QCPAxis::AxisType type,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& location,
			const QStringList& labels,
			int tick_length,
			QFont label_font,
			int label_angle
		);

		void add_round_legend(
			QCustomPlot* draw_area,
			QCPLayoutGrid* legend_layout,
			const QStringList& levels,
			const QList<QColor>& colors,
			const QString& legend_title,
			int legend_column_width,
			int legend_row_width,
			QFont legend_title_font,
			QFont legend_label_font
		);

		void single_violin_plot(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& values,
			int unit,
			QColor color,
			const QString& bottom_title,
			const QString& left_title);

		void single_box_plot(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& values,
			QColor color,
			const QString& bottom_title,
			const QString& left_title,
			bool draw_outlier = true,
			int outlier_scatter_size = 1);
	};


};

