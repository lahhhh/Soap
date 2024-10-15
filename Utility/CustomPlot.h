#pragma once

#include "Identifier.h"

#include <qCustomPlot.h>

#include "Custom.h"
#include "CustomPlotObject.h"
#include "GraphSettings.h"
#include "CustomPlotPatch.h"
#include "CustomPlotUtility.h"
#include "CustomPlotColor.h"

#include "FeatureHandler.h"

namespace custom_plot {

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> prepare(const GraphSettings& gs);

	std::pair<QCustomPlot*, QCPAxisRect*> prepare_ar(const GraphSettings& gs);

	std::pair<QCustomPlot*, QCPLayoutGrid*> prepare_lg(const GraphSettings& gs);

	std::tuple<QCustomPlot*, QCPLayoutGrid*, QCPLayoutGrid*> prepare_lg_lg(const GraphSettings& gs);

	void set_scatter_plot_axis_style(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		const QString& bottom_title,
		const QString& left_title,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const GraphSettings& gs);

	QCustomPlot* initialize_plot(const GraphSettings& gs);

	void set_axis_label(
		QCPAxis::AxisType type,
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	);

	void set_left_axis_label(
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	);

	void set_top_axis_label(
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	);

	void set_bottom_axis_label(
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	);

	void set_left_title(QCPAxisRect* axis_rect, const QString& label, const GraphSettings& gs, bool use_setting_title = false);

	void set_bottom_title(QCPAxisRect* axis_rect, const QString& label, const GraphSettings& gs, bool use_setting_title = false);

	void set_simple_axis(
		QCPAxisRect* axis_rect,
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs,
		bool use_setting_title = true
	);

	void set_simple_axis_no_title(QCPAxisRect* axis_rect, const GraphSettings& gs);

	void add_title(QCustomPlot* draw_area, const QString& title, const GraphSettings& gs);

	void set_arrow_axis(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		double margin,
		int type,
		const GraphSettings& gs
	);

	void add_gradient_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		double minimum_value,
		double maximum_value,
		const QString& title,
		const GraphSettings& gs,
		QColor low_color,
		QColor middle_color,
		QColor high_color
	);

	void add_gradient_legend(
		QCustomPlot* draw_area, 
		QCPLayoutGrid* legend_layout,
		double minimum_value, 
		double maximum_value,
		const QString& title, 
		const GraphSettings& gs
	);

	void bar_plot(
		QCustomPlot* draw_area, 
		QCPAxisRect* axis_rect, 
		QColor color,
		const Eigen::ArrayXd& bar_location, 
		const Eigen::ArrayXd& bar_length, 
		int bar_width, 
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs,
		bool vertical = true
	);

	void bar_plot_enrichment(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& bar_location,
		const Eigen::ArrayXd& bar_length,
		const Eigen::ArrayXd& values,
		int bar_width,
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs,
		bool vertical = true
	);

	void proportion_legend(
		QCustomPlot* draw_area, 
		QCPLayoutGrid* legend_layout,
		const QString& legend_title, 
		const GraphSettings& gs
	);

	void set_heatmap_color(
		QCustomPlot* draw_area,
		const GraphSettings& gs,
		QCPColorMap* heatmap,
		double min_val,
		double max_val);

	QCustomPlot* heatmap_plot(
		const HEATMAP_PLOT_ELEMENTS& elements, 
		const GraphSettings& gs);

	void heatmap_plot(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		QCPLayoutGrid* legend_layout,
		int cell_width,
		int cell_height,
		int border_width,
		const QStringList& left_labels,
		const QStringList& bottom_labels,
		const Eigen::MatrixXd& values,
		bool normalize_by_row,
		const QString& group_legend_title,
		bool annotation_at_top,
		const GraphSettings& gs
	);

	void heatmap_plot2(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		int cell_width,
		int cell_height,
		int border_width,
		const QStringList& left_labels,
		const QStringList& bottom_labels,
		const Eigen::MatrixXd& values,
		const GraphSettings& gs
	);

	void bubble_plot(
		QCustomPlot* draw_area, 
		QCPAxisRect* axis_rect, 
		QCPLayoutGrid* legend_layout, 
		const QStringList& left_labels, 
		const QStringList& bottom_labels, 
		const Eigen::MatrixXd& values, 
		const Eigen::MatrixXi& dot_size,
		const QString& value_legend_title, 
		const QString& size_legend_title,
		const GraphSettings& gs
	);

	void histogram_plot(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		const QVector<double>& x,
		int unit,
		QColor color,
		int width,
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs);

	void add_square_legend(
		QCustomPlot* draw_area, 
		QCPLayoutGrid* legend_layout,
		const QStringList& levels, 
		const QVector<QColor>& colors,
		QFont font, 
		const QString& legend_title,
		QFont legend_title_font, 
		int legend_column_width,
		int legend_row_width
	);

	void add_round_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		const QStringList& levels,
		const QList<QColor>& colors,
		const QString& legend_title,
		const GraphSettings& gs,
		int legend_index = 0
	);

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> embedding_single_color_plot(
		const QString& title,
		const Eigen::MatrixXd& embedding_matrix,
		const QColor& color,
		const QStringList& embedding_names,
		const GraphSettings& graph_settings
	);

	QCustomPlot* coverage_plot(COVERAGE_PLOT_ELEMENTS res, const GraphSettings& gs);

	std::pair<QCPAxisRect*, QCPLayoutGrid*> __feature_plot(
		QCustomPlot* draw_area,
		QCPLayoutGrid* sub_layout,
		const QUERY_DATA& data,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		bool scale,
		const GraphSettings& gs);

	QCustomPlot* feature_plot(
		const QList<QUERY_DATA>& data,
		Embedding* embedding,
		bool scale, 
		int nrow, 
		const GraphSettings& gs);

	QCustomPlot* feature_plot(
		const QList<QUERY_DATA>& data,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		bool scale,
		int nrow,
		const GraphSettings& gs);

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> feature_plot(
		const QUERY_DATA& data,
		Embedding* embedding,
		bool scale,
		const GraphSettings& gs);

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> feature_plot(
		const QUERY_DATA& data,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		bool scale,
		const GraphSettings& gs);

	QCustomPlot* monocle3_feature_plot(
		const QList<QUERY_DATA>& data,
		const Eigen::ArrayXd& x,
		const QUERY_DATA& f,
		const QString& bottom_title,
		int nrow,
		const GraphSettings& gs);
}
