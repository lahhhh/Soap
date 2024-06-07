#pragma once

#include "Identifier.h"

#include "qCustomPlot.h"

#include "GraphSettingsMacro.h"

class GraphSettings
{

public:

	GraphSettings();
	GraphSettings(const GraphSettings&) = default;
	GraphSettings(GraphSettings&&) = default;
	GraphSettings& operator=(const GraphSettings&) = default;
	GraphSettings& operator=(GraphSettings&&) = default;
	~GraphSettings() = default;

	G_SET_IDENTIFIER("Graph Settings");

	QMap<QString, QString> string_info_;

	QMap<QString, QStringList> string_list_info_;

	QMap<QString, QFont> font_info_;

	QMap<QString, int> integer_info_;

	QMap<QString, bool> bool_info_;

	QMap<QString, QColor> color_info_;

	using Palette = QMap<QString, QColor>;
	Palette palette_;

	enum class AxisStyle : int {
		Simple = 0,
		NoAxis = 1,
		Arrow = 2,
		Arrow2 = 3
	};

	void set_activated(bool activated);
	bool active() const;

	QString get_title(const QString& title) const;
	QFont get_title_font(const QFont& title_font = DEFAULT_TITLE_FONT) const;

	AxisStyle get_axis_style(AxisStyle style = DEFAULT_AXIS_STYLE) const;
		
	QString get_bottom_title(const QString& title) const;
	QString get_left_title(const QString& title) const;
	QFont get_bottom_title_font(const QFont& bottom_title_font = DEFAULT_BOTTOM_TITLE_FONT) const;
	QFont get_left_title_font(const QFont& left_title_font = DEFAULT_LEFT_TITLE_FONT) const;

	QFont get_bottom_label_font(const QFont& bottom_label_font = DEFAULT_BOTTOM_LABEL_FONT) const;
	QFont get_left_label_font(const QFont& left_label_font = DEFAULT_LEFT_LABEL_FONT) const;
	QFont get_label_font(QCPAxis::AxisType type) const;
	int get_bottom_label_angle(const int bottom_label_angle = DEFAULT_BOTTOM_LABEL_ANGLE) const;
	int get_left_label_angle(const int left_label_angle = DEFAULT_LEFT_LABEL_ANGLE) const;
	int get_label_angle(QCPAxis::AxisType type) const;

	QColor get_gradient_low_color(const QColor& color = DEFAULT_GRADIENT_LOW_COLOR) const;
	QColor get_gradient_middle_color(const QColor& color = DEFAULT_GRADIENT_MIDDLE_COLOR) const;
	QColor get_gradient_high_color(const QColor& color = DEFAULT_GRADIENT_HIGH_COLOR) const;

	QString get_legend_title(const QString& title, unsigned int index = 0) const;
	QFont get_legend_title_font(const QFont& legend_title_font = DEFAULT_LEGEND_TITLE_FONT) const;
	QFont get_legend_label_font(const QFont& legend_label_font = DEFAULT_LEGEND_LABEL_FONT) const;
	int get_legend_column_width(const int legend_column_width = DEFAULT_LEGEND_COLUMN_WIDTH) const;
	int get_legend_row_width(const int legend_row_width = DEFAULT_LEGEND_ROW_WIDTH) const;

	int get_scatter_point_size(int scatter_size = DEFAULT_SCATTER_POINT_SIZE) const;
	bool is_scatter_labeled() const;
	QFont get_scatter_label_font() const;

	bool is_legend_tick_shown() const;
	bool is_transparent_background() const;
	bool use_facet_violin_plot() const;
	bool use_boxplot() const;
	bool boxplot_draw_outlier() const;

	QFont get_legend_tick_label_font(const QFont& legend_tick_label_font = DEFAULT_LEGEND_TICK_LABEL_FONT) const;

	void apply_palette(const QStringList& names, QVector<QColor>& colors) const;
	void apply_palette(const QString& name, QColor& color) const;

	QList<QColor> palette(const QStringList& levels) const;
	QMap<QString, QColor> palette_map(const QStringList& levels) const;
};
