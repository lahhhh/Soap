#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include "Switch.h"
#include "MultipleLineEditLayout.h"

#include "GraphSettings.h"

class GraphSettingDialog 
	: public QDialog
{
	Q_OBJECT
public:
	GraphSettingDialog(const GraphSettings& gs);

	static void set_graph_setting(GraphSettings& gs);

	bool is_accepted_ = false;

	GraphSettings gs_;

	QVBoxLayout* full_layout_;

	QHBoxLayout* main_layout_;

	QLabel* title_area_label_;

	QLabel* title_label_;
	QLineEdit* title_line_edit_;
	Switch* title_active_switch_;

	QLabel* title_font_label_;
	QLabel* title_font_view_label_;
	QPushButton* title_font_button_;
	Switch* title_font_active_switch_;

	QLabel* axis_area_label_;

	QLabel* axis_style_label_;
	QComboBox* axis_style_box_;
	Switch* axis_style_active_switch_;

	QLabel* bottom_axis_area_label_;

	QLabel* bottom_title_label_;
	QLineEdit* bottom_title_line_edit_;
	Switch* bottom_title_active_switch_;

	QLabel* bottom_title_font_label_;
	QLabel* bottom_title_font_view_label_;
	QPushButton* bottom_title_font_button_;
	Switch* bottom_title_font_active_switch_;

	QLabel* bottom_label_font_label_;
	QLabel* bottom_label_font_view_label_;
	QPushButton* bottom_label_font_button_;
	Switch* bottom_label_font_active_switch_;

	QLabel* bottom_label_angle_label_;
	QLineEdit* bottom_label_angle_line_edit_;
	Switch* bottom_label_angle_active_switch_;

	QLabel* left_axis_area_label_;

	QLabel* left_title_label_;
	QLineEdit* left_title_line_edit_;
	Switch* left_title_active_switch_;

	QLabel* left_title_font_label_;
	QLabel* left_title_font_view_label_;
	QPushButton* left_title_font_button_;
	Switch* left_title_font_active_switch_;

	QLabel* left_label_font_label_;
	QLabel* left_label_font_view_label_;
	QPushButton* left_label_font_button_;
	Switch* left_label_font_active_switch_;

	QLabel* left_label_angle_label_;
	QLineEdit* left_label_angle_line_edit_;
	Switch* left_label_angle_active_switch_;

	QLabel* legend_area_label_;

	MultipleLineEditLayout* legend_title_layout_;
	Switch* legend_title_active_switch_;

	QLabel* legend_title_font_label_;
	QLabel* legend_title_font_view_label_;
	QPushButton* legend_title_font_button_;
	Switch* legend_title_font_active_switch_;

	QLabel* legend_label_font_label_;
	QLabel* legend_label_font_view_label_;
	QPushButton* legend_label_font_button_;
	Switch* legend_label_font_active_switch_;

	QLabel* legend_tick_label_font_label_;
	QLabel* legend_tick_label_font_view_label_;
	QPushButton* legend_tick_label_font_button_;
	Switch* legend_tick_label_font_active_switch_;

	QLabel* legend_column_width_label_;
	QLineEdit* legend_column_width_line_edit_;
	Switch* legend_column_width_active_switch_;

	QLabel* legend_row_width_label_;
	QLineEdit* legend_row_width_line_edit_;
	Switch* legend_row_width_active_switch_;

	QLabel* show_legend_tick_label_;
	Switch* show_legend_tick_switch_;
	Switch* show_legend_tick_active_switch_;

	QLabel* gradient_area_label_;

	QLabel* gradient_low_color_label_;
	QLabel* gradient_low_color_view_label_;
	QPushButton* gradient_low_color_button_;
	Switch* gradient_low_color_active_switch_;

	QLabel* gradient_middle_color_label_;
	QLabel* gradient_middle_color_view_label_;
	QPushButton* gradient_middle_color_button_;
	Switch* gradient_middle_color_active_switch_;

	QLabel* gradient_high_color_label_;
	QLabel* gradient_high_color_view_label_;
	QPushButton* gradient_high_color_button_;
	Switch* gradient_high_color_active_switch_;

	QLabel* scatter_area_label_;

	QLabel* scatter_labeled_label_;
	Switch* scatter_labeled_switch_;
	Switch* scatter_labeled_active_switch_;

	QLabel* scatter_point_size_label_;
	QLineEdit* scatter_point_size_line_edit_;
	Switch* scatter_point_size_active_switch_;

	QLabel* scatter_label_font_label_;
	QLabel* scatter_label_font_view_label_;
	QPushButton* scatter_label_font_button_;
	Switch* scatter_label_font_active_switch_;

	QLabel* others_area_label_;

	QPushButton* palette_settings_button_;
	Switch* palette_active_switch_;
	QPushButton* export_settings_button_;
	QPushButton* import_settings_button_;

	QLabel* transparent_background_label_;
	Switch* transparent_background_switch_;
	Switch* transparent_background_active_switch_;

	QLabel* use_facet_violin_plot_label_;
	Switch* use_facet_violin_plot_switch_;
	Switch* use_facet_violin_plot_active_switch_;

	QLabel* use_boxplot_label_;
	Switch* use_boxplot_switch_;
	Switch* use_boxplot_active_switch_;

	QLabel* boxplot_draw_outlier_label_;
	Switch* boxplot_draw_outlier_switch_;
	Switch* boxplot_draw_outlier_active_switch_;

	QPushButton* finish_button_;
	QPushButton* cancel_button_;

	void refresh_interface();

private slots:

	void s_set_title_font();

	void s_set_bottom_title_font();

	void s_set_left_title_font();

	void s_set_bottom_axis_label_font();

	void s_set_left_axis_label_font();

	void s_set_legend_title_font();

	void s_set_legend_label_font();

	void s_set_legend_tick_label_font();

	void s_set_scatter_label_font();

	void s_set_gradient_low_color();

	void s_set_gradient_middle_color();

	void s_set_gradient_high_color();

	void s_set_palette();

	void s_export_settings();
	void s_import_settings();

	void accept();
	void reject();

};
