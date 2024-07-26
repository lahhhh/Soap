#include "GraphSettingDialog.h"

#include "Identifier.h"
#include "ItemDatabase.h"
#include "MessageDialog.h"
#include "PaletteSettingDialog.h"
#include "CommonDialog.h"

#include "SoapGUI.h"

inline QString font2string(QFont font) {
	return QString("%1 %2").arg(font.family()).arg(font.pointSize());
}

GraphSettingDialog::GraphSettingDialog(const GraphSettings& gs) :
	gs_(gs)
{
	this->full_layout_ = new QVBoxLayout;
	this->main_layout_ = new QHBoxLayout;
	QVBoxLayout* column_layout = new QVBoxLayout;
	QHBoxLayout* row_layout;

	G_SET_LABEL(this->title_area_label_, "Title", soap::LargeSize);
	G_SET_LARGE_LABEL(this->title_area_label_);
	column_layout->addWidget(this->title_area_label_);

	G_SET_LABEL(this->title_label_, "title", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(this->title_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(this->title_active_switch_, false, title_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, title_label_, title_line_edit_, title_active_switch_);

	G_SET_LABEL(this->title_font_label_, "font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->title_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(this->title_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(this->title_font_active_switch_, false, title_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, title_font_label_, title_font_view_label_, title_font_button_, title_font_active_switch_);
	connect(this->title_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_title_font);

	G_SET_LABEL(axis_area_label_, "Axis", soap::LargeSize);
	G_SET_LARGE_LABEL(axis_area_label_);
	column_layout->addWidget(axis_area_label_);

	G_SET_LABEL_FIXED_HEIGHT(axis_style_label_, "style (scatter plot only)", 30);
	G_SET_COMBOBOX(axis_style_box_, QStringList() << "Simple" << "No Axis" << "Arrow" << "Arrow2", 30);
	G_SET_SWITCH(axis_style_active_switch_, false, axis_style_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, axis_style_label_, axis_style_box_, axis_style_active_switch_);

	G_SET_LABEL(bottom_axis_area_label_, "Bottom Axis", soap::LargeSize);
	G_SET_LARGE_LABEL(bottom_axis_area_label_);
	column_layout->addWidget(bottom_axis_area_label_);

	G_SET_LABEL(bottom_title_label_, "title", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(bottom_title_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(bottom_title_active_switch_, false, bottom_title_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, bottom_title_label_, bottom_title_line_edit_, bottom_title_active_switch_);

	G_SET_LABEL(bottom_title_font_label_, "title font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(bottom_title_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(bottom_title_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(bottom_title_font_active_switch_, false, bottom_title_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, bottom_title_font_label_, bottom_title_font_view_label_, bottom_title_font_button_, bottom_title_font_active_switch_);
	connect(bottom_title_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_bottom_title_font);

	G_SET_LABEL(bottom_label_font_label_, "label font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(bottom_label_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(bottom_label_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(bottom_label_font_active_switch_, false, bottom_label_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, bottom_label_font_label_, bottom_label_font_view_label_, bottom_label_font_button_, bottom_label_font_active_switch_);
	connect(bottom_label_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_bottom_axis_label_font);

	G_SET_LABEL(bottom_label_angle_label_, "label angle", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(bottom_label_angle_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(bottom_label_angle_active_switch_, false, bottom_label_angle_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, bottom_label_angle_label_, bottom_label_angle_line_edit_, bottom_label_angle_active_switch_);

	G_SET_LABEL(left_axis_area_label_, "Left Axis", soap::MiddleSize);
	G_SET_LARGE_LABEL(left_axis_area_label_);
	column_layout->addWidget(left_axis_area_label_);

	G_SET_LABEL(left_title_label_, "title", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(left_title_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(left_title_active_switch_, false, left_title_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, left_title_label_, left_title_line_edit_, left_title_active_switch_);

	G_SET_LABEL(left_title_font_label_, "title font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(left_title_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(left_title_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(left_title_font_active_switch_, false, left_title_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, left_title_font_label_, left_title_font_view_label_, left_title_font_button_, left_title_font_active_switch_);
	connect(left_title_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_left_title_font);

	G_SET_LABEL(left_label_font_label_, "label font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(left_label_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(left_label_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(left_label_font_active_switch_, false, left_label_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, left_label_font_label_, left_label_font_view_label_, left_label_font_button_, left_label_font_active_switch_);
	connect(left_label_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_left_axis_label_font);

	G_SET_LABEL(left_label_angle_label_, "label angle", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(left_label_angle_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(left_label_angle_active_switch_, false, left_label_angle_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, left_label_angle_label_, left_label_angle_line_edit_, left_label_angle_active_switch_);

	G_SET_LABEL(legend_area_label_, "Legend", soap::LargeSize);
	G_SET_LARGE_LABEL(legend_area_label_);
	column_layout->addWidget(legend_area_label_);

	this->legend_title_layout_ = new MultipleLineEditLayout("Title", nullptr, {}, this);
	G_SET_SWITCH(legend_title_active_switch_, false, static_cast<QLabel*>(nullptr), soap::MiddleSize);
	G_ADD_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, legend_title_layout_, legend_title_active_switch_);

	G_SET_LABEL(legend_title_font_label_, "title font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(legend_title_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(legend_title_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(legend_title_font_active_switch_, false, legend_title_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, legend_title_font_label_, legend_title_font_view_label_, legend_title_font_button_, legend_title_font_active_switch_);
	connect(legend_title_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_legend_title_font);

	G_SET_LABEL(legend_label_font_label_, "label font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(legend_label_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(legend_label_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(legend_label_font_active_switch_, false, legend_label_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, legend_label_font_label_, legend_label_font_view_label_, legend_label_font_button_, legend_label_font_active_switch_);
	connect(legend_label_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_legend_label_font);

	G_SET_LABEL(legend_tick_label_font_label_, "tick label font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(legend_tick_label_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(legend_tick_label_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(legend_tick_label_font_active_switch_, false, legend_tick_label_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, legend_tick_label_font_label_, legend_tick_label_font_view_label_, legend_tick_label_font_button_, legend_tick_label_font_active_switch_);
	connect(legend_tick_label_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_legend_tick_label_font);

	G_SET_LABEL(legend_column_width_label_, "column width", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(legend_column_width_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(legend_column_width_active_switch_, false, legend_column_width_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, legend_column_width_label_, legend_column_width_line_edit_, legend_column_width_active_switch_);
	legend_column_width_line_edit_->setValidator(new QIntValidator(this));

	G_SET_LABEL(legend_row_width_label_, "row width", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(legend_row_width_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(legend_row_width_active_switch_, false, legend_row_width_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, legend_row_width_label_, legend_row_width_line_edit_, legend_row_width_active_switch_);
	legend_row_width_line_edit_->setValidator(new QIntValidator(this));

	G_SET_LABEL(show_legend_tick_label_, "show legend ticks", soap::MiddleSize);
	G_SET_SWITCH_NO_LINK(show_legend_tick_switch_, false, soap::MiddleSize);
	G_SET_SWITCH(show_legend_tick_active_switch_, false, show_legend_tick_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, show_legend_tick_label_, show_legend_tick_switch_, show_legend_tick_active_switch_);

	this->main_layout_->addLayout(column_layout);

	column_layout = new QVBoxLayout;

	G_SET_LABEL(gradient_area_label_, "Gradient", soap::LargeSize);
	G_SET_LARGE_LABEL(gradient_area_label_);
	column_layout->addWidget(gradient_area_label_);

	G_SET_LABEL(gradient_low_color_label_, "gradient low color", soap::MiddleSize);
	G_SET_EMPTY_LABEL(gradient_low_color_view_label_, soap::MiddleSize);
	G_SET_BUTTON(gradient_low_color_button_, "choose color", soap::MiddleSize);
	G_SET_SWITCH(gradient_low_color_active_switch_, false, gradient_low_color_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, gradient_low_color_label_, gradient_low_color_view_label_, gradient_low_color_button_, gradient_low_color_active_switch_);
	connect(gradient_low_color_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_gradient_low_color);

	G_SET_LABEL(gradient_middle_color_label_, "gradient middle color", soap::MiddleSize);
	G_SET_EMPTY_LABEL(gradient_middle_color_view_label_, soap::MiddleSize);
	G_SET_BUTTON(gradient_middle_color_button_, "choose color", soap::MiddleSize);
	G_SET_SWITCH(gradient_middle_color_active_switch_, false, gradient_middle_color_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, gradient_middle_color_label_, gradient_middle_color_view_label_, gradient_middle_color_button_, gradient_middle_color_active_switch_);
	connect(gradient_middle_color_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_gradient_middle_color);

	G_SET_LABEL(gradient_high_color_label_, "gradient high color", soap::MiddleSize);
	G_SET_EMPTY_LABEL(gradient_high_color_view_label_, soap::MiddleSize);
	G_SET_BUTTON(gradient_high_color_button_, "choose color", soap::MiddleSize);
	G_SET_SWITCH(gradient_high_color_active_switch_, false, gradient_high_color_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, gradient_high_color_label_, gradient_high_color_view_label_, gradient_high_color_button_, gradient_high_color_active_switch_);
	connect(gradient_high_color_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_gradient_high_color);

	G_SET_LABEL(scatter_area_label_, "Scatter", soap::LargeSize);
	G_SET_LARGE_LABEL(scatter_area_label_);
	column_layout->addWidget(scatter_area_label_);

	G_SET_LABEL(scatter_labeled_label_, "scatter label", soap::MiddleSize);
	G_SET_SWITCH_NO_LINK(scatter_labeled_switch_, false, soap::MiddleSize);
	G_SET_SWITCH(scatter_labeled_active_switch_, false, scatter_labeled_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, scatter_labeled_label_, scatter_labeled_switch_, scatter_labeled_active_switch_);

	G_SET_LABEL(scatter_point_size_label_, "point size", soap::MiddleSize);
	G_SET_EMPTY_LINEEDIT(scatter_point_size_line_edit_, soap::MiddleSize);
	G_SET_SWITCH(scatter_point_size_active_switch_, false, scatter_point_size_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, scatter_point_size_label_, scatter_point_size_line_edit_, scatter_point_size_active_switch_);
	scatter_point_size_line_edit_->setValidator(new QIntValidator(this));

	G_SET_LABEL(scatter_label_font_label_, "label font", soap::MiddleSize);
	G_SET_EMPTY_LABEL(scatter_label_font_view_label_, soap::MiddleSize);
	G_SET_BUTTON(scatter_label_font_button_, "choose font", soap::MiddleSize);
	G_SET_SWITCH(scatter_label_font_active_switch_, false, scatter_label_font_label_, soap::MiddleSize);
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, scatter_label_font_label_, scatter_label_font_view_label_, scatter_label_font_button_, scatter_label_font_active_switch_);

	connect(scatter_label_font_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_scatter_label_font);

	G_SET_LABEL(others_area_label_, "Other Settings", soap::LargeSize);
	G_SET_LARGE_LABEL(others_area_label_);
	column_layout->addWidget(others_area_label_);

	G_SET_BUTTON_FINISH_STYLE(this->palette_settings_button_, "Palette Settings");
	G_SET_SWITCH(palette_active_switch_, false, this->palette_settings_button_, soap::MiddleSize);
	G_SET_BUTTON_FINISH_STYLE(this->import_settings_button_, "Import Settings");
	G_SET_BUTTON_FINISH_STYLE(this->export_settings_button_, "Export Settings");
	G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, this->palette_settings_button_, this->palette_active_switch_, this->import_settings_button_, this->export_settings_button_);
	connect(this->palette_settings_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_set_palette);
	connect(this->import_settings_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_import_settings);
	connect(this->export_settings_button_, &QPushButton::clicked, this, &GraphSettingDialog::s_export_settings);

	G_SET_LABEL_FIXED_HEIGHT(transparent_background_label_, "use transparent background", 30);
	G_SET_SWITCH_NO_LINK(transparent_background_switch_, false, soap::MiddleSize);
	G_SET_SWITCH(transparent_background_active_switch_, false, transparent_background_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, transparent_background_label_, transparent_background_switch_, transparent_background_active_switch_);

	G_SET_LABEL_FIXED_HEIGHT(use_facet_violin_plot_label_, "use facet violin plot", 30);
	G_SET_SWITCH_NO_LINK(use_facet_violin_plot_switch_, false, soap::MiddleSize);
	G_SET_SWITCH(use_facet_violin_plot_active_switch_, false, use_facet_violin_plot_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, use_facet_violin_plot_label_, use_facet_violin_plot_switch_, use_facet_violin_plot_active_switch_);

	G_SET_LABEL(use_boxplot_label_, "use boxplot", soap::MiddleSize);
	G_SET_SWITCH_NO_LINK(use_boxplot_switch_, false, soap::MiddleSize);
	G_SET_SWITCH(use_boxplot_active_switch_, false, use_boxplot_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, use_boxplot_label_, use_boxplot_switch_, use_boxplot_active_switch_);

	G_SET_LABEL_FIXED_HEIGHT(boxplot_draw_outlier_label_, "draw outlier (Boxplot)", 30);
	G_SET_SWITCH_NO_LINK(boxplot_draw_outlier_switch_, false, soap::MiddleSize);
	G_SET_SWITCH(boxplot_draw_outlier_active_switch_, false, boxplot_draw_outlier_label_, soap::MiddleSize);
	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(column_layout, row_layout, boxplot_draw_outlier_label_, boxplot_draw_outlier_switch_, boxplot_draw_outlier_active_switch_);


	this->main_layout_->addLayout(column_layout);

	this->full_layout_->addLayout(this->main_layout_);

	G_SET_FINISH_BUTTON;

	G_SET_CANCEL_BUTTON;

	G_ADD_DOUBLE_ITEM_ROWLAYOUT(this->full_layout_, row_layout, this->finish_button_, this->cancel_button_);

	connect(this->finish_button_, &QPushButton::clicked, this, &GraphSettingDialog::accept);
	connect(this->cancel_button_, &QPushButton::clicked, this, &GraphSettingDialog::reject);

	setLayout(this->full_layout_);

	this->refresh_interface();

	G_SET_ICON;
	this->setWindowTitle("Graph Settings");
	this->exec();

}

void GraphSettingDialog::s_set_palette() {
	PaletteSettingDialog::get_response(this->gs_.palette_);
};

void GraphSettingDialog::s_export_settings() {
	QString file_path = QFileDialog::getSaveFileName(this, "Set File Name for Graph Settings", "", "soap item(*.sif)");
	if (file_path.isEmpty())return;

	bool success = ItemDatabase::write_item("Graph Settings", this->gs_, file_path);
	if (!success) {
		MessageDialog::get_response("Error", "Graph Settings Exporting Failed.");
	}
};

void GraphSettingDialog::s_import_settings() {
	QString file_path = QFileDialog::getOpenFileName(this, "Select Settings", "", "soap item(*.sif)");
	if (file_path.isEmpty())return;

	try {
		bool success = ItemDatabase::read_item(file_path, this->gs_);
		if (!success) {
			MessageDialog::get_response("Error", "Graph Settings Reading Failed.");
			return;
		}
	}
	catch (...) {
		MessageDialog::get_response("Error", "Graph Settings Reading Failed.");
		return;
	}
	this->refresh_interface();
};

void GraphSettingDialog::refresh_interface() {
	this->title_line_edit_->setText(this->gs_.string_info_[TITLE]);
	this->title_active_switch_->setChecked(this->gs_.bool_info_[IS_TITLE_ACTIVE]);

	this->title_font_view_label_->setText(font2string(this->gs_.font_info_.value(TITLE_FONT, DEFAULT_TITLE_FONT)));
	this->title_font_active_switch_->setChecked(this->gs_.bool_info_[IS_TITLE_FONT_ACTIVE]);

	int axis_style = this->gs_.integer_info_[AXIS_STYLE];
	this->axis_style_box_->setCurrentIndex(axis_style);
	this->axis_style_active_switch_->setChecked(this->gs_.bool_info_[IS_AXIS_STYLE_ACTIVE]);

	this->bottom_title_line_edit_->setText(this->gs_.string_info_[BOTTOM_TITLE]);
	this->bottom_title_active_switch_->setChecked(this->gs_.bool_info_[IS_BOTTOM_TITLE_ACTIVE]);

	this->bottom_title_font_view_label_->setText(font2string(this->gs_.font_info_.value(BOTTOM_TITLE_FONT, DEFAULT_BOTTOM_TITLE_FONT)));
	this->bottom_title_font_active_switch_->setChecked(this->gs_.bool_info_[IS_BOTTOM_TITLE_FONT_ACTIVE]);

	this->bottom_label_font_view_label_->setText(font2string(this->gs_.font_info_.value(BOTTOM_LABEL_FONT, DEFAULT_BOTTOM_LABEL_FONT)));
	this->bottom_label_font_active_switch_->setChecked(this->gs_.bool_info_[IS_BOTTOM_LABEL_FONT_ACTIVE]);

	this->bottom_label_angle_line_edit_->setText(QString::number(this->gs_.integer_info_[BOTTOM_LABEL_ANGLE]));
	this->bottom_label_angle_active_switch_->setChecked(this->gs_.bool_info_[IS_BOTTOM_LABEL_ANGLE_ACTIVE]);

	this->left_title_line_edit_->setText(this->gs_.string_info_[LEFT_TITLE]);
	this->left_title_active_switch_->setChecked(this->gs_.bool_info_[IS_LEFT_TITLE_ACTIVE]);

	this->left_title_font_view_label_->setText(font2string(this->gs_.font_info_.value(LEFT_TITLE_FONT, DEFAULT_LEFT_TITLE_FONT)));
	this->left_title_font_active_switch_->setChecked(this->gs_.bool_info_[IS_LEFT_TITLE_FONT_ACTIVE]);

	this->left_label_font_view_label_->setText(font2string(this->gs_.font_info_.value(LEFT_LABEL_FONT, DEFAULT_LEFT_LABEL_FONT)));
	this->left_label_font_active_switch_->setChecked(this->gs_.bool_info_[IS_LEFT_LABEL_FONT_ACTIVE]);

	this->left_label_angle_line_edit_->setText(QString::number(this->gs_.integer_info_[LEFT_LABEL_ANGLE]));
	this->left_label_angle_active_switch_->setChecked(this->gs_.bool_info_[IS_LEFT_LABEL_ANGLE_ACTIVE]);

	this->legend_title_layout_->set_items(this->gs_.string_list_info_[LEGEND_TITLES]);
	this->legend_title_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_TITLES_ACTIVE]);

	this->legend_title_font_view_label_->setText(font2string(this->gs_.font_info_.value(LEGEND_TITLE_FONT, DEFAULT_LEGEND_TITLE_FONT)));
	this->legend_title_font_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_TITLE_FONT_ACTIVE]);

	this->legend_label_font_view_label_->setText(font2string(this->gs_.font_info_.value(LEGEND_LABEL_FONT, DEFAULT_LEGEND_LABEL_FONT)));
	this->legend_label_font_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_LABEL_FONT_ACTIVE]);

	this->legend_tick_label_font_view_label_->setText(font2string(this->gs_.font_info_.value(LEGEND_TICK_LABEL_FONT, DEFAULT_LEGEND_TICK_LABEL_FONT)));
	this->legend_tick_label_font_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_TICK_LABEL_FONT_ACTIVE]);

	this->legend_column_width_line_edit_->setText(QString::number(this->gs_.integer_info_[LEGEND_COLUMN_WIDTH]));
	this->legend_column_width_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_COLUMN_WIDTH_ACTIVE]);

	this->legend_row_width_line_edit_->setText(QString::number(this->gs_.integer_info_[LEGEND_ROW_WIDTH]));
	this->legend_row_width_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_ROW_WIDTH_ACTIVE]);

	this->show_legend_tick_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_TICK_SHOWN]);
	this->show_legend_tick_active_switch_->setChecked(this->gs_.bool_info_[IS_LEGEND_TICK_SHOWN_ACTIVE]);

	this->gradient_low_color_view_label_->setText(this->gs_.color_info_[GRADIENT_LOW_COLOR].name());
	G_SET_LABEL_COLOR(this->gradient_low_color_view_label_, this->gs_.color_info_.value(GRADIENT_LOW_COLOR, DEFAULT_GRADIENT_LOW_COLOR));
	this->gradient_low_color_active_switch_->setChecked(this->gs_.bool_info_[IS_GRADIENT_LOW_COLOR_ACTIVE]);

	this->gradient_middle_color_view_label_->setText(this->gs_.color_info_[GRADIENT_MIDDLE_COLOR].name());
	G_SET_LABEL_COLOR(this->gradient_middle_color_view_label_, this->gs_.color_info_.value(GRADIENT_MIDDLE_COLOR, DEFAULT_GRADIENT_MIDDLE_COLOR));
	this->gradient_middle_color_active_switch_->setChecked(this->gs_.bool_info_[IS_GRADIENT_MIDDLE_COLOR_ACTIVE]);

	this->gradient_high_color_view_label_->setText(this->gs_.color_info_[GRADIENT_HIGH_COLOR].name());
	G_SET_LABEL_COLOR(this->gradient_high_color_view_label_, this->gs_.color_info_.value(GRADIENT_HIGH_COLOR, DEFAULT_GRADIENT_HIGH_COLOR));
	this->gradient_high_color_active_switch_->setChecked(this->gs_.bool_info_[IS_GRADIENT_HIGH_COLOR_ACTIVE]);

	this->scatter_labeled_switch_->setChecked(this->gs_.bool_info_[IS_SCATTER_LABELED]);
	this->scatter_labeled_active_switch_->setChecked(this->gs_.bool_info_[IS_SCATTER_LABELED_ACTIVE]);

	this->scatter_point_size_line_edit_->setText(QString::number(this->gs_.integer_info_[SCATTER_POINT_SIZE]));
	this->scatter_point_size_active_switch_->setChecked(this->gs_.bool_info_[IS_SCATTER_POINT_SIZE_ACTIVE]);

	this->scatter_label_font_view_label_->setText(font2string(this->gs_.font_info_.value(SCATTER_LABEL_FONT, DEFAULT_SCATTER_LABEL_FONT)));
	this->scatter_label_font_active_switch_->setChecked(this->gs_.bool_info_[IS_SCATTER_LABEL_FONT_ACTIVE]);

	this->palette_active_switch_->setChecked(this->gs_.bool_info_[IS_PALATTE_ACTIVE]);

	this->transparent_background_switch_->setChecked(this->gs_.bool_info_[IS_TRANSPARENT_BACKGROUND]);
	this->transparent_background_active_switch_->setChecked(this->gs_.bool_info_[IS_TRANSPARENT_BACKGROUND_ACTIVE]);

	this->use_facet_violin_plot_switch_->setChecked(this->gs_.bool_info_[USE_FACET_VIOLIN_PLOT]);
	this->use_facet_violin_plot_active_switch_->setChecked(this->gs_.bool_info_[USE_FACET_VIOLIN_PLOT_ACTIVE]);

	this->use_boxplot_switch_->setChecked(this->gs_.bool_info_[USE_BOXPLOT]);
	this->use_boxplot_active_switch_->setChecked(this->gs_.bool_info_[USE_BOXPLOT_ACTIVE]);

	this->boxplot_draw_outlier_switch_->setChecked(this->gs_.bool_info_[BOXPLOT_DRAW_OUTLIER]);
	this->boxplot_draw_outlier_active_switch_->setChecked(this->gs_.bool_info_[BOXPLOT_DRAW_OUTLIER_ACTIVE]);
};

void GraphSettingDialog::s_set_gradient_low_color() {
	QColor color = QColorDialog::getColor(QColor(this->gradient_low_color_view_label_->text()), this, "Choose Color for Gradient Lower End");
	if (color.isValid()) {
		gradient_low_color_view_label_->setText(color.name());
		gradient_low_color_view_label_->setStyleSheet("QLabel{color:" + color.name() + "}");
	}
};

void GraphSettingDialog::s_set_gradient_middle_color() {
	QColor color = QColorDialog::getColor(QColor(this->gradient_middle_color_view_label_->text()), this, "Choose Color for Gradient Middle End");
	if (color.isValid()) {
		gradient_middle_color_view_label_->setText(color.name());
		gradient_middle_color_view_label_->setStyleSheet("QLabel{color:" + color.name() + "}");
	}
};

void GraphSettingDialog::s_set_gradient_high_color() {
	QColor color = QColorDialog::getColor(QColor(this->gradient_high_color_view_label_->text()), this, "Choose Color for Gradient Higher End");
	if (color.isValid()) {
		gradient_high_color_view_label_->setText(color.name());
		gradient_high_color_view_label_->setStyleSheet("QLabel{color:" + color.name() + "}");
	}
};

void GraphSettingDialog::s_set_title_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[TITLE_FONT] = font;
		title_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_bottom_title_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[BOTTOM_TITLE_FONT] = font;
		bottom_title_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_left_title_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[LEFT_TITLE_FONT] = font;
		left_title_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_bottom_axis_label_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[BOTTOM_LABEL_FONT] = font;
		bottom_label_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_left_axis_label_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[LEFT_LABEL_FONT] = font;
		left_label_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_legend_title_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[LEGEND_TITLE_FONT] = font;
		legend_title_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_legend_label_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[LEGEND_LABEL_FONT] = font;
		legend_label_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_legend_tick_label_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[LEGEND_TICK_LABEL_FONT] = font;
		legend_tick_label_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::s_set_scatter_label_font() {
	bool accepted = false;
	QFont font = QFontDialog::getFont(&accepted, this);
	if (accepted) {
		this->gs_.font_info_[SCATTER_LABEL_FONT] = font;
		scatter_label_font_view_label_->setText(font2string(font));
	}
};

void GraphSettingDialog::set_graph_setting(GraphSettings& gs) {
	GraphSettingDialog dlg(gs);
	if (!dlg.is_accepted_)return;

	gs = dlg.gs_; // assign font

	gs.string_info_[TITLE] = dlg.title_line_edit_->text();
	gs.bool_info_[IS_TITLE_ACTIVE] = dlg.title_active_switch_->value_;

	QStringList legend_titles = multiple_line_edit_to_list(dlg.legend_title_layout_->current_value());
	std::ranges::for_each(legend_titles, [](QString& title) { title = title.trimmed(); });
	gs.string_list_info_[LEGEND_TITLES] = legend_titles;
	gs.bool_info_[IS_LEGEND_TITLES_ACTIVE] = dlg.legend_title_active_switch_->value_;

	gs.string_info_[BOTTOM_TITLE] = dlg.bottom_title_line_edit_->text();
	gs.bool_info_[IS_BOTTOM_TITLE_ACTIVE] = dlg.bottom_title_active_switch_->value_;

	gs.string_info_[LEFT_TITLE] = dlg.left_title_line_edit_->text();
	gs.bool_info_[IS_LEFT_TITLE_ACTIVE] = dlg.left_title_active_switch_->value_;

	gs.bool_info_[IS_TITLE_FONT_ACTIVE] = dlg.title_font_active_switch_->value_;

	gs.integer_info_[AXIS_STYLE] = dlg.axis_style_box_->currentIndex();
	gs.bool_info_[IS_AXIS_STYLE_ACTIVE] = dlg.axis_style_active_switch_->value_;

	gs.bool_info_[IS_BOTTOM_LABEL_FONT_ACTIVE] = dlg.bottom_label_font_active_switch_->value_;

	gs.bool_info_[IS_LEFT_LABEL_FONT_ACTIVE] = dlg.left_label_font_active_switch_->value_;

	gs.integer_info_[BOTTOM_LABEL_ANGLE] = dlg.bottom_label_angle_line_edit_->text().toInt();
	gs.bool_info_[IS_BOTTOM_LABEL_ANGLE_ACTIVE] = dlg.bottom_label_angle_active_switch_->value_;

	gs.integer_info_[LEFT_LABEL_ANGLE] = dlg.left_label_angle_line_edit_->text().toInt();
	gs.bool_info_[IS_LEFT_LABEL_ANGLE_ACTIVE] = dlg.left_label_angle_active_switch_->value_;

	gs.bool_info_[IS_BOTTOM_TITLE_FONT_ACTIVE] = dlg.bottom_title_font_active_switch_->value_;

	gs.bool_info_[IS_LEFT_TITLE_FONT_ACTIVE] = dlg.left_title_font_active_switch_->value_;

	gs.bool_info_[IS_LEGEND_TICK_LABEL_FONT_ACTIVE] = dlg.legend_tick_label_font_active_switch_->value_;

	gs.color_info_[GRADIENT_LOW_COLOR] = QColor(dlg.gradient_low_color_view_label_->text());
	gs.bool_info_[IS_GRADIENT_LOW_COLOR_ACTIVE] = dlg.gradient_low_color_active_switch_->value_;

	gs.color_info_[GRADIENT_MIDDLE_COLOR] = QColor(dlg.gradient_middle_color_view_label_->text());
	gs.bool_info_[IS_GRADIENT_MIDDLE_COLOR_ACTIVE] = dlg.gradient_middle_color_active_switch_->value_;

	gs.color_info_[GRADIENT_HIGH_COLOR] = QColor(dlg.gradient_high_color_view_label_->text());
	gs.bool_info_[IS_GRADIENT_HIGH_COLOR_ACTIVE] = dlg.gradient_high_color_active_switch_->value_;

	gs.bool_info_[IS_LEGEND_TITLE_FONT_ACTIVE] = dlg.legend_title_font_active_switch_->value_;

	gs.bool_info_[IS_LEGEND_LABEL_FONT_ACTIVE] = dlg.legend_label_font_active_switch_->value_;

	gs.integer_info_[LEGEND_COLUMN_WIDTH] = dlg.legend_column_width_line_edit_->text().toInt();
	gs.bool_info_[IS_LEGEND_COLUMN_WIDTH_ACTIVE] = dlg.legend_column_width_active_switch_->value_;

	gs.integer_info_[LEGEND_ROW_WIDTH] = dlg.legend_row_width_line_edit_->text().toInt();
	gs.bool_info_[IS_LEGEND_ROW_WIDTH_ACTIVE] = dlg.legend_row_width_active_switch_->value_;

	gs.integer_info_[SCATTER_POINT_SIZE] = dlg.scatter_point_size_line_edit_->text().toInt();
	gs.bool_info_[IS_SCATTER_POINT_SIZE_ACTIVE] = dlg.scatter_point_size_active_switch_->value_;

	gs.bool_info_[IS_SCATTER_LABELED] = dlg.scatter_labeled_switch_->value_;
	gs.bool_info_[IS_SCATTER_LABELED_ACTIVE] = dlg.scatter_labeled_active_switch_->value_;

	gs.bool_info_[IS_SCATTER_LABEL_FONT_ACTIVE] = dlg.scatter_label_font_active_switch_->value_;

	gs.bool_info_[IS_LEGEND_TICK_SHOWN] = dlg.show_legend_tick_switch_->value_;
	gs.bool_info_[IS_LEGEND_TICK_SHOWN_ACTIVE] = dlg.show_legend_tick_active_switch_->value_;

	gs.bool_info_[IS_PALATTE_ACTIVE] = dlg.palette_active_switch_->value_;

	gs.bool_info_[IS_TRANSPARENT_BACKGROUND] = dlg.transparent_background_switch_->value_;
	gs.bool_info_[IS_TRANSPARENT_BACKGROUND_ACTIVE] = dlg.transparent_background_active_switch_->value_;

	gs.bool_info_[USE_FACET_VIOLIN_PLOT] = dlg.use_facet_violin_plot_switch_->value_;
	gs.bool_info_[USE_FACET_VIOLIN_PLOT_ACTIVE] = dlg.use_facet_violin_plot_active_switch_->value_;

	gs.bool_info_[USE_BOXPLOT] = dlg.use_boxplot_switch_->value_;
	gs.bool_info_[USE_BOXPLOT_ACTIVE] = dlg.use_boxplot_active_switch_->value_;

	gs.bool_info_[BOXPLOT_DRAW_OUTLIER] = dlg.boxplot_draw_outlier_switch_->value_;
	gs.bool_info_[BOXPLOT_DRAW_OUTLIER_ACTIVE] = dlg.boxplot_draw_outlier_active_switch_->value_;
};


void GraphSettingDialog::accept() {
	this->is_accepted_ = true;
	QDialog::accept();
};

void GraphSettingDialog::reject() {
	QDialog::reject();
};
