#include "MetadataViewWindow.h"
#include "CommonDialog.h"
#include "PlotWindow.h"
#include "CustomPlot.h"
#include "Identifier.h"
#include "GraphSettingDialog.h"

#include "SoapGUI.h"

MetadataViewWindow::MetadataViewWindow(
	SingleCellRna* single_cell_rna,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	signal_emitter_(signal_emitter),
	handler_(single_cell_rna)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &MetadataViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &MetadataViewWindow::s_check_data);

	set_layout();

	set_property();
}

MetadataViewWindow::MetadataViewWindow(
	SingleCellAtac* single_cell_atac,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	signal_emitter_(signal_emitter),
	handler_(single_cell_atac)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &MetadataViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &MetadataViewWindow::s_check_data);

	set_layout();

	set_property();
}

MetadataViewWindow::MetadataViewWindow(
	DataFrame* data_frame,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	signal_emitter_(signal_emitter),
	handler_(data_frame)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &MetadataViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &MetadataViewWindow::s_check_data);

	set_layout();

	set_property();
}

MetadataViewWindow::MetadataViewWindow(
	SingleCellMultiome* single_cell_multiome,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	signal_emitter_(signal_emitter),
	handler_(single_cell_multiome)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &MetadataViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &MetadataViewWindow::s_check_data);

	set_layout();

	set_property();
}

void MetadataViewWindow::s_check_data(void* data, soap::VariableType type, void* item) {

	if (data == this->handler_.data_) {
		this->close();
	}
};

void MetadataViewWindow::set_property() {
	this->setAttribute(Qt::WA_DeleteOnClose);
	G_SET_ICON;
	this->resize(1000, 800);
	this->setWindowTitle("Metadata View");
	this->show();
};

void MetadataViewWindow::set_layout() {

	auto feature_names = this->handler_.get_feature_names();

	this->valid_features_ = feature_names.all_names;

	G_SET_LABEL(this->feature_label_, "Feature", soap::MiddleSize);
	G_SET_LINEEDIT_WITH_COMPLETER(this->feature_line_edit_, "", this->valid_features_, soap::MiddleSize);

	G_SET_LABEL(this->normalize_label_, "Normalized", soap::MiddleSize);
	G_SET_SWITCH(this->normalize_switch_, true, this->normalize_label_, soap::MiddleSize);

	G_SET_LABEL(this->main_group_label_, "Main Group", soap::MiddleSize);
	G_SET_COMBOBOX(this->main_group_box_, QStringList() << "" << feature_names.factor_names, 30);

	G_SET_LABEL(this->sub_group_label_, "Sub Group", soap::MiddleSize);
	G_SET_COMBOBOX(this->sub_group_box_, QStringList() << "" << feature_names.factor_names, 30);

	G_SET_BUTTON(this->graph_setting_button_, "Graph Settings", soap::MiddleSize);
	G_SET_SWITCH(this->graph_setting_switch_, false, this->graph_setting_button_, soap::MiddleSize);
	connect(this->graph_setting_button_, &QPushButton::clicked, this, &MetadataViewWindow::s_set_graph_settings);

	G_SET_BUTTON(this->refresh_picture_button_, "Refresh", soap::MiddleSize);

	G_SET_PLOTSUITE_BUTTON;

	this->left_layout_ = new QVBoxLayout;
	this->left_layout_->addWidget(this->feature_label_);
	this->left_layout_->addWidget(this->feature_line_edit_);
	this->left_layout_->addWidget(this->normalize_label_);
	this->left_layout_->addWidget(this->normalize_switch_);

	if (this->handler_.type_ == FeatureHandler::DataType::SingleCellAtac || this->handler_.type_ == FeatureHandler::DataType::SingleCellMultiome) {
		G_SET_LABEL(this->gene_activity_label_, "Gene Activity", soap::MiddleSize);
		G_SET_SWITCH(this->gene_activity_switch_, false, this->gene_activity_label_, soap::MiddleSize);

		this->left_layout_->addWidget(this->gene_activity_label_);
		this->left_layout_->addWidget(this->gene_activity_switch_);
	}
	this->left_layout_->addWidget(this->main_group_label_);
	this->left_layout_->addWidget(this->main_group_box_);
	this->left_layout_->addWidget(this->sub_group_label_);
	this->left_layout_->addWidget(this->sub_group_box_);

	this->left_layout_->addStretch();

	this->left_layout_->addWidget(this->graph_setting_button_);
	this->left_layout_->addWidget(this->graph_setting_switch_);

	this->left_layout_->addStretch();

	this->left_layout_->addWidget(this->refresh_picture_button_);
	this->left_layout_->addWidget(this->previous_picture_button_);
	this->left_layout_->addWidget(this->next_picture_button_);
	this->left_layout_->addWidget(this->pop_picture_button_);
	this->left_layout_->addWidget(this->clear_picture_button_);
	this->left_layout_->addWidget(this->save_picture_button_);

	this->information_area_ = new InformationTextBrowser(this);
	this->information_area_->setFixedWidth(150);
	this->information_area_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);

	this->left_layout_->addStretch();
	this->left_layout_->addWidget(this->information_area_);

	this->menu_save_picture_ = new QMenu();
	this->action_save_png_ = new QAction("PNG", this);
	this->action_save_jpg_ = new QAction("JPG", this);
	this->action_save_bmp_ = new QAction("BMP", this);
	this->action_save_pdf_ = new QAction("PDF", this);
	this->action_save_pdf_and_png_ = new QAction("PDF and PNG", this);
	this->menu_save_picture_->addAction(this->action_save_png_);
	this->menu_save_picture_->addAction(this->action_save_jpg_);
	this->menu_save_picture_->addAction(this->action_save_bmp_);
	this->menu_save_picture_->addAction(this->action_save_pdf_);
	this->menu_save_picture_->addAction(this->action_save_pdf_and_png_);


	connect(this->action_save_jpg_, &QAction::triggered, this, &MetadataViewWindow::s_save_jpg);
	connect(this->action_save_png_, &QAction::triggered, this, &MetadataViewWindow::s_save_png);
	connect(this->action_save_bmp_, &QAction::triggered, this, &MetadataViewWindow::s_save_bmp);
	connect(this->action_save_pdf_, &QAction::triggered, this, &MetadataViewWindow::s_save_pdf);
	connect(this->action_save_pdf_and_png_, &QAction::triggered, this, &MetadataViewWindow::s_save_pdf_and_png);

	this->save_picture_button_->setMenu(this->menu_save_picture_);

	connect(this->feature_line_edit_, &QLineEdit::editingFinished, this, &MetadataViewWindow::s_check_feature_name);

	connect(this->refresh_picture_button_, &QPushButton::clicked, this, &MetadataViewWindow::s_refresh_plot);
	connect(this->previous_picture_button_, &QPushButton::clicked, this, &MetadataViewWindow::s_previous_plot);
	connect(this->next_picture_button_, &QPushButton::clicked, this, &MetadataViewWindow::s_next_plot);
	connect(this->clear_picture_button_, &QPushButton::clicked, this, &MetadataViewWindow::s_clear_plot);
	connect(this->pop_picture_button_, &QPushButton::clicked, this, &MetadataViewWindow::s_pop_plot);

	this->main_layout_ = new QHBoxLayout;
	this->main_layout_->addLayout(this->left_layout_);

	this->main_interface_ = new QWidget(this);
	this->main_interface_->setLayout(this->main_layout_);
	this->setCentralWidget(this->main_interface_);

	this->draw_suite_ = new PlotsSuite();
	connect(this->draw_suite_, &PlotsSuite::x_plot_prepared, this, &MetadataViewWindow::s_new_plot);
	connect(this->graph_setting_switch_, &Switch::toggled, this->draw_suite_, &PlotsSuite::s_setting_activate);
	this->draw_suite_->prepare();
};

void MetadataViewWindow::view(SingleCellRna* data, SignalEmitter* signal_emitter) {
	MetadataViewWindow* window = new MetadataViewWindow(data, signal_emitter);
};

void MetadataViewWindow::view(SingleCellAtac* data, SignalEmitter* signal_emitter) {
	MetadataViewWindow* window = new MetadataViewWindow(data, signal_emitter);
};

void MetadataViewWindow::view(SingleCellMultiome* data, SignalEmitter* signal_emitter) {
	MetadataViewWindow* window = new MetadataViewWindow(data, signal_emitter);
};

void MetadataViewWindow::view(DataFrame* data, SignalEmitter* signal_emitter) {
	MetadataViewWindow* window = new MetadataViewWindow(data, signal_emitter);
};

MetadataViewWindow::~MetadataViewWindow() {
	delete this->draw_suite_;
}

void MetadataViewWindow::s_check_feature_name() {
	QString currentFeature = this->feature_line_edit_->text();
	if (currentFeature.isEmpty()) {
		this->feature_line_edit_->setStyleSheet("border:0px");
	}
	else if (this->valid_features_.contains(currentFeature)) {
		this->feature_line_edit_->setStyleSheet("border:4px solid green");
	}
	else {
		this->feature_line_edit_->setStyleSheet("border:4px solid red");
	}
};

void MetadataViewWindow::s_set_graph_settings() {
	GraphSettingDialog::set_graph_setting(this->draw_suite_->graph_settings_);
	this->graph_setting_switch_->set_status(this->draw_suite_->graph_settings_.active());
};

void MetadataViewWindow::s_save_png() {

	QStringList picture_size = CommonDialog::get_response(
		nullptr,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "PNG(*.png)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->savePng(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MetadataViewWindow::s_save_bmp() {
	QStringList picture_size = CommonDialog::get_response(
		nullptr,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "BMP(*.bmp)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->saveBmp(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MetadataViewWindow::s_save_jpg() {
	QStringList picture_size = CommonDialog::get_response(
		nullptr,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "JPG(*.jpg)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->saveJpg(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MetadataViewWindow::s_save_pdf_and_png() {
	QStringList picture_size = CommonDialog::get_response(
		nullptr,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set PNG Name", "", "PNG(*.PNG)");
	if (picture_name.isEmpty())return;

	QString pdf_name = QFileDialog::getSaveFileName(this, "Set PDF Name", "", "PDF(*.pdf)");
	if (pdf_name.isEmpty())return;

	bool success = this->draw_area_->savePdf(pdf_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (!success) {
		G_LOG("PDF saving failed")
			return;
	}
	else {
		G_LOG("PDF saved.")
	}
	std::string pdf_file_name = pdf_name.toStdString();
	std::string picture_file_name = picture_name.toStdString();

	if (!custom::save_pdf_page_as_png(pdf_file_name, 0, picture_file_name)) {
		G_WARN("PNG Saving Failed.");
	}
	else {
		G_LOG("PNG saved.");
	}
}

void MetadataViewWindow::s_save_pdf() {
	QStringList picture_size = CommonDialog::get_response(
		nullptr,
		"PDF Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set PDF Name", "", "PDF(*.pdf)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->savePdf(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MetadataViewWindow::s_previous_plot() {
	if (this->draw_suite_->current_plot_id_ <= 1)return;
	this->main_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[-- this->draw_suite_->current_plot_id_];
	this->main_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
};

void MetadataViewWindow::s_next_plot() {
	if (this->draw_suite_->current_plot_id_ == this->draw_suite_->maximum_plot_id_)return;
	this->main_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[++ this->draw_suite_->current_plot_id_];
	this->main_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
};

void MetadataViewWindow::s_clear_plot() {
	this->main_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);
	this->draw_area_ = nullptr;
	this->draw_suite_->clear();
};

void MetadataViewWindow::s_pop_plot() {
	PlotWindow::show_plot(this->draw_suite_, "Figure", this->draw_area_->width(), this->draw_area_->height(), this);
};

void MetadataViewWindow::s_new_plot() {
	if (this->draw_area_ == nullptr) {
		this->draw_area_ = this->draw_suite_->current_plot_;
		this->main_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();
	}
	else {
		this->main_layout_->removeWidget(this->draw_area_);
		this->draw_area_->setVisible(false);
		this->draw_area_ = this->draw_suite_->current_plot_;
		this->main_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();
	}

};


void MetadataViewWindow::no_group_plot() {
	
	QString feature = this->feature_line_edit_->text();
	const auto& gs = this->draw_suite_->graph_settings_;

	bool normalize = this->normalize_switch_->value_;

	bool gene_activity{ false };
	if (this->gene_activity_switch_ != nullptr) {
		gene_activity = this->gene_activity_switch_->value_;
	}

	auto feature_data = this->handler_.get_data({ feature, normalize, gene_activity });

	if (feature_data.type == QUERY_DATA::DataType::notype) {
		G_WARN("No Data Found.");
		return;
	}

	auto colors = gs.palette({ feature });

	if (feature_data.type == QUERY_DATA::DataType::integer) {

		if (feature_data.dil.isEmpty()) {

			auto [imin, imax] = std::ranges::minmax(feature_data.di);
			double min = imin, max = imax;
			auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);

			if (gs.use_boxplot()) {
				custom_plot::patch::box(
					draw_area,
					axis_rect,
					colors[0],
					custom::cast<QVector>(custom::cast<double>(feature_data.di)),
					1.0,
					gs.boxplot_draw_outlier(),
					gs.get_scatter_point_size()
				);
			}
			else {
				std::tie(min, max) = custom_plot::patch::violin(
					draw_area,
					axis_rect,
					colors[0],
					custom::cast<QVector>(custom::cast<double>(feature_data.di)),
					1.0
				);
			}

			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::set_simple_axis_no_title(axis_rect, gs);
			custom_plot::patch::remove_bottom_ticks(axis_rect);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, 2), custom_plot::utility::get_range(min, max));
			custom_plot::add_title(draw_area, feature, gs);
			this->draw_suite_->update(draw_area);
		}
		else {
			auto distribution = custom::table(custom::cast<QString>(feature_data.di));
			QStringList levels = custom::cast<QString>(custom::sorted(feature_data.dil));

			auto colors = gs.palette(levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			custom_plot::patch::pie(draw_area, axis_rect, distribution, levels, colors, 0, 0, 1);
			custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, feature, gs);
			custom_plot::patch::set_range(axis_rect, QCPRange(-1.2, 1.2), QCPRange(-1.2, 1.2));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::patch::remove_left_bottom_axis(axis_rect);
			this->draw_suite_->update(draw_area);
		}
	}

	if (feature_data.type == QUERY_DATA::DataType::numeric) {

		auto [imin, imax] = std::ranges::minmax(feature_data.dd);
		double min = imin, max = imax;
		auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);

		if (gs.use_boxplot()) {
			custom_plot::patch::box(
				draw_area,
				axis_rect,
				colors[0],
				custom::cast<QVector>(feature_data.dd),
				1.0,
				gs.boxplot_draw_outlier(),
				gs.get_scatter_point_size()
			);
		}
		else {
			std::tie(min, max) = custom_plot::patch::violin(
				draw_area,
				axis_rect,
				colors[0],
				custom::cast<QVector>(feature_data.dd),
				1.0
			);
		}

		custom_plot::set_left_title(axis_rect, feature, gs, true);
		custom_plot::set_simple_axis_no_title(axis_rect, gs);
		custom_plot::patch::remove_bottom_ticks(axis_rect);
		custom_plot::patch::set_range(axis_rect, QCPRange(0, 2), custom_plot::utility::get_range(min, max));
		custom_plot::add_title(draw_area, feature, gs);
		this->draw_suite_->update(draw_area);
	}

	if (feature_data.type == QUERY_DATA::DataType::string) {

		if (!feature_data.dsl.isEmpty()) {

			auto distribution = custom::table(feature_data.ds);
			QStringList levels = custom::sorted(feature_data.dsl);

			auto colors = gs.palette(levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			custom_plot::patch::pie(draw_area, axis_rect, distribution, levels, colors, 0, 0, 1);
			custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, feature, gs);
			custom_plot::patch::set_range(axis_rect, QCPRange(-1.2, 1.2), QCPRange(-1.2, 1.2));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::patch::remove_left_bottom_axis(axis_rect);
			this->draw_suite_->update(draw_area);
		}
	}
};

void MetadataViewWindow::mainplot() {

	QString main_group_feature = this->main_group_box_->currentText();
	QString feature = this->feature_line_edit_->text();

	const auto& gs = this->draw_suite_->graph_settings_;

	bool normalize = this->normalize_switch_->value_;

	bool gene_activity{ false };
	if (this->gene_activity_switch_ != nullptr) {
		gene_activity = this->gene_activity_switch_->value_;
	}

	auto feature_data = this->handler_.get_data({ feature, normalize, gene_activity });
	auto group_data = this->handler_.get_data({ main_group_feature });

	if (feature_data.type == QUERY_DATA::DataType::notype || group_data.type == QUERY_DATA::DataType::notype) {
		G_WARN("No Data Found.");
		return;
	}

	QStringList group;
	QStringList levels;

	if (group_data.type == QUERY_DATA::DataType::integer) {
		if (group_data.dil.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		group = custom::cast<QString>(group_data.di);
		levels = custom::cast<QString>(custom::sorted(group_data.dil));
	}
	else if (group_data.type == QUERY_DATA::DataType::string) {
		if (group_data.dsl.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		group = group_data.ds;
		levels = custom::sorted(group_data.dsl);
	}
	else {
		G_WARN("Meeting error when querying data.");
		return;
	}

	int n_level = levels.size();
	auto colors = gs.palette(levels);

	if (feature_data.type == QUERY_DATA::DataType::integer) {
		if (feature_data.dsl.isEmpty()) {

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			double min{ 0.0 }, max{ 0.0 };

			if (gs.use_boxplot()) {
				auto [imin, imax] = std::ranges::minmax(feature_data.di);
				min = imin;
				max = imax;

				custom_plot::patch::box_batch(
					draw_area,
					axis_rect,
					group,
					levels,
					colors,
					custom::cast<Eigen::ArrayX>(custom::cast<double>(feature_data.di)),
					1.0,
					2.0,
					gs.boxplot_draw_outlier(),
					gs.get_scatter_point_size()
				);
			}
			else {

				std::tie(min, max) = custom_plot::patch::violin_batch(
					draw_area,
					axis_rect,
					group,
					levels,
					colors,
					custom::cast<Eigen::ArrayX>(custom::cast<double>(feature_data.di)),
					1.0,
					2.0);
			}

			custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, main_group_feature, gs);

			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * n_level), custom_plot::utility::get_range(min, max));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::set_simple_axis_no_title(axis_rect, gs);
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1),
				levels,
				6,
				gs
			);

			this->draw_suite_->update(draw_area);
		}
		else {

			QStringList feature_group = custom::cast<QString>(feature_data.di);
			QStringList feature_levels = custom::cast<QString>(custom::sorted(feature_data.dil));
			auto feature_colors = gs.palette(feature_levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			custom_plot::patch::bar_stack(draw_area, axis_rect, group, levels, feature_group, feature_levels, feature_colors, 1.0, 0.4, 1.0);
			custom_plot::add_round_legend(draw_area, legend_layout, feature_levels, feature_colors, feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_proportion_left_axis(axis_rect, gs.get_left_label_font());
			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_level + 1), QCPRange(0, 1));
			custom_plot::add_title(draw_area, feature, gs);

			custom_plot::patch::remove_bottom_axis(axis_rect);
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_level, 1, n_level),
				levels,
				6,
				gs
			);
			this->draw_suite_->update(draw_area);
		}
	}

	if (feature_data.type == QUERY_DATA::DataType::numeric) {

		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		double min{ 0.0 }, max{ 0.0 };

		if (gs.use_boxplot()) {
			auto [imin, imax] = std::ranges::minmax(feature_data.dd);
			min = imin;
			max = imax;

			custom_plot::patch::box_batch(
				draw_area,
				axis_rect,
				group,
				levels,
				colors,
				custom::cast<Eigen::ArrayX>(feature_data.dd),
				1.0,
				2.0,
				gs.boxplot_draw_outlier(),
				gs.get_scatter_point_size()
			);
		}
		else {

			std::tie(min, max) = custom_plot::patch::violin_batch(
				draw_area,
				axis_rect,
				group,
				levels,
				colors,
				custom::cast<Eigen::ArrayX>(feature_data.dd),
				1.0,
				2.0);
		}

		custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, main_group_feature, gs);

		custom_plot::set_left_title(axis_rect, feature, gs, true);
		custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * n_level), custom_plot::utility::get_range(min, max));
		custom_plot::add_title(draw_area, feature, gs);
		custom_plot::set_simple_axis_no_title(axis_rect, gs);
		custom_plot::set_bottom_axis_label(
			axis_rect,
			Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1),
			levels,
			6,
			gs
		);

		this->draw_suite_->update(draw_area);
	}

	if (feature_data.type == QUERY_DATA::DataType::string) {
		if (!feature_data.dsl.isEmpty()) {
			QStringList feature_group = feature_data.ds;
			QStringList feature_levels = custom::sorted(feature_data.dsl);

			auto feature_colors = gs.palette(feature_levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			custom_plot::patch::bar_stack(draw_area, axis_rect, group, levels, feature_group, feature_levels, feature_colors, 1.0, 0.4, 1.0);
			custom_plot::add_round_legend(draw_area, legend_layout, feature_levels, feature_colors, feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_proportion_left_axis(axis_rect, gs.get_left_label_font());
			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_level + 1), QCPRange(0, 1));
			custom_plot::add_title(draw_area, feature, gs);

			custom_plot::patch::remove_bottom_axis(axis_rect);
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_level, 1, n_level),
				levels,
				6,
				gs
			);
			this->draw_suite_->update(draw_area);
		}
	}
};

void MetadataViewWindow::subplot() {
	
	QString sub_group_feature = this->sub_group_box_->currentText();
	QString feature = this->feature_line_edit_->text();

	const auto& gs = this->draw_suite_->graph_settings_;

	bool normalize = this->normalize_switch_->value_;

	bool gene_activity{ false };
	if (this->gene_activity_switch_ != nullptr) {
		gene_activity = this->gene_activity_switch_->value_;
	}

	auto feature_data = this->handler_.get_data({ feature, normalize, gene_activity });
	auto group_data = this->handler_.get_data({ sub_group_feature });

	if (feature_data.type == QUERY_DATA::DataType::notype || group_data.type == QUERY_DATA::DataType::notype) {
		G_WARN("No Data Found.");
		return;
	}

	QStringList group;
	QStringList levels;

	if (group_data.type == QUERY_DATA::DataType::integer) {
		if (group_data.dil.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		group = custom::cast<QString>(group_data.di);
		levels = custom::cast<QString>(custom::sorted(group_data.dil));
	}
	else if (group_data.type == QUERY_DATA::DataType::string) {
		if (group_data.dsl.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		group = group_data.ds;
		levels = custom::sorted(group_data.dsl);
	}
	else {
		G_WARN("Meeting error when querying data.");
		return;
	}
	int n_level = levels.size();
	auto colors = gs.palette(levels);

	if (feature_data.type == QUERY_DATA::DataType::integer) {
		if (feature_data.dsl.isEmpty()) {

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

			custom_plot::set_simple_axis_no_title(axis_rect, gs);

			if (gs.use_boxplot()) {
				auto [imin, imax] = std::ranges::minmax(feature_data.di);
				double min = imin;
				double max = imax;

				custom_plot::patch::box_batch(
					draw_area,
					axis_rect,
					group,
					levels,
					colors,
					custom::cast<Eigen::ArrayX>(custom::cast<double>(feature_data.di)),
					1.0,
					2.0,
					gs.boxplot_draw_outlier(),
					gs.get_scatter_point_size()
				);
				custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * n_level), custom_plot::utility::get_range(min, max));
				custom_plot::set_bottom_axis_label(
					axis_rect,
					Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1),
					levels,
					6,
					gs
				);
			}
			else {

				if (n_level == 2 && gs.use_facet_violin_plot()) {
					auto [min, max] = custom_plot::patch::violin_facet(
						draw_area,
						axis_rect,
						group,
						levels,
						colors,
						custom::cast<Eigen::ArrayX>(custom::cast<double>(feature_data.di)),
						1.0);
					custom_plot::patch::set_range(axis_rect, QCPRange(0, 2), custom_plot::utility::get_range(min, max));
					custom_plot::set_bottom_title(axis_rect, sub_group_feature, gs, true);
				}
				else {
					auto [min, max] = custom_plot::patch::violin_batch(
						draw_area,
						axis_rect,
						group,
						levels,
						colors,
						custom::cast<Eigen::ArrayX>(custom::cast<double>(feature_data.di)),
						1.0,
						2.0);
					custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * n_level), custom_plot::utility::get_range(min, max));
					custom_plot::set_bottom_axis_label(
						axis_rect,
						Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1),
						levels,
						6,
						gs
					);
				}
			}

			custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, sub_group_feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::add_title(draw_area, feature, gs);

			this->draw_suite_->update(draw_area);
		}
		else {

			QStringList feature_group = custom::cast<QString>(feature_data.di);
			QStringList feature_levels = custom::cast<QString>(custom::sorted(feature_data.dil));

			auto feature_colors = gs.palette(feature_levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			custom_plot::patch::bar_stack(draw_area, axis_rect, group, levels, feature_group, feature_levels, feature_colors, 1.0, 0.4, 1.0);
			custom_plot::patch::set_proportion_left_axis(axis_rect, gs.get_left_label_font());
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_level, 1, n_level),
				levels,
				6,
				gs
			);
			custom_plot::add_round_legend(draw_area, legend_layout, feature_levels, feature_colors, feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_level + 1), QCPRange(0, 1));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::patch::remove_bottom_axis(axis_rect);
			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			this->draw_suite_->update(draw_area);
		}
	}

	if (feature_data.type == QUERY_DATA::DataType::numeric) {

		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::set_simple_axis_no_title(axis_rect, gs);

		if (gs.use_boxplot()) {
			auto [min, max] = std::ranges::minmax(feature_data.dd);

			custom_plot::patch::box_batch(
				draw_area,
				axis_rect,
				group,
				levels,
				colors,
				custom::cast<Eigen::ArrayX>(feature_data.dd),
				1.0,
				2.0,
				gs.boxplot_draw_outlier(),
				gs.get_scatter_point_size()
			);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * n_level), custom_plot::utility::get_range(min, max));
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1),
				levels,
				6,
				gs
			);
		}
		else {

			if (n_level == 2 && gs.use_facet_violin_plot()) {
				auto [min, max] = custom_plot::patch::violin_facet(
					draw_area,
					axis_rect,
					group,
					levels,
					colors,
					custom::cast<Eigen::ArrayX>(feature_data.dd),
					1.0);
				custom_plot::patch::set_range(axis_rect, QCPRange(0, 2), custom_plot::utility::get_range(min, max));
				custom_plot::set_bottom_title(axis_rect, sub_group_feature, gs, true);
			}
			else {

				auto [min, max] = custom_plot::patch::violin_batch(
					draw_area,
					axis_rect,
					group,
					levels,
					colors,
					custom::cast<Eigen::ArrayX>(feature_data.dd),
					1.0,
					2.0);
				custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * n_level), custom_plot::utility::get_range(min, max));
				custom_plot::set_bottom_axis_label(
					axis_rect,
					Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1),
					levels,
					6,
					gs
				);
			}
		}

		custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, sub_group_feature, gs);
		custom_plot::set_left_title(axis_rect, feature, gs, true);
		custom_plot::add_title(draw_area, feature, gs);

		this->draw_suite_->update(draw_area);
	}

	if (feature_data.type == QUERY_DATA::DataType::string) {
		if (!feature_data.dsl.isEmpty()) {
			QStringList feature_group = feature_data.ds;
			QStringList feature_levels = custom::sorted(feature_data.dsl);

			auto feature_colors = gs.palette(feature_levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
			custom_plot::patch::bar_stack(draw_area, axis_rect, group, levels, feature_group, feature_levels, feature_colors, 1.0, 0.4, 1.0);
			custom_plot::patch::set_proportion_left_axis(axis_rect, gs.get_left_label_font());
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_level, 1, n_level),
				levels,
				6,
				gs
			);
			custom_plot::add_round_legend(draw_area, legend_layout, feature_levels, feature_colors, feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_level + 1), QCPRange(0, 1));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::patch::remove_bottom_axis(axis_rect);
			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			this->draw_suite_->update(draw_area);
		}
	}
};

void MetadataViewWindow::plot() {

	QString main_group_feature = this->main_group_box_->currentText();
	QString sub_group_feature = this->sub_group_box_->currentText();
	QString feature = this->feature_line_edit_->text();

	const auto& gs = this->draw_suite_->graph_settings_;

	bool normalize = this->normalize_switch_->value_;

	bool gene_activity{ false };
	if (this->gene_activity_switch_ != nullptr) {
		gene_activity = this->gene_activity_switch_->value_;
	}

	auto feature_data = this->handler_.get_data({ feature, normalize, gene_activity });
	auto sub_group_data = this->handler_.get_data({ sub_group_feature });
	auto main_group_data = this->handler_.get_data({ main_group_feature });

	if (feature_data.type == QUERY_DATA::DataType::notype 
		|| main_group_data.type == QUERY_DATA::DataType::notype
		|| sub_group_data.type == QUERY_DATA::DataType::notype) {
		G_WARN("No Data Found.");
		return;
	}

	QStringList main_group;
	QStringList main_levels;

	if (main_group_data.type == QUERY_DATA::DataType::integer) {
		if (main_group_data.dil.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		main_group = custom::cast<QString>(main_group_data.di);
		main_levels = custom::cast<QString>(custom::sorted(main_group_data.dil));
	}
	else if (main_group_data.type == QUERY_DATA::DataType::string) {
		if (main_group_data.dsl.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		main_group = main_group_data.ds;
		main_levels = custom::sorted(main_group_data.dsl);
	}
	else {
		G_WARN("Meeting error when querying data.");
		return;
	}

	QStringList sub_group;
	QStringList sub_levels;

	if (sub_group_data.type == QUERY_DATA::DataType::integer) {
		if (sub_group_data.dil.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		sub_group = custom::cast<QString>(sub_group_data.di);
		sub_levels = custom::cast<QString>(custom::sorted(sub_group_data.dil));
	}
	else if (sub_group_data.type == QUERY_DATA::DataType::string) {
		if (sub_group_data.dsl.isEmpty()) {
			G_WARN("Meeting error when querying data.");
			return;
		}

		sub_group = sub_group_data.ds;
		sub_levels = custom::sorted(sub_group_data.dsl);
	}
	else {
		G_WARN("Meeting error when querying data.");
		return;
	}

	int n_main_level = main_levels.size();
	int n_sub_level = sub_levels.size();

	auto sub_colors = gs.palette(sub_levels);

	if (feature_data.type == QUERY_DATA::DataType::integer) {
		if (feature_data.dsl.isEmpty()) {

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

			custom_plot::set_simple_axis_no_title(axis_rect, gs);

			Eigen::ArrayXd data = custom::cast<Eigen::ArrayX>(custom::cast<double>(feature_data.di));

			auto [min, max] = std::ranges::minmax(data);

			for (int i = 0; i < n_main_level; ++i) {
				QString level = main_levels[i];
				auto filter = custom::equal(main_group, level);
				if (filter.count() == 0) {
					continue;
				}
				auto sub_feature_data = custom::sliced(data, filter);
				auto sub_sub_data = custom::sliced(sub_group, filter);

				if (gs.use_boxplot()) {

					custom_plot::patch::box_batch(
						draw_area,
						axis_rect,
						sub_sub_data,
						sub_levels,
						sub_colors,
						sub_feature_data,
						i* (n_sub_level * 2 + 1) + 1.0,
						2.0,
						gs.boxplot_draw_outlier(),
						gs.get_scatter_point_size()
					);
				}
				else {


					if (n_sub_level == 2 && gs.use_facet_violin_plot()) {
						auto [sub_min, sub_max] = custom_plot::patch::violin_facet(
							draw_area,
							axis_rect,
							sub_sub_data,
							sub_levels,
							sub_colors,
							sub_feature_data,
							2 * i + 1.0);

						if (sub_min < min) {
							min = sub_min;
						}

						if (sub_max > max) {
							max = sub_max;
						}
					}
					else {
						auto [sub_min, sub_max] = custom_plot::patch::violin_batch(
							draw_area,
							axis_rect,
							sub_sub_data,
							sub_levels,
							sub_colors,
							sub_feature_data,
							i * (n_sub_level * 2 + 1) + 1.0,
							2.0);

						if (sub_min < min) {
							min = sub_min;
						}

						if (sub_max > max) {
							max = sub_max;
						}
					}
				}
			}
			
			if (gs.use_boxplot()) {
				custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level* (n_sub_level * 2 + 1) - 1), custom_plot::utility::get_range(min, max));
				custom_plot::set_bottom_axis_label(
					axis_rect,
					Eigen::ArrayXd::LinSpaced(n_main_level, n_sub_level, n_main_level* (n_sub_level * 2 + 1) - 1 - n_sub_level),
					main_levels,
					6,
					gs
				);
			}
			else {
				if (n_sub_level == 2 && gs.use_facet_violin_plot()) {
					custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level * 2), custom_plot::utility::get_range(min, max));
					custom_plot::set_bottom_axis_label(
						axis_rect,
						Eigen::ArrayXd::LinSpaced(n_main_level, 1, 2 * n_main_level - 1),
						main_levels,
						6,
						gs
					);
				}
				else {
					custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level * (n_sub_level * 2 + 1) - 1), custom_plot::utility::get_range(min, max));
					custom_plot::set_bottom_axis_label(
						axis_rect,
						Eigen::ArrayXd::LinSpaced(n_main_level, n_sub_level, n_main_level * (n_sub_level * 2 + 1) - 1 - n_sub_level),
						main_levels,
						6,
						gs
					);
				}
			}

			custom_plot::add_round_legend(draw_area, legend_layout, sub_levels, sub_colors, sub_group_feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::add_title(draw_area, feature, gs);

			this->draw_suite_->update(draw_area);
		}
		else {

			QStringList data = custom::cast<QString>(feature_data.di);
			QStringList feature_levels = custom::cast<QString>(custom::sorted(feature_data.dil));

			auto feature_colors = gs.palette(feature_levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

			for (int i = 0; i < n_main_level; ++i) {
				QString level = main_levels[i];
				auto filter = custom::equal(main_group, level);
				if (filter.count() == 0) {
					continue;
				}
				auto sub_feature_data = custom::sliced(data, filter);
				auto sub_sub_data = custom::sliced(sub_group, filter);
				custom_plot::patch::bar_stack(
					draw_area, 
					axis_rect, 
					sub_sub_data, 
					sub_levels, 
					sub_feature_data, 
					feature_levels, 
					feature_colors,
					i * (n_sub_level + 1) + 1.0,
					0.4, 
					1.0
				);
			}

			custom_plot::add_round_legend(draw_area, legend_layout, feature_levels, feature_colors, feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level * (n_sub_level + 1)), QCPRange(0, 1));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::patch::set_proportion_left_axis(axis_rect, gs.get_left_label_font());
			custom_plot::patch::remove_bottom_axis(axis_rect);
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_main_level, (n_sub_level + 1) / 2., n_main_level * (n_sub_level + 1) - (n_sub_level + 1) / 2.),
				main_levels,
				6,
				gs
			);
			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			this->draw_suite_->update(draw_area);
		}
	}

	if (feature_data.type == QUERY_DATA::DataType::numeric) {

		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::set_simple_axis_no_title(axis_rect, gs);

		Eigen::ArrayXd data = custom::cast<Eigen::ArrayX>(feature_data.dd);

		auto [min, max] = std::ranges::minmax(data);
		
		auto colors = gs.palette(sub_levels);

		for (int i = 0; i < n_main_level; ++i) {
			QString level = main_levels[i];
			auto filter = custom::equal(main_group, level);
			if (filter.count() == 0) {
				continue;
			}
			auto sub_feature_data = custom::sliced(data, filter);
			auto sub_sub_data = custom::sliced(sub_group, filter);

			if (gs.use_boxplot()) {

				custom_plot::patch::box_batch(
					draw_area,
					axis_rect,
					sub_sub_data,
					sub_levels,
					sub_colors,
					sub_feature_data,
					i * (n_sub_level * 2 + 1) + 1.0,
					2.0,
					gs.boxplot_draw_outlier(),
					gs.get_scatter_point_size()
				);
			}
			else {

				if (n_sub_level == 2 && gs.use_facet_violin_plot()) {
					auto [sub_min, sub_max] = custom_plot::patch::violin_facet(
						draw_area,
						axis_rect,
						sub_sub_data,
						sub_levels,
						colors,
						sub_feature_data,
						2 * i + 1.0);

					if (sub_min < min) {
						min = sub_min;
					}

					if (sub_max > max) {
						max = sub_max;
					}
				}
				else {
					auto [sub_min, sub_max] = custom_plot::patch::violin_batch(
						draw_area,
						axis_rect,
						sub_sub_data,
						sub_levels,
						colors,
						sub_feature_data,
						i * (n_sub_level * 2 + 1) + 1.0,
						2.0);

					if (sub_min < min) {
						min = sub_min;
					}

					if (sub_max > max) {
						max = sub_max;
					}
				}
			}
		}

		if (gs.use_boxplot()) {

			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level* (n_sub_level * 2 + 1) - 1), custom_plot::utility::get_range(min, max));
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_main_level, n_sub_level, n_main_level* (n_sub_level * 2 + 1) - 1 - n_sub_level),
				main_levels,
				6,
				gs
			);
		}
		else {
			if (n_sub_level == 2 && gs.use_facet_violin_plot()) {
				custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level * 2), custom_plot::utility::get_range(min, max));
				custom_plot::set_bottom_axis_label(
					axis_rect,
					Eigen::ArrayXd::LinSpaced(n_main_level, 1, 2 * n_main_level - 1),
					main_levels,
					6,
					gs
				);
			}
			else {
				custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level * (n_sub_level * 2 + 1) - 1), custom_plot::utility::get_range(min, max));
				custom_plot::set_bottom_axis_label(
					axis_rect,
					Eigen::ArrayXd::LinSpaced(n_main_level, n_sub_level, n_main_level * (n_sub_level * 2 + 1) - 1 - n_sub_level),
					main_levels,
					6,
					gs
				);
			}
		}

		custom_plot::add_round_legend(draw_area, legend_layout, sub_levels, colors, sub_group_feature, gs);
		custom_plot::set_left_title(axis_rect, feature, gs, true);
		custom_plot::add_title(draw_area, feature, gs);

		this->draw_suite_->update(draw_area);
	}

	if (feature_data.type == QUERY_DATA::DataType::string) {
		if (!feature_data.dsl.isEmpty()) {

			QStringList data = feature_data.ds;
			QStringList feature_levels = custom::sorted(feature_data.dsl);

			int n_feature_level = feature_levels.size();
			auto feature_colors = gs.palette(feature_levels);

			auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

			for (int i = 0; i < n_main_level; ++i) {
				QString level = main_levels[i];
				auto filter = custom::equal(main_group, level);
				if (filter.count() == 0) {
					continue;
				}
				auto sub_feature_data = custom::sliced(data, filter);
				auto sub_sub_data = custom::sliced(sub_group, filter);
				custom_plot::patch::bar_stack(
					draw_area,
					axis_rect,
					sub_sub_data,
					sub_levels,
					sub_feature_data,
					feature_levels,
					feature_colors,
					i * (n_sub_level + 1) + 1.0,
					0.4,
					1.0
				);
			}

			custom_plot::add_round_legend(draw_area, legend_layout, feature_levels, feature_colors, feature, gs);
			custom_plot::set_left_title(axis_rect, feature, gs, true);
			custom_plot::patch::set_range(axis_rect, QCPRange(0, n_main_level * (n_sub_level + 1)), QCPRange(0, 1));
			custom_plot::add_title(draw_area, feature, gs);
			custom_plot::patch::set_proportion_left_axis(axis_rect, gs.get_left_label_font());
			custom_plot::patch::remove_bottom_axis(axis_rect);
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(n_main_level, (n_sub_level + 1) / 2., n_main_level * (n_sub_level + 1) - (n_sub_level + 1) / 2.),
				main_levels,
				6,
				gs
			);
			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			this->draw_suite_->update(draw_area);
		}
	}
}

void MetadataViewWindow::s_refresh_plot() {
	QString feature = this->feature_line_edit_->text();
	if (!this->valid_features_.contains(feature)) {
		G_LOG("Feature is not valid!");
		return;
	}

	QString main_group = this->main_group_box_->currentText();
	QString sub_group = this->sub_group_box_->currentText();

	if (main_group == "") {
		if (sub_group == "") {
			no_group_plot();
		}
		else {
			subplot();
		}
	}
	else {
		if (sub_group == "") {
			mainplot();
		}
		else {
			plot();
		}
	}
};
