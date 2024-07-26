#include "CoverageTrackWindow.h"

#include "CommonDialog.h"
#include "PlotWindow.h"
#include "GraphSettingDialog.h"
#include "GenomeUtility.h"
#include "CustomPlot.h"

#include "SoapGUI.h"

CoverageTrackWindow::CoverageTrackWindow(
	SingleCellMultiome* single_cell_multiome,
	CoverageTrack* coverage_track,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	coverage_track_(coverage_track),
	signal_emitter_(signal_emitter)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &CoverageTrackWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &CoverageTrackWindow::s_check_data);

	QStringList peak_names = single_cell_multiome->atac_counts()->rownames_;
	const qsizetype size = peak_names.size();
	for (qsizetype i = 0; i < size; ++i) {
		auto [sequence_name, start, end, success] = _Cs string_to_peak(peak_names[i]);
		if (success) {
			this->available_peak_locations_ << std::make_tuple(sequence_name, start, end);
		}
	}

	this->set_layout();

	this->set_property();
};

CoverageTrackWindow::CoverageTrackWindow(
	SingleCellAtac* single_cell_atac,
	CoverageTrack* coverage_track,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	coverage_track_(coverage_track),
	signal_emitter_(signal_emitter)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &CoverageTrackWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &CoverageTrackWindow::s_check_data);

	QStringList peak_names = single_cell_atac->counts()->rownames_;
	const qsizetype size = peak_names.size();
	for (qsizetype i = 0; i < size; ++i) {
		auto [sequence_name, start, end, success] = _Cs string_to_peak(peak_names[i]);
		if (success) {
			this->available_peak_locations_ << std::make_tuple(sequence_name, start, end);
		}
	}

	this->set_layout();

	this->set_property();
};

CoverageTrackWindow::CoverageTrackWindow(
	CoverageTrack* coverage_track,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	coverage_track_(coverage_track),
	signal_emitter_(signal_emitter)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &CoverageTrackWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &CoverageTrackWindow::s_check_data);

	this->set_layout();

	this->set_property();
};

void CoverageTrackWindow::view(
	SingleCellMultiome* single_cell_multiome,
	CoverageTrack* coverage_track,
	SignalEmitter* signal_emitter
) {
	CoverageTrackWindow* window = new CoverageTrackWindow(single_cell_multiome, coverage_track, signal_emitter);
};

void CoverageTrackWindow::view(
	SingleCellAtac* single_cell_atac,
	CoverageTrack* coverage_track,
	SignalEmitter* signal_emitter
) {
	CoverageTrackWindow* window = new CoverageTrackWindow(single_cell_atac, coverage_track, signal_emitter);
};

void CoverageTrackWindow::view(
	CoverageTrack* coverage_track,
	SignalEmitter* signal_emitter
) {
	CoverageTrackWindow* window = new CoverageTrackWindow(coverage_track, signal_emitter);
};

CoverageTrackWindow::~CoverageTrackWindow() {

	delete this->draw_suite_;
};

void CoverageTrackWindow::s_set_graph_settings() {
	GraphSettingDialog::set_graph_setting(this->draw_suite_->graph_settings_);
	this->graph_setting_switch_->set_status(this->draw_suite_->graph_settings_.active());
};

void CoverageTrackWindow::s_check_data(void* data, soap::VariableType type, void* item) {

	if (this->coverage_track_ == data) {
		this->close();
	}
};

void CoverageTrackWindow::set_layout() {

	this->main_interface_ = new QWidget(this);
	this->setCentralWidget(this->main_interface_);

	this->main_layout_ = new QVBoxLayout;
	this->main_interface_->setLayout(this->main_layout_);

	this->set_top_layout();
	this->set_middle_layout();
	this->set_bottom_layout();

	this->draw_suite_ = new PlotsSuite();

	connect(this->graph_setting_switch_, &Switch::toggled, this->draw_suite_, &PlotsSuite::s_setting_activate);
	connect(this->draw_suite_, &PlotsSuite::x_plot_prepared, this, &CoverageTrackWindow::s_new_plot);

	this->draw_suite_->prepare();
};

void CoverageTrackWindow::set_property() {
	this->setAttribute(Qt::WA_DeleteOnClose);
	G_SET_ICON;
	this->resize(800, 600);
	this->setWindowTitle("Coverage Track View");
	this->show();
};

void CoverageTrackWindow::set_top_layout() {
	G_SET_BUTTON(this->graph_setting_button_, "Graph Settings", soap::MiddleSize);
	connect(this->graph_setting_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_set_graph_settings);

	G_SET_SWITCH(this->graph_setting_switch_, false, this->graph_setting_button_, soap::MiddleSize);

	G_SET_PLOTSUITE_BUTTON;

	this->top_layout_ = new QHBoxLayout;

	this->top_layout_->addWidget(this->graph_setting_button_);
	this->top_layout_->addWidget(this->graph_setting_switch_);

	this->top_layout_->addStretch();

	this->top_layout_->addWidget(this->previous_picture_button_);
	this->top_layout_->addWidget(this->next_picture_button_);
	this->top_layout_->addWidget(this->pop_picture_button_);
	this->top_layout_->addWidget(this->clear_picture_button_);
	this->top_layout_->addWidget(this->save_picture_button_);

	this->save_picture_menu_ = new QMenu();

	this->save_png_action_ = new QAction("PNG", this);
	this->save_jpg_action_ = new QAction("JPG", this);
	this->save_bmp_action_ = new QAction("BMP", this);
	this->save_pdf_action_ = new QAction("PDF", this);
	this->save_pdf_and_png_action_ = new QAction("PDF and PNG", this);

	this->save_picture_menu_->addAction(this->save_png_action_);
	this->save_picture_menu_->addAction(this->save_jpg_action_);
	this->save_picture_menu_->addAction(this->save_bmp_action_);
	this->save_picture_menu_->addAction(this->save_pdf_action_);
	this->save_picture_menu_->addAction(this->save_pdf_and_png_action_);

	connect(this->save_jpg_action_, &QAction::triggered, this, &CoverageTrackWindow::s_save_jpg);
	connect(this->save_png_action_, &QAction::triggered, this, &CoverageTrackWindow::s_save_png);
	connect(this->save_bmp_action_, &QAction::triggered, this, &CoverageTrackWindow::s_save_bmp);
	connect(this->save_pdf_action_, &QAction::triggered, this, &CoverageTrackWindow::s_save_pdf);
	connect(this->save_pdf_and_png_action_, &QAction::triggered, this, &CoverageTrackWindow::s_save_pdf_and_png);

	this->save_picture_button_->setMenu(this->save_picture_menu_);

	this->main_layout_->addLayout(this->top_layout_);

	connect(this->previous_picture_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_previous_plot);
	connect(this->next_picture_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_next_plot);
	connect(this->pop_picture_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_pop_plot);
	connect(this->clear_picture_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_clear_plot);
};

void CoverageTrackWindow::set_bottom_layout() {
	this->bottom_layout_ = new QHBoxLayout;

	G_SET_BUTTON(this->left_button_, "◄", QSize(50, 30));
	G_SET_BUTTON(this->right_button_, "►", QSize(50, 30));

	connect(this->left_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_left);
	connect(this->right_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_right);

	G_SET_LABEL(this->region_label_, "Region", QSize(80, 30));
	G_SET_LINEEDIT(this->region_line_, "", QSize(200, 30));
	G_SET_BUTTON(this->go_button_, "go", QSize(100, 30));

	this->connect(this->go_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_go);

	this->bottom_layout_->addStretch();
	this->bottom_layout_->addWidget(this->left_button_);
	this->bottom_layout_->addWidget(this->right_button_);

	this->bottom_layout_->addStretch();

	this->bottom_layout_->addWidget(this->region_label_);
	this->bottom_layout_->addWidget(this->region_line_);
	this->bottom_layout_->addWidget(this->go_button_);

	this->bottom_layout_->addStretch();

	this->main_layout_->addLayout(this->bottom_layout_);
};

void CoverageTrackWindow::set_middle_layout() {
	this->middle_layout_ = new QHBoxLayout;

	this->set_figure_layout();
	this->set_information_layout();

	this->main_layout_->addLayout(this->middle_layout_);
}

void CoverageTrackWindow::set_information_layout() {
	this->information_layout_ = new QVBoxLayout;

	this->information_area_ = new InformationTextBrowser(this);
	this->information_area_->setFixedWidth(200);
	this->information_layout_->addWidget(this->information_area_);

	G_SET_LABEL(this->search_gene_location_label_, "Search Gene Location", soap::MiddleSize);
	G_SET_LINEEDIT(this->search_gene_location_line_, "", soap::MiddleSize);
	G_SET_BUTTON(this->search_gene_location_button_, "Search", soap::MiddleSize);
	G_SET_LABEL_FIXED_HEIGHT(this->gene_location_label_, "", 30);

	this->connect(this->search_gene_location_button_, &QPushButton::clicked, this, &CoverageTrackWindow::s_search_gene_location);

	this->information_layout_->addWidget(this->search_gene_location_label_);
	this->information_layout_->addWidget(this->search_gene_location_line_);
	this->information_layout_->addWidget(this->search_gene_location_button_);
	this->information_layout_->addWidget(this->gene_location_label_);

	this->middle_layout_->addLayout(this->information_layout_);
}

void CoverageTrackWindow::set_figure_layout() {
	this->figure_layout_ = new QVBoxLayout;
	this->middle_layout_->addLayout(this->figure_layout_);

};

void CoverageTrackWindow::s_save_png() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_, 
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

void CoverageTrackWindow::s_save_bmp() {

	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_, 
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

void CoverageTrackWindow::s_save_jpg() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_, 
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

void CoverageTrackWindow::s_save_pdf_and_png() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_, 
		"PDF-PNG Settings",
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
		G_LOG("PDF saving failed");
		return;
	}
	else {
		G_LOG("PDF saved.");
	}
	std::string pdf_file_name = pdf_name.toStdString();
	std::string picture_file_name = picture_name.toStdString();

	if (!_Cs save_pdf_page_as_png(pdf_file_name, 0, picture_file_name)) {
		G_WARN("PNG Saving Failed.");
	}
	else {
		G_LOG("PNG saved.");
	}
};

void CoverageTrackWindow::s_save_pdf() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_, 
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

void CoverageTrackWindow::s_previous_plot() {
	if (this->draw_suite_->current_plot_id_ <= 1)return;
	this->figure_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[-- this->draw_suite_->current_plot_id_];
	this->figure_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
};

void CoverageTrackWindow::s_next_plot() {
	if (this->draw_suite_->current_plot_id_ == this->draw_suite_->maximum_plot_id_)return;
	this->figure_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[++ this->draw_suite_->current_plot_id_];
	this->figure_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
};

void CoverageTrackWindow::s_clear_plot() {
	this->main_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);
	this->draw_area_ = nullptr;
	this->draw_suite_->clear();
	this->plot_location_.clear();
};

void CoverageTrackWindow::s_pop_plot() {
	PlotWindow::show_plot(this->draw_suite_, "Figure", this->draw_area_->width(), this->draw_area_->height(), this);
};

void CoverageTrackWindow::s_new_plot() {
	if (this->draw_area_ == nullptr) {
		this->draw_area_ = this->draw_suite_->current_plot_;
		this->figure_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();
	}
	else {
		this->figure_layout_->removeWidget(this->draw_area_);
		this->draw_area_->setVisible(false);
		this->draw_area_ = this->draw_suite_->current_plot_;
		this->figure_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();
	}

};

void CoverageTrackWindow::expand_region() {
	int extend_length = 0;
	if (this->get_location_by_gene_name_) {
		extend_length = (this->location_.end - this->location_.start) / 6;
		if (extend_length < 100) {
			extend_length = 100;
		}
	}
	else if (this->location_.end - this->location_.start < 200) {
		G_NOTICE("Given region is too narrow for visualization. Region is expanded.");
		extend_length = 100;
	}
	else {
		return;
	}
	this->location_.start -= extend_length;
	if (this->location_.start < 1) {
		this->location_.start = 1;
	}
	this->location_.end += extend_length;
};

bool CoverageTrackWindow::get_location() {
	auto [sequence_name, start, end, success] = _Cs string_to_peak(this->plot_elements_.region_name);
	if (!success) {
		return this->get_location_by_gene_name();
	}
	else {
		this->location_.sequence_name = sequence_name;
		this->location_.start = start;
		this->location_.end = end;
		this->get_location_by_gene_name_ = false;
		return true;
	}
};

bool CoverageTrackWindow::get_location_by_gene_name() {
	auto [seq_name, start, end, strand, success] = 
		_Cs find_gene_in_genome(this->plot_elements_.region_name, this->coverage_track_->annotation_);
	if (!success) {
		return false;
	}
	this->location_.sequence_name = seq_name;
	this->location_.start = start;
	this->location_.end = end;
	this->get_location_by_gene_name_ = true;
	return true;
};

bool CoverageTrackWindow::find_gene_in_region() {

	auto start = this->coverage_track_->annotation_.ranges_.start_;
	auto end = this->coverage_track_->annotation_.get_sequence_end();
	QVector<char> strand = this->coverage_track_->annotation_.strand_.to_qvector();
	auto sequence_filter = this->coverage_track_->annotation_.sequence_names_ == this->location_.sequence_name;
	if (sequence_filter.count() == 0) {
		G_NOTICE("No Gene Found in region : " + this->location_.sequence_name + ":" + QString::number(this->location_.start) + "-" + QString::number(this->location_.end) + ".");
		return false;
	}

	start = _Cs sliced(start, sequence_filter);
	end = _Cs sliced(end, sequence_filter);
	strand = _Cs sliced(strand, sequence_filter);
	QStringList type = _Cs sliced(this->coverage_track_->annotation_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_TYPE), sequence_filter);
	QStringList gene_name = _Cs sliced(this->coverage_track_->annotation_.metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME), sequence_filter);

	const qsizetype sequence_size = start.size();
	int last_start = 0, last_end = 0;
	for (qsizetype i = 0; i < sequence_size; ++i) {
		int local_start = start[i], local_end = end[i];
		if (local_start < this->location_.end && local_end > this->location_.start) {
			if (local_start == last_start && local_end == last_end) {
				continue;
			}
			if (type[i] != "gap") {
				this->plot_elements_.gene_structure[gene_name[i]].first = strand[i];
			}
			last_start = local_start;
			last_end = local_end;
			this->plot_elements_.gene_structure[gene_name[i]].second.append(std::make_tuple(local_start, local_end, type[i]));
		}
	}
	if (_Cs sum(_Cs sapply(this->plot_elements_.gene_structure.values(), [](auto&& value) {return value.second.size(); })) == 0) {
		G_NOTICE("No Gene Found in region : " + this->location_.sequence_name + ":" + QString::number(this->location_.start) + "-" + QString::number(this->location_.end) + ".");
		return false;
	}
	return true;
};

bool CoverageTrackWindow::find_peak_in_region() {
	if (this->available_peak_locations_.isEmpty()) {
		return false;
	}

	for (auto&& [sequence_name, start, end] : this->available_peak_locations_) {
		if (sequence_name == this->location_.sequence_name) {
			if (start < this->location_.end && end > this->location_.start) {
				this->plot_elements_.peak_locations << std::make_pair(start, end);
			}
		}
	}
	return this->plot_elements_.peak_locations.size() != 0;
};

void CoverageTrackWindow::smooth_matrix() {
	const Eigen::Index width = this->plot_elements_.normalized_matrix.cols();
	const Eigen::Index n_level = this->plot_elements_.normalized_matrix.rows();

	Eigen::MatrixXd average_matrix = this->plot_elements_.normalized_matrix / 10;

	for (Eigen::Index row = 0; row < n_level; ++row) {
		for (Eigen::Index col = 0; col < 6; ++col) {
			this->plot_elements_.normalized_matrix(row, col) = average_matrix.row(row).segment(0, col + 5).sum() / (col + 5) * 10;
		}
		double average = this->plot_elements_.normalized_matrix(row, 5);
		for (Eigen::Index col = 6; col < width - 5; ++col) {
			average += (average_matrix(row, col + 4) - average_matrix(row, col - 6));
			this->plot_elements_.normalized_matrix(row, col) = average;
		}
		for (Eigen::Index col = width - 4; col < width; ++col) {
			this->plot_elements_.normalized_matrix(row, col) = average_matrix.row(row).segment(col - 5, width - col + 5).sum() / (width - col + 5) * 10;
		}
	}

	double max_value = this->plot_elements_.normalized_matrix.maxCoeff();

	if (max_value == 0) {
		return;
	}

	this->plot_elements_.normalized_matrix /= max_value;
};

bool CoverageTrackWindow::calculate_matrix() {

	int start = (this->location_.start - 1) / 10;
	int end = (this->location_.end - 1) / 10;


	int width = end - start;
	int n_level = this->coverage_track_->levels_.size();

	this->plot_elements_.normalized_matrix = Eigen::MatrixXd::Zero(n_level, width);

	const auto& matrix = this->coverage_track_->insertion_matrix_[this->location_.sequence_name];

	for (int i = 0; i < n_level; ++i) {
		const auto& track = matrix[i];

		this->plot_elements_.normalized_matrix.row(i) = _Cs cast<Eigen::ArrayX>(track.segment(start, width)).cast<double>();
	}
	return this->plot_elements_.normalized_matrix.maxCoeff() != 0;
};

void CoverageTrackWindow::draw_plot() {

	auto draw_area = _Cp coverage_plot(this->plot_elements_, this->draw_suite_->graph_settings_);

	this->draw_suite_->update(draw_area);

	this->plot_elements_.clear();
}

void CoverageTrackWindow::s_search_gene_location() {
	this->plot_elements_.region_name = this->search_gene_location_line_->text();

	if (!this->get_location()) {
		G_WARN("Can not determine the gene location.");
		return;
	}

	this->gene_location_label_->setText(this->location_.sequence_name + ':'
		+ QString::number(this->location_.start) + "-" + QString::number(this->location_.end));
	this->gene_location_label_->adjustSize();
};

void CoverageTrackWindow::s_go() {

	this->plot_elements_.region_name = this->plot_elements_.plot_title = this->region_line_->text();

	if (!this->get_location()) {
		G_WARN("Can not determine the plot region.");
		return;
	}

	if (!this->coverage_track_->insertion_matrix_.contains(this->location_.sequence_name)) {
		G_WARN("Did not find sequence : " + this->location_.sequence_name + " in data.");
		return;
	}

	this->expand_region();

	this->update();
}

void CoverageTrackWindow::s_left() {
	if (!this->plot_location_.contains(this->draw_area_)) {
		return;
	}
	this->location_ = this->plot_location_[this->draw_area_];
	int width = this->location_.end - this->location_.start;

	this->location_.start -= width / 2;
	if (this->location_.start < 1) {
		this->location_.start = 1;
	}
	this->location_.end -= width / 2;

	this->plot_elements_.region_name = this->plot_elements_.plot_title =
		this->location_.sequence_name + ":" + QString::number(this->location_.start) +
		"-" + QString::number(this->location_.end);

	this->update();
};

void CoverageTrackWindow::s_right() {
	if (!this->plot_location_.contains(this->draw_area_)) {
		return;
	}
	this->location_ = this->plot_location_[this->draw_area_];
	int width = this->location_.end - this->location_.start;

	this->location_.start += width / 2;
	if (this->location_.start < 1) {
		this->location_.start = 1;
	}
	this->location_.end += width / 2;

	this->plot_elements_.region_name = this->plot_elements_.plot_title =
		this->location_.sequence_name + ":" + QString::number(this->location_.start) +
		"-" + QString::number(this->location_.end);

	this->update();
};

void CoverageTrackWindow::update() {

	this->plot_elements_.region = this->location_;
	this->plot_elements_.group_factors = this->coverage_track_->levels_;

	if (this->location_.end - this->location_.start > 1e7) {
		G_WARN("Region are two broad for visualization.");
		return;
	}

	this->find_gene_in_region();

	this->find_peak_in_region();

	if (!this->calculate_matrix()) {
		G_WARN("Did not find insertion in query region.");
		this->plot_elements_.clear();
		return;
	}

	this->smooth_matrix();

	this->draw_plot();

	this->plot_location_[this->draw_area_] = this->location_;
};