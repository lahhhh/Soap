#include "GSEAItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CustomPlot.h"
#include "FileWritingWorker.h"

void GSEAItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " ]");
}

void GSEAItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->mat_, 
		this->title_, 
		this->signal_emitter_, 
		false,
		this->data_);
}

void GSEAItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Mountain Plot", s_mountain_plot);

	ADD_MAIN_ACTION("Show Significant", s_show_significant);

	ADD_MAIN_ACTION("Show Selected Pathway", s_show_selected);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void GSEAItem::s_show_significant() {
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Significance Setting", 
		{ "P Threshold(<):0.05", "P Type", "Regulation Type" },
		{ soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox, soap::InputStyle::ComboBox},
		{ { "Nominal", "FDR", "FWER" }, { "All", "Upregulated", "Downregulated" } }
	);

	if (settings.isEmpty())return;

	double threshold = settings[0].toDouble();
	if (threshold <= 0 || threshold > 1) {
		G_LOG("Invalid p threshold. Reset to 0.05");
		threshold = 0.05;
	}

	Eigen::ArrayX<bool> filter;
	if (settings[1] == "Nominal") {
		filter = _Cs less_than(this->data()->mat_.get_const_double_reference(METADATA_GSEA_P_VALUE), threshold);
	}
	else if (settings[1] == "FDR") {
		filter = _Cs less_than(this->data()->mat_.get_const_double_reference(METADATA_GSEA_FALSE_DISCOVERY_RATE), threshold);
	}
	else {
		filter = _Cs less_than(this->data()->mat_.get_const_double_reference(METADATA_GSEA_FAMILY_WISE_ERROR_RATE), threshold);
	}

	QString regulate_type = settings[2];
	if (regulate_type == "Upregulated") {
		filter *= _Cs greater_than(this->data()->mat_.get_const_double_reference(METADATA_GSEA_ENRICHMENT_SCORE), 0.);
	}
	else if (regulate_type == "Downregulated") {
		filter *= _Cs less_than(this->data()->mat_.get_const_double_reference(METADATA_GSEA_ENRICHMENT_SCORE), 0.);
	}

	if (filter.count() == 0) {
		G_NOTICE("No Result is Significant");
		return;
	}

	auto tmp = this->data()->mat_.row_sliced(filter);
	MatrixWindow::show_matrix(
		&tmp, 
		"Significant Result", 
		this->signal_emitter_);
}

void GSEAItem::single_mountain_plot(const QString& path_name, const QString& custom_name) {

	if (!this->data()->point_x_.contains(path_name)) {
		return;
	}

	QString label = custom_name.isEmpty() ? path_name : custom_name;

	int gene_number = this->data()->correlations_.size();

	const auto& gs = this->draw_suite_->graph_settings_;
	QCustomPlot* draw_area = _Cp initialize_plot(gs);
	draw_area->plotLayout()->setRowSpacing(-10);
		
	draw_area->setBackground(QColor(242, 242, 242));

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	QPen pen;
	pen.setColor(Qt::green);
	pen.setWidth(4);

	const auto& x = this->data()->point_x_[path_name];
	const auto& y = this->data()->point_y_[path_name];

	QCPAxisRect* top_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(0, 0, top_rect);
	draw_area->addGraph(top_rect->axis(QCPAxis::atBottom), top_rect->axis(QCPAxis::atLeft));
	draw_area->graph()->setName("Enrichment Profile");
	draw_area->graph()->setPen(pen);
	draw_area->graph()->setData(x, y);

	auto [ymin, ymax] = std::ranges::minmax(y);
	_CpPatch set_range(top_rect, QCPRange(-1, gene_number - 1), QCPRange(ymin, ymax));

	top_rect->axis(QCPAxis::atLeft)->ticker()->setTickCount(6);
	top_rect->axis(QCPAxis::atBottom)->ticker()->setTickCount(11);
	top_rect->axis(QCPAxis::atLeft)->setSubTicks(false);
	top_rect->axis(QCPAxis::atLeft)->setTickLengthIn(0);
	top_rect->axis(QCPAxis::atLeft)->setTickLengthOut(5);
	top_rect->axis(QCPAxis::atBottom)->setTicks(false);
	top_rect->axis(QCPAxis::atTop)->setVisible(true);
	top_rect->axis(QCPAxis::atTop)->setTicks(false);
	top_rect->axis(QCPAxis::atRight)->setVisible(true);
	top_rect->axis(QCPAxis::atRight)->setTicks(false);
	top_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	top_rect->axis(QCPAxis::atLeft)->setLabel("Enrichment score (ES)");
	top_rect->axis(QCPAxis::atLeft)->setLabelFont(gs.get_left_title_font());
	top_rect->setBackground(Qt::white);
	top_rect->setMargins({ 0, 0, 0, 0 });
	top_rect->setMinimumMargins({ 0, 0, 0, 0 });

	QCPAxisRect* middle_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(1, 0, middle_rect);
	draw_area->addGraph(middle_rect->axis(QCPAxis::atBottom), middle_rect->axis(QCPAxis::atLeft));
	draw_area->graph()->setName("Hits");
	draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
	draw_area->graph()->setPen(QPen(Qt::black));
	draw_area->graph()->setData(_Cs cast<double>(this->data()->gene_location_[path_name]), QVector<double>(this->data()->gene_location_[path_name].size(), 1));
	middle_rect->axis(QCPAxis::atTop)->setVisible(true);
	middle_rect->axis(QCPAxis::atTop)->setTicks(false);
	middle_rect->axis(QCPAxis::atRight)->setVisible(true);
	middle_rect->axis(QCPAxis::atRight)->setTicks(false);
	middle_rect->axis(QCPAxis::atBottom)->setVisible(true);
	middle_rect->axis(QCPAxis::atBottom)->setTicks(false);
	middle_rect->axis(QCPAxis::atLeft)->setVisible(true);
	middle_rect->axis(QCPAxis::atLeft)->setTicks(false);
	middle_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
	middle_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	ymin = 0;
	ymax = 1;
	_CpPatch set_range(middle_rect, QCPRange(-1, gene_number - 1), QCPRange(ymin, ymax));
	middle_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	middle_rect->setBackground(Qt::white);
	middle_rect->setMargins({ 0, 0, 0, 0 });
	middle_rect->setMinimumMargins({ 0, 0, 0, 0 });

	QCPAxisRect* bottom_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(2, 0, bottom_rect);
	draw_area->addGraph(bottom_rect->axis(QCPAxis::atBottom), bottom_rect->axis(QCPAxis::atLeft));
	draw_area->graph()->setName("Signal to Noise");
	draw_area->graph()->setData(_Cs linspaced(gene_number, 0, gene_number - 1), this->data()->correlations_);
	draw_area->graph()->setPen(Qt::NoPen);
	draw_area->graph()->setBrush(QBrush(Qt::gray));

	ymin = std::ranges::min(this->data()->correlations_);
	ymax = std::ranges::max(this->data()->correlations_);

	bottom_rect->axis(QCPAxis::atTop)->setVisible(true);
	bottom_rect->axis(QCPAxis::atTop)->setTicks(false);
	bottom_rect->axis(QCPAxis::atRight)->setVisible(true);
	bottom_rect->axis(QCPAxis::atRight)->setTicks(false);
	bottom_rect->axis(QCPAxis::atBottom)->setRange(-1, gene_number - 1);
	bottom_rect->axis(QCPAxis::atLeft)->setRange(ymin, ymax);
	bottom_rect->axis(QCPAxis::atBottom)->setSubTicks(false);
	bottom_rect->axis(QCPAxis::atLeft)->setSubTicks(false);
	bottom_rect->axis(QCPAxis::atLeft)->setTickPen(Qt::NoPen);
	bottom_rect->axis(QCPAxis::atBottom)->setTickLengthIn(0);
	bottom_rect->axis(QCPAxis::atBottom)->setTickLengthOut(5);
	bottom_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	bottom_rect->axis(QCPAxis::atLeft)->setLabel("Signal to Noise");
	bottom_rect->axis(QCPAxis::atLeft)->setLabelFont(gs.get_left_title_font());
	bottom_rect->axis(QCPAxis::atBottom)->setLabel("Rank in Ordered Dataset");
	bottom_rect->axis(QCPAxis::atBottom)->setLabelFont(gs.get_left_title_font());
	bottom_rect->axis(QCPAxis::atLeft)->ticker()->setTickCount(6);
	bottom_rect->axis(QCPAxis::atBottom)->ticker()->setTickCount(11);
	bottom_rect->setBackground(Qt::white);
	bottom_rect->setMargins({ 0, 0, 0, 0 });
	bottom_rect->setMinimumMargins({ 0, 0, 0, 0 });

	QCPItemText* label_tracer_text = new QCPItemText(draw_area);
	label_tracer_text->setClipToAxisRect(false);
	label_tracer_text->position->setAxisRect(bottom_rect);
	label_tracer_text->position->setAxes(bottom_rect->axis(QCPAxis::atBottom), bottom_rect->axis(QCPAxis::atLeft));
	label_tracer_text->position->setType(QCPItemPosition::ptAxisRectRatio);
	label_tracer_text->setPositionAlignment(Qt::AlignTop | Qt::AlignLeft);
	label_tracer_text->position->setCoords(0, 0);
	label_tracer_text->setText(this->data()->comparison_[1] + " (positively correlated)");
	label_tracer_text->setFont(QFont("Arial", 10));
	label_tracer_text->setColor(Qt::red);

	QCPItemText* label_tracer_text2 = new QCPItemText(draw_area);
	label_tracer_text2->setClipToAxisRect(false);
	label_tracer_text2->position->setAxisRect(bottom_rect);
	label_tracer_text2->position->setAxes(bottom_rect->axis(QCPAxis::atBottom), bottom_rect->axis(QCPAxis::atLeft));
	label_tracer_text2->position->setType(QCPItemPosition::ptAxisRectRatio);
	label_tracer_text2->setPositionAlignment(Qt::AlignBottom | Qt::AlignLeft);
	label_tracer_text2->position->setCoords(0.6, 1);
	label_tracer_text2->setText(this->data()->comparison_[2] + " (negatively correlated)");
	label_tracer_text2->setFont(QFont("Arial", 10));
	label_tracer_text2->setColor(Qt::blue);

	int length = this->data()->bar_colors_.size();
	for (int i = 0; i < length; ++i) {
		draw_area->addGraph(middle_rect->axis(QCPAxis::atBottom), middle_rect->axis(QCPAxis::atLeft));
		draw_area->graph()->setData(QVector<double>() << this->data()->bar_points_[i] << this->data()->bar_points_[i + 1], QVector<double>() << 0.2 << 0.2);
		draw_area->graph()->setBrush(this->data()->bar_colors_[i]);
		draw_area->graph()->setPen(Qt::NoPen);
	}

	QCPAxisRect* legend_rect = new QCPAxisRect(draw_area, false);

	QCPLegend* legend = new QCPLegend;
	legend_rect->setMargins(QMargins(0, 35, 0, 0));
	legend_rect->insetLayout()->addElement(legend, Qt::AlignBottom | Qt::AlignCenter);
	legend->setVisible(true);
	legend->setLayer("legend");
	draw_area->setAutoAddPlottableToLegend(false);
	legend->addElement(0, 0, new QCPPlottableLegendItem(legend, draw_area->graph(0)));
	legend->addElement(0, 1, new QCPPlottableLegendItem(legend, draw_area->graph(1)));
	legend->addElement(0, 2, new QCPPlottableLegendItem(legend, draw_area->graph(2)));
	legend->setFont(QFont("Arial", 16));

	draw_area->plotLayout()->insertRow(0);

	SoapTextElement* plot_title = new SoapTextElement(
		draw_area, 
		gs.get_title(label), 
		gs.get_title_font()
	);

	plot_title->setMargins(QMargins(0, 0, 0, 30));
	draw_area->plotLayout()->addElement(0, 0, plot_title);
	draw_area->plotLayout()->addElement(4, 0, legend_rect);
	draw_area->plotLayout()->setRowStretchFactor(0, 1);
	draw_area->plotLayout()->setRowStretchFactor(1, 4);
	draw_area->plotLayout()->setRowStretchFactor(2, 2);
	draw_area->plotLayout()->setRowStretchFactor(3, 3);
	draw_area->plotLayout()->setRowStretchFactor(4, 1);
	this->draw_suite_->update(draw_area);
};

void GSEAItem::mountain_plot_patch(
	QCustomPlot* draw_area,
	QCPLayoutGrid* layout,
	const QString& path_name,
	const QString& path_label
) {
	QString label = path_label.isEmpty() ? path_name : path_label;

	int gene_number = this->data()->correlations_.size();
	const auto& gs = this->draw_suite_->graph_settings_;

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	QPen pen;
	pen.setColor(Qt::green);
	pen.setWidth(4);

	const auto& x = this->data()->point_x_[path_name];
	const auto& y = this->data()->point_y_[path_name];

	QCPAxisRect* top_rect = new QCPAxisRect(draw_area, true);
	layout->addElement(0, 0, top_rect);
	draw_area->addGraph(top_rect->axis(QCPAxis::atBottom), top_rect->axis(QCPAxis::atLeft));
	draw_area->graph()->setPen(pen);
	draw_area->graph()->setData(x, y);

	auto [ymin, ymax] = std::ranges::minmax(y);
	_CpPatch set_range(top_rect, QCPRange(-1, gene_number - 1), QCPRange(ymin, ymax));

	top_rect->axis(QCPAxis::atLeft)->ticker()->setTickCount(6);
	top_rect->axis(QCPAxis::atBottom)->ticker()->setTickCount(11);
	top_rect->axis(QCPAxis::atLeft)->setSubTicks(false);
	top_rect->axis(QCPAxis::atLeft)->setTickLengthIn(0);
	top_rect->axis(QCPAxis::atLeft)->setTickLengthOut(5);
	top_rect->axis(QCPAxis::atBottom)->setTicks(false);
	top_rect->axis(QCPAxis::atTop)->setVisible(true);
	top_rect->axis(QCPAxis::atTop)->setTicks(false);
	top_rect->axis(QCPAxis::atRight)->setVisible(true);
	top_rect->axis(QCPAxis::atRight)->setTicks(false);
	top_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	top_rect->axis(QCPAxis::atLeft)->setLabel("Enrichment score (ES)");
	top_rect->axis(QCPAxis::atLeft)->setLabelFont(gs.get_left_title_font());
	top_rect->setMargins(QMargins(0, 0, 0, 0));

	QCPAxisRect* middle_rect = new QCPAxisRect(draw_area, true);
	layout->addElement(1, 0, middle_rect);
	draw_area->addGraph(middle_rect->axis(QCPAxis::atBottom), middle_rect->axis(QCPAxis::atLeft));
	draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
	draw_area->graph()->setPen(QPen(Qt::black));
	draw_area->graph()->setData(_Cs cast<double>(this->data()->gene_location_[path_name]), QVector<double>(this->data()->gene_location_[path_name].size(), 1));
	middle_rect->axis(QCPAxis::atTop)->setVisible(true);
	middle_rect->axis(QCPAxis::atTop)->setTicks(false);
	middle_rect->axis(QCPAxis::atRight)->setVisible(true);
	middle_rect->axis(QCPAxis::atRight)->setTicks(false);
	middle_rect->axis(QCPAxis::atBottom)->setVisible(true);
	middle_rect->axis(QCPAxis::atBottom)->setTicks(false);
	middle_rect->axis(QCPAxis::atLeft)->setVisible(true);
	middle_rect->axis(QCPAxis::atLeft)->setTicks(false);
	middle_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
	middle_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	ymin = 0;
	ymax = 1;
	_CpPatch set_range(middle_rect, QCPRange(-1, gene_number - 1), QCPRange(ymin, ymax));
	middle_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	middle_rect->setMargins(QMargins(0, 0, 0, 0));

	int length = this->data()->bar_colors_.size();
	for (int i = 0; i < length; ++i) {
		draw_area->addGraph(middle_rect->axis(QCPAxis::atBottom), middle_rect->axis(QCPAxis::atLeft));
		draw_area->graph()->setData(QVector<double>() << this->data()->bar_points_[i] << this->data()->bar_points_[i + 1], QVector<double>() << 0.2 << 0.2);
		draw_area->graph()->setBrush(this->data()->bar_colors_[i]);
		draw_area->graph()->setPen(Qt::NoPen);
	}

	layout->insertRow(0);

	SoapTextElement* plot_title = new SoapTextElement(
		draw_area,
		gs.get_title(label),
		gs.get_title_font()
	);

	plot_title->setMargins(QMargins(0, 0, 0, 30));
	layout->addElement(0, 0, plot_title);
	layout->setRowStretchFactor(0, 1);
	layout->setRowStretchFactor(1, 4);
	layout->setRowStretchFactor(2, 2);
	layout->setRowSpacing(-30);
}

void GSEAItem::s_mountain_plot() {

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Mountain Plot Setting", 
		{ "Pathway:Custom Pathway Name", "Number of row:1" },
		{ soap::InputStyle::MultipleDoubleLineEditWithCompleterLayout, soap::InputStyle::IntegerLineEdit},
		{ this->data()->mat_.rownames_}
	);

	if (settings.isEmpty())return;

	auto pathways = multiple_double_line_edit_with_completer_layout_to_pair(settings[0]);

	if (pathways.first.isEmpty()) {
		return;
	}

	if (pathways.first.size() == 1) {
		this->single_mountain_plot(pathways.first[0], pathways.second[0]);
		return;
	}

	int nrow = settings[1].toInt();
	if (nrow < 1) {
		G_WARN("Number of rows can not be less than 1.");
		return;
	}


	QStringList path_names = pathways.first;
	QStringList labels = pathways.second;

	int n_path = path_names.size();

	QStringList valid_path;
	QStringList valid_labels;

	for (int i = 0; i < n_path; ++i) {
		if (this->data()->point_x_.contains(path_names[i])) {
			valid_path << path_names[i];
			valid_labels << labels[i];
		}
		else {
			G_WARN(path_names[i] + " is not found in GSEA data.");
		}
	}

	if (valid_path.isEmpty()) {
		G_WARN("No Valid Pathway Found.");
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, layout] = _Cp prepare_lg(gs);

	const int n_valid_path = valid_path.size();

	for (int i = 0; i < n_valid_path; ++i) {
		QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
		int col = i / nrow;
		int row = i - col * nrow;
		layout->addElement(row, col, sub_layout);
		this->mountain_plot_patch(
			draw_area,
			sub_layout,
			valid_path[i],
			valid_labels[i]
		);
	}

	this->draw_suite_->update(draw_area);

}

void GSEAItem::s_show_selected() {

	QStringList pathways = this->data()->mat_.rownames_;

	if (pathways.isEmpty()) {
		G_WARN("Empty Objecy");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Pathway",
		{ "Choose Pathway" },
		{soap::InputStyle::LineEditWithCompleter2},
		{pathways}
	);

	if (settings.isEmpty()) {
		return;
	}

	auto filter = _Cs valid_index_of(QStringList{ settings[0] }, pathways);

	if (filter.isEmpty()) {
		return;
	}

	auto tmp = this->data()->mat_.row_reordered(filter);
	MatrixWindow::show_matrix(
		&tmp,
		settings[0],
		this->signal_emitter_);
};