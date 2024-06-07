#include "EnrichmentItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "CustomPlot.h"
#include "ItemIOWorker.h"
#include "FileWritingWorker.h"

void EnrichmentItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " | " + QString::number(this->data()->mat_.cols()) + " ]");
};

void EnrichmentItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->mat_,
		this->title_,
		this->signal_emitter_);
}

void EnrichmentItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Show Significant", s_show_significant);

	ADD_MAIN_MENU("Visualize");

	ADD_MENU("Visualize | Barplot", "Barplot", "Visualize");
	ADD_ACTION("Default", "Visualize | Barplot", s_barplot_default);
	ADD_ACTION("Custom", "Visualize | Barplot", s_barplot);
	ADD_ACTION("Show Selected", "Visualize | Barplot", s_barplot_show_selected);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void EnrichmentItem::s_show_significant() {
	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Set P Threshold",
		{ "P Threhsold(<):0.05" },
		{ soap::InputStyle::StringLineEdit}
	);

	if (setting.isEmpty())return;

	double threshold = setting[0].toDouble();
	if (threshold <= 0 || threshold > 1) {
		G_WARN("Invalid p threshold.")
			return;
	}
	Eigen::ArrayX<bool> filter = _Cs less_than(this->data()->mat_.get_double(METADATA_ENRICHMENT_ADJUSTED_P_VALUE), threshold);
	if (filter.count() == 0) {
		G_WARN("No result meets requirement");
		return;
	}
	auto tmp = this->data()->mat_.row_sliced(filter);
	MatrixWindow::show_matrix(
		&tmp,
		"Significant Result",
		this->signal_emitter_);
};

void EnrichmentItem::s_barplot_default() {

	barplot();
}

void EnrichmentItem::s_barplot_show_selected() {
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Pathways",
		{ "Pathways" },
		{ soap::InputStyle::SimpleChoice},
		{ _Cs sorted(this->data()->mat_.get_const_qstring_reference(METADATA_ENRICHMENT_PATHWAY_NAMES))}
	);

	if (settings.isEmpty()) {
		return;
	}

	QStringList pathways = simple_choice_to_list(settings[0]);
	if (pathways.isEmpty()) {
		return;
	}
	barplot(pathways);
};

void EnrichmentItem::s_barplot() {
	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Barplot Setting",
		{ "Number of Pathway:10" },
		{ soap::InputStyle::IntegerLineEdit}
	);

	if (setting.isEmpty())return;

	int path_number = setting[0].toInt();
	if (path_number <= 0) {
		G_LOG("Number of Pathway should be positive!");
		return;
	}
	if (path_number < 1 || path_number > 100) {
		G_LOG("Invalid pathway number! Reset to 10.");
		path_number = 10;
	}
	barplot(path_number);
}

void EnrichmentItem::barplot(const QStringList& pathways) {
	const QStringList& all_pathways = this->data()->mat_.get_const_qstring_reference(METADATA_ENRICHMENT_PATHWAY_NAMES);
	auto index = _Cs valid_index_of(pathways, all_pathways);
	int path_number = index.size();
	if (path_number < 1) {
		return;
	}

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect, legend_layout] = _Cp prepare(gs);
	axis_rect->axis(QCPAxis::atLeft)->setBasePen(QPen(Qt::NoPen));
	axis_rect->axis(QCPAxis::atLeft)->setTickLength(0);
	axis_rect->axis(QCPAxis::atLeft)->setSubTickPen(QPen(Qt::NoPen));
	axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);

	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);

	QVector<double> x = _Cs linspaced(path_number, path_number, 1);
	QVector<double> y = _Cs reordered(this->data()->mat_.get_double(METADATA_ENRICHMENT_COUNT), index);
	QStringList labels = _Cs reordered(all_pathways, index);
	labels = _Cs sapply(labels, _Cs capitalize_first);
	Eigen::ArrayXd p_adjusted = _Cs cast<Eigen::ArrayX>(_Cs reordered(this->data()->mat_.get_double(METADATA_ENRICHMENT_ADJUSTED_P_VALUE), index));
	for (int i = 0; i < path_number; ++i) {
		if (p_adjusted[i] < 1e-100) {
			p_adjusted[i] = 1e-100;
		}
	}
	p_adjusted = -log10(p_adjusted);

	_Cp bar_plot_enrichment(
		draw_area,
		axis_rect,
		_Cs cast<Eigen::ArrayX>(x),
		_Cs cast<Eigen::ArrayX>(y),
		p_adjusted,
		24,
		"-log<sub>10</sub>(P<sub>adj</sub>)",
		"",
		gs,
		false
	);

	_Cp add_gradient_legend(
		draw_area,
		legend_layout,
		p_adjusted.minCoeff(),
		p_adjusted.maxCoeff(),
		"-log<sub>10</sub>(P<sub>adj</sub>)",
		gs
	);

	_Cp set_left_axis_label(axis_rect, Eigen::ArrayXd::LinSpaced(labels.size(), labels.size(), 1),
		labels, 0, gs);

	QPen pen(Qt::black);
	pen.setWidth(3);
	axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
	_Cp set_bottom_title(axis_rect, "Gene Number", gs, true);

	this->draw_suite_->update(draw_area);
};

void EnrichmentItem::barplot(int path_number) {
	if (this->data()->mat_.rows() == 0) {
		G_WARN("It's an empty enrichment.");
		return;
	}
	path_number = path_number > this->data()->mat_.rows() ? this->data()->mat_.rows() : path_number;

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect, legend_layout] = _Cp prepare(gs);

	_CpPatch remove_left_axis(axis_rect);

	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);

	QVector<double> x = _Cs linspaced(path_number, path_number, 1);
	QVector<double> y = this->data()->mat_.get_double(METADATA_ENRICHMENT_COUNT).sliced(0, path_number);
	auto labels = this->data()->mat_.get_const_qstring_reference(METADATA_ENRICHMENT_PATHWAY_NAMES).sliced(0, path_number);
	labels = _Cs sapply(labels, _Cs capitalize_first);
	Eigen::ArrayXd p_adjusted = _Cs cast<Eigen::ArrayX>(this->data()->mat_.get_double(METADATA_ENRICHMENT_ADJUSTED_P_VALUE).sliced(0, path_number));
	for (int i = 0; i < path_number; ++i) {
		if (p_adjusted[i] < 1e-100) {
			p_adjusted[i] = 1e-100;
		}
	}
	p_adjusted = -log10(p_adjusted);

	_Cp bar_plot_enrichment(
		draw_area,
		axis_rect,
		_Cs cast<Eigen::ArrayX>(x),
		_Cs cast<Eigen::ArrayX>(y),
		p_adjusted,
		24,
		"-log<sub>10</sub>(P<sub>adj</sub>)",
		"",
		gs,
		false
	);

	_Cp add_gradient_legend(
		draw_area,
		legend_layout,
		p_adjusted.minCoeff(),
		p_adjusted.maxCoeff(),
		"-log<sub>10</sub>(P<sub>adj</sub>)",
		gs
	);

	_Cp set_left_axis_label(axis_rect, Eigen::ArrayXd::LinSpaced(labels.size(), labels.size(), 1),
		labels, 0, gs);

	QPen pen(Qt::black);
	pen.setWidth(3);
	axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
	_Cp set_bottom_title(axis_rect, "Gene Number", gs, true);

	this->draw_suite_->update(draw_area);
};
