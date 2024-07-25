#include "CoverageTrackItem.h"

#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CoverageTrackWindow.h"
#include "GenomeUtility.h"
#include "CustomPlot.h"

void CoverageTrackItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->levels_.size()) + " ]");
};

void CoverageTrackItem::__set_menu() {
	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Show Coverage", s_show_coverage);

	ADD_MAIN_MENU("Export");

	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void CoverageTrackItem::__show_this() {

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

		CoverageTrackWindow::view(single_cell_multiome, this->data(), this->signal_emitter_);
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);

		CoverageTrackWindow::view(single_cell_atac, this->data(), this->signal_emitter_);
	}
	else {

		CoverageTrackWindow::view(this->data(), this->signal_emitter_);
	}
}

void CoverageTrackItem::s_show_coverage() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This data is not part of single-cell multiome data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Coverage Plot Settings",
		QStringList{ "Region" },
		QList<soap::InputStyle>{  soap::InputStyle::MultipleLineEdit}
	);
	if (settings.isEmpty()) {
		return;
	}

	auto regions = multiple_line_edit_to_list(settings[0]);

	if (regions.isEmpty()) {
		return;
	}
	QList<Location> locs;
	QList<bool> is_gene_name;
	QStringList valid_regions;
	for (auto&& region : regions) {
		auto [seq_name, start, end, success] = _Cs string_to_peak(region);
		if (success) {
			if (end - start > 1e6) {
				G_WARN("Region: " + region + " is too broad for visualization.");
				continue;
			}
			if (this->data()->insertion_matrix_.contains(seq_name)) {
				locs << Location{ seq_name, start, end };
				is_gene_name << false;
				valid_regions << region;
			}
		}
		else {
			auto [seq_name, start, end, strand, success] = _Cs find_gene_in_genome(region, this->data()->annotation_);
			if (success) {
				if (end - start > 1e6) {
					G_WARN("Region: " + region + " is too broad for visualization.");
					continue;
				}
				if (this->data()->insertion_matrix_.contains(seq_name)) {
					locs << Location{ seq_name, start, end };
					is_gene_name << true;
					valid_regions << region;
				}
			}
		}
	}

	if (locs.isEmpty()) {
		G_WARN("No Valid Region.");
		return;
	}

	int n_region = locs.size();

	for (int i = 0; i < n_region; ++i) {
		int extend_length = 0;
		if (is_gene_name[i]) {
			extend_length = (locs[i].end - locs[i].start) / 6;
			if (extend_length < 100) {
				extend_length = 100;
			}
		}
		else if (locs[i].end - locs[i].start < 200) {
			extend_length = 100;
		}
		else {
			continue;
		}
		locs[i].start -= extend_length;
		if (locs[i].start < 1) {
			locs[i].start = 1;
		}
		locs[i].end += extend_length;
	}

	auto& gs = this->draw_suite_->graph_settings_;
	auto draw_area = _Cp initialize_plot(gs);

	draw_area->setNoAntialiasingOnDrag(true);

	QCPLayoutGrid* main_layout = new QCPLayoutGrid;
	draw_area->plotLayout()->addElement(0, 0, main_layout);

	QCPLayoutGrid* legend_layout = new QCPLayoutGrid;
	QCPAxisRect* legend_top = new QCPAxisRect(draw_area, false);
	QCPAxisRect* legend_bottom = new QCPAxisRect(draw_area, false);

	QCPLayoutGrid* right_layout = new QCPLayoutGrid;
	draw_area->plotLayout()->addElement(0, 1, right_layout);

	right_layout->addElement(0, 0, legend_top);
	right_layout->addElement(1, 0, legend_layout);
	right_layout->addElement(2, 0, legend_bottom);

	legend_layout->setRowSpacing(15);

	right_layout->setMargins(QMargins(0, 0, 20, 0));

	const int n_level = this->data()->levels_.size();

	QCPLayoutGrid* title_layout = new QCPLayoutGrid;
	main_layout->addElement(0, 1, title_layout);

	for (int j = 0; j < n_region; ++j) {
		SoapTextElement* title = new SoapTextElement(draw_area, valid_regions[j], gs.get_title_font());
		title_layout->addElement(0, j, title);
	}

	QPen annotation_axis_pen(Qt::black), coverage_bottom_pen(Qt::black);
	annotation_axis_pen.setWidth(3);
	coverage_bottom_pen.setWidth(1);

	QCPMarginGroup* annotation_column_margin = new QCPMarginGroup(draw_area);
	QList<QCPMarginGroup*> column_margins(n_region);
	for (int i = 0; i < n_region; ++i) {
		column_margins[i] = new QCPMarginGroup(draw_area);
	}

	QCPAxisRect* coverage_annotation_rect = new QCPAxisRect(draw_area);
	_CpPatch remove_bottom_axis(coverage_annotation_rect);
	_CpPatch clear_left_axis(coverage_annotation_rect);
	coverage_annotation_rect->setMarginGroup(QCP::msLeft | QCP::msRight, annotation_column_margin);
	main_layout->addElement(1, 0, coverage_annotation_rect);

	coverage_annotation_rect->axis(QCPAxis::atLeft)->setBasePen(annotation_axis_pen);
	_Cp set_left_title(coverage_annotation_rect, "Coverage", gs);

	QCPLayoutGrid* coverage_layout = new QCPLayoutGrid;
	coverage_layout->setRowSpacing(0);
	main_layout->addElement(1, 1, coverage_layout);

	main_layout->setColumnStretchFactors(QList<double>{ 1, 5.0 * n_level});

	main_layout->setRowStretchFactor(0, std::round(0.7 * n_level));

	auto colors = gs.palette(this->data()->levels_);

	for (int j = 0; j < n_region; ++j) {

		int start = (locs[j].start - 1) / 10;
		int end = (locs[j].end - 1) / 10;
		int width = end - start;

		QVector<double> bottom_axis_value = _Cs linspaced(width, locs[j].start, locs[j].start + (width - 1) * 10);

		Eigen::MatrixXd insertion_matrix = Eigen::MatrixXd::Zero(n_level, width);

		const auto& matrix = this->data()->insertion_matrix_[locs[j].sequence_name];

		for (int i = 0; i < n_level; ++i) {
			const auto& track = matrix[i];

			insertion_matrix.row(i) = _Cs cast<Eigen::ArrayX>(track.segment(start, width)).cast<double>();
		}

		Eigen::MatrixXd average_matrix = insertion_matrix / 10;

		for (Eigen::Index row = 0; row < n_level; ++row) {
			for (Eigen::Index col = 0; col < 6; ++col) {
				insertion_matrix(row, col) = average_matrix.row(row).segment(0, col + 5).sum() / (col + 5) * 10;
			}
			double average = insertion_matrix(row, 5);
			for (Eigen::Index col = 6; col < width - 5; ++col) {
				average += (average_matrix(row, col + 4) - average_matrix(row, col - 6));
				insertion_matrix(row, col) = average;
			}
			for (Eigen::Index col = width - 4; col < width; ++col) {
				insertion_matrix(row, col) = average_matrix.row(row).segment(col - 5, width - col + 5).sum() / (width - col + 5) * 10;
			}
		}

		double max_value = insertion_matrix.maxCoeff();

		if (max_value != 0) {
			insertion_matrix /= max_value;
		}

		double threshold = insertion_matrix.maxCoeff() / 100;

		for (int i = 0; i < n_level; ++i) {
			QCPAxisRect* group_rect = new QCPAxisRect(draw_area);
			coverage_layout->addElement(i, j, group_rect);

			_CpPatch remove_left_axis(group_rect);
			_CpPatch clear_bottom_axis(group_rect);
			group_rect->axis(QCPAxis::atBottom)->setBasePen(coverage_bottom_pen);
			group_rect->setMarginGroup(QCP::msLeft | QCP::msRight, column_margins[j]);

			QVector<double> normalized_value = _Cs cast<QVector>(insertion_matrix.row(i));
			const int n_val = normalized_value.size();

			Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(n_val, true);

			for (int j = 1; j < n_val - 1; ++j) {
				if (normalized_value[j - 1] < threshold && normalized_value[j] < threshold && normalized_value[j + 1] < threshold) {
					filter[j] = false;
				}
			}

			QCPGraph* graph = draw_area->addGraph(group_rect->axis(QCPAxis::atBottom), group_rect->axis(QCPAxis::atLeft));

			graph->setPen(Qt::NoPen);
			graph->setBrush(QBrush(colors[i]));

			graph->setData(_Cs sliced(bottom_axis_value, filter), _Cs sliced(normalized_value, filter), true);

			_CpPatch set_range(group_rect, QCPRange(locs[j].start, locs[j].end + 1), QCPRange(0, 1));
		}
	}

	_Cp add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->level_name_,
		gs
	);
	this->draw_suite_->update(draw_area);
};
