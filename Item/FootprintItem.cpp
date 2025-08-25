#include "FootprintItem.h"

#include "GeneNameItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "FileWritingWorker.h"
#include "ItemIOWorker.h"
#include "CalculateCountsByGenomicRangeWorker.h"
#include "CustomPlot.h"

#include "WilcoxTest.h"

void FootprintItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->insertion_matrix_.mat_,
		this->data()->insertion_matrix_.rownames_,
		this->data()->insertion_matrix_.colnames_,
		this->title_,
		this->signal_emitter_);
}

void FootprintItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Show", s_show);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void FootprintItem::s_show() {


	Metadata* metadata{ nullptr };

	QMap<QString, QStringList> map;

	if (this->stem_from(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();
		metadata = single_cell_multiome->metadata();
		map = metadata->mat_.get_factor_information();

		if (map.isEmpty()) {
			G_NOTICE("No suitable metadata detected.");
			return;
		}
	}
	else if (this->stem_from(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->get_root<SingleCellAtac>();
		metadata = single_cell_atac->metadata();

		map = metadata->mat_.get_factor_information();

		if (map.isEmpty()) {
			G_NOTICE("No suitable metadata detected.");
			return;
		}
	}
	else {
		G_WARN("Illegal data status.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Footprint Visualization Settings",
		{ "Factor:Factor", "Show P Value:no", "Show in group:no", "Show Group"},
		{ soap::InputStyle::FactorChoice, soap::InputStyle::SwitchButton,
		soap::InputStyle::SwitchButton, soap::InputStyle::FactorChoice},
		QList<QStringList>(), { map, map}
	);
	if (settings.isEmpty()) {
		return;
	}

	bool show_p_value = switch_to_bool(settings[1]);

	auto [factor_name, levels] = factor_choice_to_pair(settings[0]);
	if (levels.isEmpty())return;

	bool show_in_group = switch_to_bool(settings[2]);

	auto [show_group_name, show_levels] = factor_choice_to_pair(settings[3]);
	if (show_in_group) {
		if (show_levels.isEmpty()) {
			return;
		}
	}

	QStringList factors = metadata->mat_.get_qstring(factor_name);
	QStringList show_group = metadata->mat_.get_qstring(show_group_name);
	auto show_filter = custom::in(show_group, show_levels);

	const int nrow = this->data()->insertion_matrix_.mat_.rows(), ncol = this->data()->insertion_matrix_.mat_.cols();

	if (nrow != factors.size() || ncol <= 495) {
		G_WARN("Illegal Footprint Object!");
		return;
	}
	Eigen::ArrayXd flank_mean = this->data()->insertion_matrix_.mat_.block(0, 0, nrow, 50).rowwise().sum().cast<double>() + this->data()->insertion_matrix_.mat_.block(0, ncol - 50, nrow, 50).rowwise().sum().cast<double>();
	flank_mean /= 100;
	double all_mean = flank_mean.mean();
	if (all_mean == 0) {
		G_WARN("Illegal Footprint Data with too many zeros!");
		return;
	}
	for (int i = 0; i < nrow; ++i) {
		if (flank_mean[i] == 0) {
			flank_mean[i] = all_mean;
		}
	}
	Eigen::MatrixXd normalized = this->data()->insertion_matrix_.mat_.cast<double>();
	normalized.array().colwise() /= flank_mean; 

	QList<Eigen::ArrayXd> locations;
	QStringList levels_use;
	for (const auto& factor : levels) {
		auto index = show_in_group ? custom::which(custom::equal(factors, factor) * show_filter) : custom::match(factors, factor);

		if (index.isEmpty()) {
			continue;
		}

		Eigen::MatrixXd sub_matrix = normalized(index, Eigen::all);

		Eigen::ArrayXd means = sub_matrix.colwise().mean();

		levels_use << factor;
		locations << means - custom::cast<Eigen::ArrayX>(this->data()->expected_insertions_);
	}

	if (locations.isEmpty()) {
		G_WARN("No Factor is included.");
		return;
	}

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, left_layout, legend_layout] = custom_plot::prepare_lg_lg(gs);
	
	QCPAxisRect* observation_axisrect = new QCPAxisRect(draw_area), * expect_axisrect = new QCPAxisRect(draw_area);
	left_layout->addElement(0, 0, observation_axisrect);
	left_layout->addElement(1, 0, expect_axisrect);
	left_layout->setRowStretchFactors(QList<double>() << 4 << 1);

	int ncolor = levels_use.size();
	auto colors = gs.palette(levels_use);

	QVector<double> x_axis = custom::cast<double>(this->data()->insertion_matrix_.colnames_);

	for (int i = 0; i < ncolor; ++i) {
		custom_plot::patch::line(draw_area, observation_axisrect, x_axis, custom::cast<QVector>(locations[i]), colors[i], 2);
	}

	if (show_p_value) {
		if (locations.size() != 2) {
			G_NOTICE("P value can only be calculated between two groups.");
		}
		else {
			double p_value = wilcox_test(locations[0], locations[1], true);
			QCPItemText* p_value_label = new QCPItemText(draw_area);
			p_value_label->setClipToAxisRect(false);
			p_value_label->position->setAxisRect(observation_axisrect);
			p_value_label->position->setAxes(observation_axisrect->axis(QCPAxis::atBottom), observation_axisrect->axis(QCPAxis::atLeft));
			p_value_label->position->setType(QCPItemPosition::ptAxisRectRatio);
			p_value_label->setPositionAlignment(Qt::AlignTop | Qt::AlignRight);
			p_value_label->position->setCoords(1, 0);
			p_value_label->setText("P Value : " + QString::number(p_value, 'f'));
			p_value_label->setFont(QFont("Arial", 15, 500, true));
			p_value_label->setColor(Qt::black);

			G_NOTICE("The Test Method for Footprint is not finally determined. Wilcoxon Paired Sign Rank Test is now used.");
		}
	}

	double min_ob = std::ranges::min(custom::sapply(locations, [](auto&& arr) {return arr.minCoeff(); }));
	double max_ob = std::ranges::max(custom::sapply(locations, [](auto&& arr) {return arr.maxCoeff(); }));

	auto [min_x, max_x] = std::ranges::minmax(x_axis);

	custom_plot::patch::set_range(observation_axisrect, QCPRange(min_x, max_x), QCPRange(1.1 * min_ob - 0.1 * max_ob, 1.1 * max_ob - 0.1 * min_ob));
	custom_plot::set_left_title(observation_axisrect, "Insertion Frequency", gs);

	QPen axis_pen(Qt::black);
	axis_pen.setWidth(3);

	custom_plot::patch::remove_bottom_axis(observation_axisrect);

	observation_axisrect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	observation_axisrect->axis(QCPAxis::atLeft)->setSubTickPen(Qt::NoPen);
	observation_axisrect->axis(QCPAxis::atLeft)->setBasePen(axis_pen);

	auto [min_expect, max_expect] = std::ranges::minmax(this->data()->expected_insertions_);
	custom_plot::patch::set_range(expect_axisrect, QCPRange(min_x, max_x), QCPRange(1.1 * min_expect - 0.1 * max_expect, 1.1 * max_expect - 0.1 * min_expect));
	custom_plot::patch::line(draw_area, expect_axisrect, x_axis, this->data()->expected_insertions_, Qt::black, 2);

	QCPMarginGroup* edge_margin_group = new QCPMarginGroup(draw_area);
	observation_axisrect->setMarginGroup(QCP::msLeft | QCP::msRight, edge_margin_group);
	expect_axisrect->setMarginGroup(QCP::msLeft | QCP::msRight, edge_margin_group);
	
	custom_plot::patch::remove_grid(expect_axisrect);
	expect_axisrect->axis(QCPAxis::atBottom)->setBasePen(axis_pen);
	expect_axisrect->axis(QCPAxis::atBottom)->setSubTickPen(Qt::NoPen);
	expect_axisrect->axis(QCPAxis::atLeft)->setBasePen(axis_pen);
	expect_axisrect->axis(QCPAxis::atLeft)->setSubTickPen(Qt::NoPen);
	expect_axisrect->axis(QCPAxis::atLeft)->ticker()->setTickCount(3);

	custom_plot::set_left_title(expect_axisrect, "Expected", gs);
	custom_plot::set_bottom_title(expect_axisrect, "Relative Position (bp)", gs, true);

	custom_plot::add_round_legend(draw_area, legend_layout, levels_use, colors, factor_name, gs);

	custom_plot::add_title(draw_area, this->data()->motif_.motif_name_, gs);
	this->draw_suite_->update(draw_area);
}

void FootprintItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(GeneName));
}

void FootprintItem::s_receive_correlated_gene(QString item_name, QStringList gene_names) {

	item_name = this->signal_emitter_->get_unique_name(item_name);

	DATA_SUBMODULES(GeneName)[item_name] = GeneName(gene_names);

	GeneNameItem* item = new GeneNameItem(
		item_name,
		this->index_tree_,
		&DATA_SUBMODULES(GeneName)[item_name],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);
};