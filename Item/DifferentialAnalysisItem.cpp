#include "DifferentialAnalysisItem.h"

#include "MatrixWindow.h"
#include "EnrichWorker.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CustomPlot.h"

#include "StringVector.h"

#include "FileWritingWorker.h"

void DifferentialAnalysisItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " | " + QString::number(this->data()->mat_.cols()) + " ]");
};

void DifferentialAnalysisItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Enrichment));

}

void DifferentialAnalysisItem::__show_this() {

	MatrixWindow::show_matrix(
		&this->data()->mat_,
		this->title_,
		this->signal_emitter_);
}

void DifferentialAnalysisItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Volcano Plot", s_volcano_plot);

	ADD_MAIN_ACTION("Heatmap Plot", s_heatmap_plot);

	ADD_MAIN_ACTION("Show Significant", s_show_significant);

	ADD_MAIN_ACTION("Extract Feature Names", s_extract_feature_names);

	ADD_MAIN_MENU("Enrich");

	ADD_ACTION("GO", "Enrich", s_enrich_go);
	ADD_ACTION("KEGG", "Enrich", s_enrich_kegg);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
};

void DifferentialAnalysisItem::s_heatmap_plot() {

	// TODO : should check whether the object has been changed

	if (this->attached_to(soap::VariableType::ChromVAR)) {
		G_WARN("Not Support in ChromVAR.");
		return;
	}

	Metadata* metadata{ nullptr };

	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		metadata = this->get_root<SingleCellRna>()->metadata();
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		metadata = this->get_root<SingleCellMultiome>()->metadata();
	}
	else {
		G_WARN("No Enough Data for heatmap visualization.");
		return;
	}

	QString metadata_name = this->data()->mat_.get_const_qstring_reference(METADATA_DE_METADATA_NAME).first();

	if (metadata->mat_.contains(metadata_name)) {
		auto data_type = metadata->mat_.data_type_.at(metadata_name);

		if (data_type != CustomMatrix::DataType::QStringFactor && data_type != CustomMatrix::DataType::IntegerFactor) {
			G_WARN("Metadata has been rewritten.");
			return;
		}
	}
	else {
		G_WARN("Metadata has been rewritten.");
		return;
	}

	auto metadata_content = metadata->mat_.get_qstring(metadata_name);

	auto& comparison_1 = this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1);

	auto comparison1s = custom::unique(comparison_1);

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Heatmap Plot Settings",
		{ "Comparisons", "log2 Fold Change(>):0.5", "Adjusted P Value(<):0.05" , 
		"Gene Number of Each Cluster:5", "Annotate"},
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::ComboBox},
		{comparison1s, { "Cell", "Gene"}}
	);

	if (settings.isEmpty()) {
		return;
	}

	comparison1s = simple_choice_to_list(settings[0]);

	if (comparison1s.isEmpty()) {
		return;
	}

	double fc_threshold = settings[1].toDouble();
	double p_threshold = settings[2].toDouble();
	int top_n = settings[3].toInt();

	if (top_n < 1) {
		G_WARN("illegal feature number.");
		return;
	}

	bool annotate_cell = settings[4] == "Cell";

	QStringList feature_list;
	QStringList comparison1_list;

	auto& pval = this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE);
	auto& fc_val = this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE);
	auto& feature_name = this->data()->mat_.get_const_qstring_reference(METADATA_DE_FEATURE_NAME);
	
	Eigen::ArrayX<bool> filter = custom::less_than(pval, p_threshold);
	filter *= custom::greater_than(fc_val, fc_threshold);
	
	for (auto&& comparison1 : comparison1s) {

		Eigen::ArrayX<bool> sub_filter = filter * custom::equal(comparison_1, comparison1);

		int n_feature_valid = sub_filter.count();

		if (n_feature_valid == 0) {
			continue;
		}

		QStringList valid_feature_names = custom::sliced(feature_name, sub_filter);

		if (n_feature_valid <= top_n) {
			feature_list << valid_feature_names;
			comparison1_list << QStringList(n_feature_valid, comparison1);
			continue;
		}

		auto valid_fcs = custom::sliced(fc_val, sub_filter);

		auto order = custom::order(valid_fcs, true);

		feature_list << custom::reordered(valid_feature_names, order(custom::seq_n(0, top_n)));
		comparison1_list << QStringList(top_n, comparison1);
	}

	if (feature_list.isEmpty()) {
		G_WARN("No valid features found.");
		return;
	}

	QVector<int> cell_order;
	QList<QPair<int, QString>> cell_types;
	for (auto&& comparison1 : comparison1s) {

		auto sub_order = custom::match(metadata_content, comparison1);

		if (sub_order.isEmpty()) {
			continue;
		}

		cell_order << sub_order;
		cell_types << qMakePair( sub_order.size(), comparison1);
	}

	if (cell_order.isEmpty()) {
		G_WARN("Invalid metadata.");
		return;
	}

	int n_cell = cell_order.size();

	Eigen::MatrixXd heatmat;

	if (this->attached_to(soap::VariableType::DataField)) {

		auto df = this->trace_back<DataField>(1);

		auto normalized = df->normalized();

		if (normalized == nullptr) {
			G_WARN("Cannot find normalized data.");
			return;
		}

		auto feature_index = custom::index_of(feature_list, normalized->rownames_);

		auto valid_elements = custom::not_equal(feature_index, -1);
		if (valid_elements.count() == 0) {
			G_WARN("No valid feature found in normalized data.");
			return;
		}

		feature_index.removeAll(-1);

		comparison1_list = custom::sliced(comparison1_list, valid_elements);

		int n_feature = feature_index.size();
		heatmat.resize(n_feature, n_cell);

	#pragma omp parallel for
		for (int i = 0; i < n_feature; ++i) {
			Eigen::ArrayXd f = normalized->mat_.row(feature_index[i]);
			heatmat.row(i) = custom::reordered(f, cell_order);
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellRna)) {

		SingleCellRna* single_cell_rna = this->trace_back<SingleCellRna>(1);

		auto normalized = single_cell_rna->normalized();

		if (normalized == nullptr) {
			G_WARN("Cannot find normalized data.");
			return;
		}

		auto feature_index = custom::index_of(feature_list, normalized->rownames_);

		auto valid_elements = custom::not_equal(feature_index, -1);
		if (valid_elements.count() == 0) {
			G_WARN("No valid feature found in normalized data.");
			return;
		}

		feature_index.removeAll(-1);

		comparison1_list = custom::sliced(comparison1_list, valid_elements);

		int n_feature = feature_index.size();
		heatmat.resize(n_feature, n_cell);

	#pragma omp parallel for
		for (int i = 0; i < n_feature; ++i) {
			Eigen::ArrayXd f = normalized->mat_.row(feature_index[i]);
			heatmat.row(i) = custom::reordered(f, cell_order);
		}
	}
	else {
		G_WARN("Cannot find normalized data.");
		return;
	}

	int n_feature = comparison1_list.size();

	heatmat.colwise() -= heatmat.rowwise().mean();

#pragma omp parallel for
	for (int i = 0; i < n_feature; ++i) {
		double sd = custom::sd(heatmat.row(i));

		if (sd != 0.0) {
			heatmat.row(i) /= sd;
		}
	}

	const auto& gs = this->draw_suite_->graph_settings_;

	QCustomPlot* draw_area = custom_plot::initialize_plot(gs);

	QCPLayoutGrid* main_layout = new QCPLayoutGrid;
	draw_area->plotLayout()->addElement(0, 0, main_layout);

	if (annotate_cell) {

		QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);
		QCPAxisRect* top_legend = new QCPAxisRect(draw_area, true);

		QCPLayoutGrid* legend_layout = new QCPLayoutGrid;

		QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);
		axis_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
		top_legend->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);

		main_layout->addElement(0, 0, top_legend);
		main_layout->addElement(1, 0, axis_rect);
		main_layout->addElement(1, 1, legend_layout);

		legend_layout = custom_plot::patch::set_legend_layout(draw_area, legend_layout);

		auto pal = gs.palette_map(comparison1_list);

		int now_start = 0;

		for (auto&& p : cell_types) {
			custom_plot::patch::rectangle_borderless(
				draw_area, top_legend, now_start, -3.0, p.first, 3.0, pal[p.second]
			);

			custom_plot::patch::add_label(
				draw_area,
				top_legend,
				p.second,
				now_start + p.first / 2,
				1.0,
				gs.get_bottom_label_font(),
				Qt::AlignBottom | Qt::AlignHCenter);

			now_start += p.first;
		}

		custom_plot::patch::remove_left_bottom_axis(top_legend);

		custom_plot::patch::set_fixed_height(top_legend,
			std::ceil(custom_plot::utility::get_max_text_height(custom::unique(comparison1_list), gs.get_left_label_font()) * 1.4));

		custom_plot::patch::set_range(top_legend, { 0.0, (double)n_cell }, {-3.0, 11.0});

		QCPColorMap* heatmap = new QCPColorMap(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		heatmap->data()->setSize(n_cell, n_feature);
		if (n_feature > 1) {
			heatmap->data()->setRange(QCPRange(0, n_cell - 1), QCPRange(0, n_feature - 1));
		}
		else {
			heatmap->data()->setRange(QCPRange(0, n_cell - 1), QCPRange(-1, 1));
		}
		custom_plot::patch::remove_left_bottom_axis(axis_rect);

		for (int i = 0; i < n_cell; ++i) {
			for (int j = 0; j < n_feature; ++j) {
				heatmap->data()->setCell(i, j, heatmat(j, i));
			}
		}
		heatmap->setInterpolate(false);
		heatmap->setTightBoundary(false);
		custom_plot::patch::set_range(axis_rect, QCPRange(-0.5, n_cell - 0.5), QCPRange(-0.5, n_feature - 0.5));

		double max_val = heatmat.cwiseAbs().maxCoeff();

		QCPColorGradient gradient;
		gradient.setColorStopAt(0.5, custom_plot::color::aquamarine3);
		gradient.setColorStopAt(1.0, custom_plot::color::gold);
		gradient.setColorStopAt(0.0, custom_plot::color::navy);

		heatmap->setGradient(gradient);
		heatmap->setDataRange({ -1.0, 1.0 });

		custom_plot::add_gradient_legend(
			draw_area,
			legend_layout,
			-1.0,
			1.0,
			"Expression",
			gs,
			custom_plot::color::navy,
			custom_plot::color::aquamarine3,
			custom_plot::color::gold,
			"Low",
			"High"
		);

		this->draw_suite_->update(draw_area);
	}
	else {

		QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);
		QCPAxisRect* left_legend = new QCPAxisRect(draw_area, true);

		QCPLayoutGrid* legend_layout = new QCPLayoutGrid;

		QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);
		axis_rect->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);
		left_legend->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);
		legend_layout->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);

		main_layout->addElement(0, 0, left_legend);
		main_layout->addElement(0, 1, axis_rect);
		main_layout->addElement(0, 2, legend_layout);

		legend_layout = custom_plot::patch::set_legend_layout(draw_area, legend_layout);

		auto pal = gs.palette_map(comparison1_list);

		auto* now_comparison1 = &comparison1_list.first();
		int now_start = 0;
		for (int i = 1; i < n_feature; ++i) {

			if (comparison1_list[i] != *now_comparison1) {

				custom_plot::patch::rectangle_borderless(
					draw_area, left_legend, 0, now_start, 10, i - now_start, pal[*now_comparison1]
				);

				custom_plot::patch::add_label(
					draw_area,
					left_legend,
					*now_comparison1,
					-0.5,
					(now_start + i) / 2.0,
					gs.get_left_label_font(),
					Qt::AlignRight | Qt::AlignVCenter);

				now_comparison1 = &comparison1_list[i];
				now_start = i;
			}
		}

		custom_plot::patch::rectangle_borderless(
			draw_area, left_legend, 0, now_start, 10, n_feature - now_start, pal[*now_comparison1]
		);

		custom_plot::patch::add_label(
			draw_area,
			left_legend,
			*now_comparison1,
			-0.5,
			(now_start + n_feature) / 2.0,
			gs.get_left_label_font(),
			Qt::AlignRight | Qt::AlignVCenter);

		custom_plot::patch::remove_left_bottom_axis(left_legend);
		custom_plot::patch::set_fixed_width(left_legend, 
			std::ceil(custom_plot::utility::get_max_text_width(custom::unique(comparison1_list), gs.get_left_label_font()) * 2.11));

		custom_plot::patch::set_range(left_legend, { -10.0, 10.0 }, { 0.0, (double)n_feature });

		QCPColorMap* heatmap = new QCPColorMap(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		heatmap->data()->setSize(n_cell, n_feature);
		if (n_feature > 1) {
			heatmap->data()->setRange(QCPRange(0, n_cell - 1), QCPRange(0, n_feature - 1));
		}
		else {
			heatmap->data()->setRange(QCPRange(0, n_cell - 1), QCPRange(-1, 1));
		}
		custom_plot::patch::remove_left_bottom_axis(axis_rect);

		for (int i = 0; i < n_cell; ++i) {
			for (int j = 0; j < n_feature; ++j) {
				heatmap->data()->setCell(i, j, heatmat(j, i));
			}
		}
		heatmap->setInterpolate(false);
		heatmap->setTightBoundary(false);
		custom_plot::patch::set_range(axis_rect, QCPRange(-0.5, n_cell - 0.5), QCPRange(-0.5, n_feature - 0.5));

		double max_val = heatmat.cwiseAbs().maxCoeff();

		QCPColorGradient gradient;
		gradient.setColorStopAt(0.5, custom_plot::color::aquamarine3);
		gradient.setColorStopAt(1.0, custom_plot::color::gold);
		gradient.setColorStopAt(0.0, custom_plot::color::navy);

		heatmap->setGradient(gradient);
		heatmap->setDataRange({ -1.0, 1.0 });

		custom_plot::add_gradient_legend(
			draw_area,
			legend_layout,
			-1.0,
			1.0,
			"Expression",			
			gs,
			custom_plot::color::navy,
			custom_plot::color::aquamarine3,
			custom_plot::color::gold,
			"Low",
			"High"
		);

		this->draw_suite_->update(draw_area);
	}
};

void DifferentialAnalysisItem::s_volcano_plot() {

	if (this->attached_to(soap::VariableType::ChromVAR)) {
		G_WARN("Not Support in ChromVAR.");
		return;
	}

	auto comparisons = custom::paste(this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1),
		this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_2), " vs. ");

	auto comp = custom::unique(comparisons);

	QStringList standard = CommonDialog::get_response(
		this->signal_emitter_,
		"Volcano Plot Setting",
		{ "Adjusted P Value(<):0.05", "log2 Fold Change(Absolute Value >):1", "Comparison", "Max Trans P:100", "Axis Style", "Show border:yes" },
		{ soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox
		, soap::InputStyle::IntegerLineEdit, soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton},
		{ comp, { "Simple", "Border" }}
	);
	if (standard.isEmpty())return;

	double feature_p_threshold = standard[0].toDouble();
	double fold_change_threshold = standard[1].toDouble();
	auto filter = custom::equal(comparisons, standard[2]);
	if (feature_p_threshold <= 0) {
		G_LOG("Invalid P Value.");
		return;
	}
	if (fold_change_threshold <= 0) {
		G_LOG("Invalid Threshold Value.");
		return;
	}
	feature_p_threshold = -log10(feature_p_threshold);
	int maximum_transformed_p = standard[3].toInt();
	if (maximum_transformed_p < 1) {
		G_WARN("Illegal plot setting for transformed P value!");
		return;
	}
	double minimum_p = std::pow(10, -maximum_transformed_p);

	int n_point = filter.count();
	Eigen::ArrayXd transformed_p(n_point);
	Eigen::ArrayXd p_adjusted = custom::cast<Eigen::ArrayX>(custom::sliced(this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), filter));

	Eigen::ArrayXd fold_change = custom::cast<Eigen::ArrayX>(custom::sliced(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), filter));

	for (int i = 0; i < n_point; ++i) {
		if (p_adjusted[i] < minimum_p) {
			transformed_p[i] = maximum_transformed_p;
		}
		else {
			transformed_p[i] = -log10(p_adjusted[i]);
		}
	}

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
	
	QString axis_style = standard[4];
	if (axis_style == "Simple") {
		custom_plot::set_simple_axis(axis_rect, "log<sub>2</sub>(Fold Change)", "-log<sub>10</sub>(P<sub>adj</sub>)", gs);
	}
	else if (axis_style == "Border") {
		custom_plot::patch::set_border_only(axis_rect, Qt::black, 3);
		custom_plot::set_left_title(axis_rect, "-log<sub>10</sub>(P<sub>adj</sub>)", gs, true);
		custom_plot::set_bottom_title(axis_rect, "log<sub>2</sub>(Fold Change)", gs, true);
		axis_rect->axis(QCPAxis::atLeft)->setTickLabelFont(gs.get_left_label_font());
		axis_rect->axis(QCPAxis::atBottom)->setTickLabelFont(gs.get_bottom_label_font());
	}

	bool show_border = switch_to_bool(standard[5]);

	auto [x_range, y_range] = custom_plot::utility::get_range(fold_change, transformed_p);
	custom_plot::patch::set_range(axis_rect, x_range, y_range);

	if (show_border) {
		custom_plot::patch::line(draw_area, axis_rect, QVector<double>{fold_change_threshold, fold_change_threshold}, QVector<double>{y_range.lower, y_range.upper}, Qt::black, 2, Qt::DashLine);
		custom_plot::patch::line(draw_area, axis_rect, QVector<double>{-fold_change_threshold, -fold_change_threshold}, QVector<double>{y_range.lower, y_range.upper}, Qt::black, 2, Qt::DashLine);
		custom_plot::patch::line(draw_area, axis_rect, QVector<double>{ x_range.lower, x_range.upper}, QVector<double>{-std::log10(feature_p_threshold), -std::log10(feature_p_threshold)}, Qt::black, 2, Qt::DashLine);
	}

	QStringList values(n_point);

	for (int i = 0; i < n_point; ++i) {
		values[i] = transformed_p[i] > feature_p_threshold ? (fold_change[i] > fold_change_threshold ? "Red" : (fold_change[i] < -fold_change_threshold ? "Blue" : "Grey")) : "Grey";
	}
	custom_plot::patch::scatter_category(
		draw_area, 
		axis_rect, 
		fold_change, 
		transformed_p, 
		values,
		QStringList() << "Red" << "Blue" << "Grey", 
		QList<QColor>{ gs.get_gradient_high_color(custom_plot::color::firebrick3), gs.get_gradient_low_color(custom_plot::color::navy), gs.get_gradient_middle_color(custom_plot::color::gray)},
		gs.get_scatter_point_size());
	this->draw_suite_->update(draw_area);
};

void DifferentialAnalysisItem::s_extract_feature_names() {

	auto comparisons = custom::paste(this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1),
		this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_2), " vs. ");

	auto comp = custom::unique(comparisons);

	QString value_filter_string{ "log2 Fold Change(Absolute Value >):1" };
	QString value_name{ METADATA_DE_LOG2_FOLD_CHANGE };

	if (this->data()->data_type_ == DifferentialAnalysis::DataType::ChromVAR) {
		value_filter_string = "Z Score Difference Threshold:2.0";
		value_name = METADATA_DE_Z_SCORE_DIFFERENCE;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Setting For Significance",
		{ "P Value(<):0.05", value_filter_string, "Type", "All Comparison:no", "Comparison", "Top N:0"},
		{ soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox, 
		soap::InputStyle::SwitchButton, soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit },
		{ { "All", "Upregulated", "Downregulated" }, comp }
	);
	if (settings.isEmpty())return;

	double feature_p_threshold = settings[0].toDouble();
	double value_threshold = settings[1].toDouble();
	if (feature_p_threshold <= 0 || feature_p_threshold > 1) {
		G_LOG("Invalid p value! Reset to 0.05");
		feature_p_threshold = 0.05;
	}

	Eigen::ArrayX<bool> filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), feature_p_threshold);

	QString regulate_type = settings[2];
	if (regulate_type == "Upregulated") {
		filter *= custom::greater_than(this->data()->mat_.get_const_double_reference(value_name), 0.);
	}
	else if (regulate_type == "Downregulated") {
		filter *= custom::less_than(this->data()->mat_.get_const_double_reference(value_name), 0.);
	}

	bool all_comparison = switch_to_bool(settings[3]);
	if (!all_comparison) {
		filter *= custom::equal(comparisons, settings[4]);
	}

	filter *= custom::greater_equal(custom::abs(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE)), std::abs(value_threshold));
	if (filter.count() == 0) {
		G_NOTICE("No feature found.");
		return;
	}

	int top_n = settings[5].toInt();
	if (top_n < 0) {
		G_WARN("Invalid settings : top N");
		return;
	}

	if (all_comparison) {

		QVector<int> order;

		for (auto&& c : comp) {
			Eigen::ArrayX<bool> sub_filter = filter * custom::equal(comparisons, c);
			auto sub_vals = custom::sliced(this->data()->mat_.get_const_double_reference(value_name), sub_filter);
			auto sub_order = custom::order(sub_vals, true);

			if (top_n > 0 && sub_order.size() > top_n) {
				sub_order = sub_order.segment(0, top_n).eval();
			}

			order << custom::cast<QVector>(custom::reordered(custom::which(sub_filter), sub_order));
		}

		QStringList feature_names = custom::reordered(this->data()->mat_.get_const_qstring_reference(METADATA_DE_FEATURE_NAME), order);
		StringVector* sv = new StringVector(feature_names);

		this->signal_emitter_->x_data_create_soon(sv, soap::VariableType::StringVector, "Extracted Feature Names");
	}
	else {

		auto vals = custom::sliced(this->data()->mat_.get_const_double_reference(value_name), filter);
		auto order = custom::order(vals, true);
		QStringList feature_names = custom::reordered(this->data()->mat_.get_const_qstring_reference(METADATA_DE_FEATURE_NAME), custom::reordered(custom::which(filter), order));
		
		if (top_n > 0 && feature_names.size() > top_n) {
			feature_names = feature_names.sliced(0, top_n);
		}
		
		StringVector* sv = new StringVector(feature_names);

		this->signal_emitter_->x_data_create_soon(sv, soap::VariableType::StringVector, "Extracted Feature Names");
	}

};

void DifferentialAnalysisItem::s_show_significant() {

	auto comparisons = custom::paste(this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1),
		this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_2), " vs. ");

	auto comp = custom::unique(comparisons);

	QString value_filter_string{ "log2 Fold Change(Absolute Value >):1" };
	QString value_name{ METADATA_DE_LOG2_FOLD_CHANGE };

	if (this->data()->data_type_ == DifferentialAnalysis::DataType::ChromVAR) {
		value_filter_string = "Z Score Difference Threshold:2.0";
		value_name = METADATA_DE_Z_SCORE_DIFFERENCE;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Setting For Significance",
		{ "P Value(<):0.05", value_filter_string, "Type", "All Comparison:no","Comparison" },
		{soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox, 
		soap::InputStyle::SwitchButton, soap::InputStyle::ComboBox},
		{{ "All", "Upregulated", "Downregulated" }, comp }
	);
	if (settings.isEmpty())return;

	double feature_p_threshold = settings[0].toDouble();
	double value_threshold = settings[1].toDouble();
	if (feature_p_threshold <= 0 || feature_p_threshold > 1) {
		G_LOG("Invalid p value! Reset to 0.05");
		feature_p_threshold = 0.05;
	}

	Eigen::ArrayX<bool> filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), feature_p_threshold);
	
	QString regulate_type = settings[2];
	if (regulate_type == "Upregulated") {
		filter *= custom::greater_than(this->data()->mat_.get_const_double_reference(value_name), 0.);
	}
	else if (regulate_type == "Downregulated") {
		filter *= custom::less_than(this->data()->mat_.get_const_double_reference(value_name), 0.);
	}

	bool all_comparison = switch_to_bool(settings[3]);
	if (!all_comparison) {
		filter *= custom::equal(comparisons, settings[4]);
	}

	filter *= custom::greater_equal(custom::abs(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE)), std::abs(value_threshold));
	if (filter.count() == 0) {
		G_NOTICE("No significant differential feature found.");
		return;
	}

	auto vals = custom::sliced(this->data()->mat_.get_const_double_reference(value_name), filter);
	auto order = custom::order(vals, true);

	auto tmp = this->data()->mat_.row_sliced(filter).row_reordered(order);
	
	MatrixWindow::show_matrix(
		&tmp,
		"Significant Result",
		this->signal_emitter_);
};

void DifferentialAnalysisItem::s_enrich_go() {

	auto data_type = this->data()->data_type_;

	if (data_type != DifferentialAnalysis::DataType::Gene && data_type != DifferentialAnalysis::DataType::GeneActivity) {
		G_WARN("The Data Type of Differential Analysis is not valid.");
		return;
	}

	G_GETLOCK;

	soap::Species species;
	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		species = (this->get_root<SingleCellRna>())->species_;
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		species = (this->get_root<SingleCellMultiome>())->species_;
	}
	else {
		QStringList res = CommonDialog::get_response(
			this->signal_emitter_,
			"Species Setting",
			{ "Species" },
			{ soap::InputStyle::ComboBox},
			{ { "Human", "Mouse" }}
		);
		if (res.isEmpty()) {
			G_UNLOCK;
			return;
		}
		species = res[0] == "Human" ? soap::Species::Human : soap::Species::Mouse;
	}
	if (species == soap::Species::Undefined) {
		G_LOG("Only human and mouse feature enrichment are supported");
		G_UNLOCK;
		return;
	}

	auto comparisons = custom::paste(this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1),
		this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_2), " vs. ");

	auto comp = custom::unique(comparisons);

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Enrichment Setting",
		{ "Type", "P threshold (feature):0.05", "Log2 fold change:0.5",
		"Ontology", "P adjust method", "P threshold (pathway):0.05", "Comparison" },
		{soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit
		, soap::InputStyle::ComboBox, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox},
		{{ "Upregulated", "Downregulated", "ALL"}, {"BP", "MF", "CC", "ALL"},
		{ "FDR", "Bonferroni" }, comp }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	double feature_p_threshold = settings[1].toDouble(), fold_change_threshold = settings[2].toDouble(), pathway_p_threshold = settings[5].toDouble();
	if (feature_p_threshold <= 0 || feature_p_threshold > 1) {
		G_LOG("Invalid p value to filter features! Reset to 0.05");
		feature_p_threshold = 0.05;
	}
	if (pathway_p_threshold <= 0 || pathway_p_threshold > 1) {
		G_LOG("Invalid p value to filter pathways! Reset to 0.05");
		pathway_p_threshold = 0.05;
	}

	Eigen::ArrayX<bool> filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), feature_p_threshold);
	filter *= custom::greater_equal(custom::abs(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE)), std::abs(fold_change_threshold));

	QString regulate_type = settings[0];
	if (regulate_type == "Upregulated") {
		filter *= custom::greater_than(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), 0.);
	}
	else if (regulate_type == "Downregulated") {
		filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), 0.);
	}
	filter *= custom::equal(comparisons, settings[6]);
	if (filter.count() == 0) {
		G_LOG("No feature meets requirement.");
		G_UNLOCK;
		return;
	}
	QStringList feature_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_DE_FEATURE_NAME), filter);
	EnrichWorker* worker = new EnrichWorker("GO", settings[6] + " " + regulate_type, feature_names, settings[3], species, settings[4], pathway_p_threshold);
	G_LINK_WORKER_THREAD(EnrichWorker, x_enrichment_ready, DifferentialAnalysisItem, s_receive_enrichment);
};

void DifferentialAnalysisItem::s_enrich_kegg() {

	auto data_type = this->data()->data_type_;

	if (data_type != DifferentialAnalysis::DataType::Gene && data_type != DifferentialAnalysis::DataType::GeneActivity) {
		G_WARN("The Data Type of Differential Analysis is not valid.");
		return;
	}

	G_GETLOCK;
	
	soap::Species species;
	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		species = (this->get_root<SingleCellRna>())->species_;
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		species = (this->get_root<SingleCellMultiome>())->species_;
	}
	else {
		QStringList res = CommonDialog::get_response(
			this->signal_emitter_,
			"Species Setting",
			{ "Species" },
			{soap::InputStyle::ComboBox},
			{ { "Human", "Mouse" }}
		);

		if (res.isEmpty()) {
			G_UNLOCK;
			return;
		}
		species = res[0] == "Human" ? soap::Species::Human : soap::Species::Mouse;
	}

	if (species == soap::Species::Undefined) {
		G_LOG("Only human and mouse feature enrichment are supported");
		G_UNLOCK;
		return;
	}

	auto comparisons = custom::paste(this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1),
		this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_2), " vs. ");

	auto comp = custom::unique(comparisons);

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Enrichment Setting",
		{ "Type", "P threshold (filter features):0.05", "Log2 fold change:0.5",
		"P adjust method", "P threshold (filter pathways):0.05", "Comparison" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit
		, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox},
		{ { "Upregulated", "Downregulated", "ALL"}, {"FDR", "Bonferroni"}, comp }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	double feature_p_threshold = settings[1].toDouble(), fold_change_threshold = settings[2].toDouble(), pathway_p_threshold = settings[4].toDouble();
	if (feature_p_threshold <= 0 || feature_p_threshold > 1) {
		G_LOG("Invalid p value to filter features! Reset to 0.05");
		feature_p_threshold = 0.05;
	}
	if (pathway_p_threshold <= 0 || pathway_p_threshold > 1) {
		G_LOG("Invalid p value to filter pathways! Reset to 0.05");
		pathway_p_threshold = 0.05;
	}
	Eigen::ArrayX<bool> filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), feature_p_threshold);
	filter *= custom::greater_equal(custom::abs(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE)), std::abs(fold_change_threshold));

	QString regulate_type = settings[0];
	if (regulate_type == "Upregulated") {
		filter *= custom::greater_than(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), 0.);
	}
	else if (regulate_type == "Downregulated") {
		filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), 0.);
	}

	filter *= custom::equal(comparisons, settings[5]);

	if (filter.count() == 0) {
		G_LOG("No feature meets requirement.");
		G_UNLOCK;
		return;
	}
	QStringList feature_names = custom::sliced(this->data()->mat_.rownames_, filter);

	EnrichWorker* worker = new EnrichWorker(
		"KEGG", 
		settings[5] + " " + regulate_type,
		feature_names, 
		"ALL", 
		species,
		settings[3],
		pathway_p_threshold
	);
	
	G_LINK_WORKER_THREAD(EnrichWorker, x_enrichment_ready, DifferentialAnalysisItem, s_receive_enrichment)
};


void DifferentialAnalysisItem::s_enrich_motif() {

	auto data_type = this->data()->data_type_;

	if (data_type != DifferentialAnalysis::DataType::Peak) {
		G_WARN("The Data Type of Differential Analysis is not valid.");
		return;
	}

	G_GETLOCK;

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This data is not attached to single cell multiome data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto motif_position = single_cell_multiome->motif_position();
	if (motif_position == nullptr) {
		G_WARN("Motif is not located.");
		G_UNLOCK;
		return;
	}

	soap::Species species = single_cell_multiome->species_;

	if (species != soap::Species::Human) {
		G_LOG("Only human gene enrichment are supported");
		G_UNLOCK;
		return;
	}

	auto comparisons = custom::paste(this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_1),
		this->data()->mat_.get_const_qstring_reference(METADATA_DE_COMPARISON_2), " vs. ");

	auto comp = custom::unique(comparisons);

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Enrichment Setting",
		{ "Type", "P threshold (filter genes):0.05", "Log2 fold change:0.5",
		"P adjust method", "P threshold (filter motifs):0.05", "Comparison" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit
		, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox},
		{ { "Upregulated", "Downregulated", "ALL"}, {"FDR", "Bonferroni"}, comp}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	double peak_p_threshold = settings[1].toDouble(), fold_change_threshold = settings[2].toDouble(), motif_p_threshold = settings[4].toDouble();
	if (peak_p_threshold <= 0 || peak_p_threshold > 1) {
		G_LOG("Invalid p value to filter genes! Reset to 0.05");
		peak_p_threshold = 0.05;
	}

	if (motif_p_threshold <= 0 || motif_p_threshold > 1) {
		G_LOG("Invalid p value to filter pathways! Reset to 0.05");
		motif_p_threshold = 0.05;
	}

	Eigen::ArrayX<bool> filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), peak_p_threshold);
	filter *= custom::greater_equal(custom::abs(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE)), std::abs(fold_change_threshold));

	QString regulate_type = settings[0];
	if (regulate_type == "Upregulated") {
		filter *= custom::greater_than(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), 0.);
	}
	else if (regulate_type == "Downregulated") {
		filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_DE_LOG2_FOLD_CHANGE), 0.);
	}

	filter *= custom::equal(comparisons, settings[5]);

	if (filter.count() == 0) {
		G_LOG("No peak meets requirement.");
		G_UNLOCK;
		return;
	}

	QStringList peak_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_DE_FEATURE_NAME), filter);
	
	EnrichWorker* worker = new EnrichWorker(
		peak_names,
		motif_position,
		regulate_type,
		settings[3],
		motif_p_threshold
	);
	G_LINK_WORKER_THREAD(EnrichWorker, x_enrichment_ready, DifferentialAnalysisItem, s_receive_enrichment);
};

void DifferentialAnalysisItem::s_receive_enrichment(const CustomMatrix& matrix, QString name) {

	QString new_title = this->signal_emitter_->get_unique_name(name);
	DATA_SUBMODULES(Enrichment)[new_title] = Enrichment(matrix);

	EnrichmentItem* item = new EnrichmentItem(
		new_title,
		this->index_tree_,
		&DATA_SUBMODULES(Enrichment)[new_title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Enrichment finished");
};
