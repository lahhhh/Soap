#include "PandoItem.h"

#include "MatrixWindow.h"
#include "CustomPlot.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"

#include "GeneNameItem.h"
#include "EmbeddingItem.h"

void PandoItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " ]");
}

void PandoItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(GeneName));
}

void PandoItem::__show_this() {

	MatrixWindow::show_matrix(
		&this->data()->mat_,
		this->title_,
		this->signal_emitter_,
		false,
		this->data_
	);
};

void PandoItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Abstract", s_abstract);
	ADD_MAIN_ACTION("Show Significant", s_show_significant);
	ADD_MAIN_ACTION("Show Transcriptional Factor", s_show_factor);
	ADD_MAIN_ACTION("Show Transcriptional Factor-Gene Coefficient", s_show_factor_coef);
	ADD_MAIN_ACTION("Show Gene", s_show_gene_name);
	ADD_MAIN_ACTION("Show Gene-Transcriptional Factor Coefficient", s_show_gene_coef);
	ADD_MAIN_ACTION("Show TF-Gene Strength", s_show_factor_gene_strength);

	ADD_MAIN_ACTION("Extract Gene Name", s_extract_gene_name);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void PandoItem::s_extract_gene_name() {

	G_GETLOCK;

	QStringList tfs = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME));

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Transcriptional Factor",
		{ "Transcriptional Factor Name" , "Correlation",  "P Threshold:0.05",  "N Variable:10" },
		{soap::InputStyle::LineEditWithCompleter, soap::InputStyle::ComboBox, 
		soap::InputStyle::NumericLineEdit, soap::InputStyle::IntegerLineEdit},
		{tfs, {"All", "Positive", "Negative"}}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString tf_name = settings[0];
	QString cor_type = settings[1];
	double p_val_threshold = settings[2].toDouble();
	int n_var_threshold = settings[3].toInt();

	auto filter = custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME), tf_name);
	filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	filter *= custom::greater_equal(this->data()->mat_.get_const_integer_reference(METADATA_PANDO_N_VARIABLE_NAME), n_var_threshold);

	if (cor_type == "Positive") {
		filter *= custom::greater_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_ESTIMATE_NAME), 0.0);
	}
	else if (cor_type == "Negative") {
		filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_ESTIMATE_NAME), 0.0);
	}

	if (filter.count() == 0) {
		G_NOTICE("No Gene found.");
		G_UNLOCK;
		return;
	}

	QStringList gene_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME), filter);

	QString new_title = this->signal_emitter_->get_unique_name(VARIABLE_GENE_NAME);
	DATA_SUBMODULES(GeneName)[new_title] = GeneName(gene_names);

	auto* item = new GeneNameItem(
		new_title,
		this->index_tree_,
		&DATA_SUBMODULES(GeneName)[new_title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_UNLOCK;
};

void PandoItem::s_show_gene_name() {

	QStringList gene_names = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME));

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Gene Name",
		{ "Gene Name" ,  "P Threshold:0.05" },
		{soap::InputStyle::LineEditWithCompleter, soap::InputStyle::NumericLineEdit},
		{gene_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString gene_name = settings[0];
	double p_val_threshold = settings[1].toDouble();

	auto filter = custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME), gene_name);
	filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	if (filter.count() == 0) {
		G_NOTICE("No significant differential feature expression found.");
		return;
	}

	auto tmp = this->data()->mat_.row_sliced(filter);
	MatrixWindow::show_matrix(
		&tmp,
		gene_name,
		this->signal_emitter_);
};

void PandoItem::s_show_gene_coef() {

	QStringList gene_names = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME));

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Gene",
		{ "Gene Name",  "P Threshold:0.05",  "N Variable:10" ,
		"Max Show Number:10" },
		{soap::InputStyle::LineEditWithCompleter,
		soap::InputStyle::NumericLineEdit, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit},
		{gene_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString gene_name = settings[0];
	double p_val_threshold = settings[1].toDouble();
	int n_variable_threshold = settings[2].toInt();
	int n_max_show = settings[3].toInt();

	if (n_max_show <= 0) {
		G_WARN("Show Number must be positive.");
		return;
	}

	auto filter = custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME), gene_name);
	filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	filter *= custom::greater_equal(this->data()->mat_.get_const_integer_reference(METADATA_PANDO_N_VARIABLE_NAME), n_variable_threshold);
	if (filter.count() == 0) {
		G_NOTICE("No significant differential feature expression found.");
		return;
	}

	QStringList tf_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME), filter);
	QVector<double> coef = custom::sliced(this->data()->mat_.get_const_double_reference(METADATA_PANDO_ESTIMATE_NAME), filter);

	QStringList unique_tfs = custom::unique(tf_names);
	QVector<double> coef_use;

	for (auto&& tf_name : unique_tfs) {
		coef_use << custom::sum(custom::sliced(coef, custom::equal(tf_names, tf_name)));
	}

	auto order = custom::order(coef_use);

	unique_tfs = custom::reordered(unique_tfs, order);
	coef_use = custom::reordered(coef_use, order);

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	custom_plot::patch::remove_left_axis(axis_rect);

	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);

	int n_tf = unique_tfs.size();
	int n_positive = std::ranges::count(custom::greater_than(coef_use, 0.0), true);
	int n_negative = n_tf - n_positive;

	filter = Eigen::ArrayX<bool>::Constant(n_tf, true);
	if (n_negative > n_max_show) {
		filter.segment(n_max_show, n_negative - n_max_show) = false;
		n_negative = n_max_show;
	}

	if (n_positive > n_max_show) {
		filter.segment(n_tf - n_positive, n_positive - n_max_show) = false;
		n_positive = n_max_show;
	}

	n_tf = n_positive + n_negative;

	if (n_tf < 1) {
		return;
	}

	if (custom::any(filter, false)) {
		coef_use = custom::sliced(coef_use, filter);
		unique_tfs = custom::sliced(unique_tfs, filter);
	}

	Eigen::ArrayXd bar_loc = Eigen::ArrayXd::LinSpaced(n_tf, 1, n_tf);

	custom_plot::patch::bar_polychrome(draw_area, axis_rect, bar_loc, custom::cast<Eigen::ArrayX>(coef_use),
		QList<QColor>() << QList<QColor>(n_negative, custom_plot::color::navy) << QList<QColor>(n_positive, custom_plot::color::firebrick3),
		24, false
	);

	for (int i = 0; i < n_tf; ++i) {
		if (coef_use[i] > 0) {
			custom_plot::patch::add_label(
				draw_area,
				axis_rect,
				unique_tfs[i] + ' ',
				0.0,
				i + 1.0,
				gs.get_left_label_font(),
				Qt::AlignVCenter | Qt::AlignRight
			);
		}
		else {
			custom_plot::patch::add_label(
				draw_area,
				axis_rect,
				' ' + unique_tfs[i],
				0.0,
				i + 1.0,
				gs.get_left_label_font(),
				Qt::AlignVCenter | Qt::AlignLeft
			);
		}
	}

	QPen pen(Qt::black);
	pen.setWidth(3);
	axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
	custom_plot::set_bottom_title(axis_rect, "Coefficient", gs, true);
	custom_plot::add_title(draw_area, gene_name, gs);

	auto [min_coef, max_coef] = std::ranges::minmax(coef_use);
	if (min_coef > -0.05) {
		min_coef = -0.05;
	}

	if (max_coef < 0.05) {
		max_coef = 0.05;
	}

	custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(min_coef, max_coef), QCPRange(0.0, n_tf + 1.0));

	this->draw_suite_->update(draw_area);
};

void PandoItem::s_show_factor_gene_strength() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This data is not a part of single-cell multiome data.")
	}

	SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

	QStringList embedding_names;

	for (auto&& [_, data_field] : SUBMODULES(*single_cell_multiome, DataField)) {
		for (auto&& [name, __] : SUBMODULES(data_field, Embedding)) {
			embedding_names << name;
		}
	}

	if (embedding_names.isEmpty()) {
		G_WARN("No Embedding for visualization.");
		return;
	}

	QStringList tfs = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME));

	QStringList gene_names = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME));

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Transcriptional Factor",
		{ "Transcriptional Factor Name",  "Gene Name", "Embedding", "Add to metadata:no"},
		{soap::InputStyle::LineEditWithCompleter, soap::InputStyle::LineEditWithCompleter,
		soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton},
		{tfs, gene_names, embedding_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString tf_name = settings[0];
	QString gene_name = settings[1];
	QString embedding_name = settings[2];
	bool add_to_metadata = switch_to_bool(settings[3]);

	auto filter = custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME), tf_name);
	filter *= custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME), gene_name);
	filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), 0.05);

	if (filter.count() == 0) {
		G_WARN("No Connection between TF and Gene");
		return;
	}

	auto peak_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_PEAK_NAME), filter);
	auto coef = custom::sliced(this->data()->mat_.get_const_double_reference(METADATA_PANDO_ESTIMATE_NAME), filter);
	int n_peak = peak_names.size();
	filter = Eigen::ArrayX<bool>::Constant(n_peak, true);

	auto rna_normalized = single_cell_multiome->rna_normalized();
	if (rna_normalized == nullptr) {
		G_WARN("No RNA Normalized Data.");
		return;
	}

	auto tf_index = rna_normalized->rownames_.indexOf(tf_name);
	if (tf_index == -1) {
		G_WARN("No Expression data of " + tf_name);
	}

	Eigen::ArrayXd tf_exp = rna_normalized->mat_.row(tf_index);

	auto atac_normalized = single_cell_multiome->atac_normalized();
	if (atac_normalized == nullptr) {
		G_WARN("No ATAC Normalized Data.");
		return;
	}

	Eigen::ArrayXd res = Eigen::ArrayXd::Zero(tf_exp.size());

	bool valid{ false };

	for (int i = 0; i < n_peak; ++i) {

		int peak_index = atac_normalized->rownames_.indexOf(peak_names[i]);
		if (peak_index == -1) {
			G_WARN("No ATAC data for " + peak_names[i]);
			continue;
		}

		Eigen::ArrayXd peak_acc = atac_normalized->mat_.row(peak_index);
		double estimate = coef[i];

		res += tf_exp * peak_acc * estimate;

		valid = true;
	}

	if (!valid) {
		G_WARN("No Peak data detected.");
		return;
	}

	if (add_to_metadata) {
		single_cell_multiome->metadata()->mat_.update(tf_name + " - " + gene_name + " Coefficient", custom::cast<QVector>(res));
		this->signal_emitter_->x_update_interface();
	}

	auto item = this->signal_emitter_->get_item(embedding_name);
	if (item == nullptr) {
		G_WARN("Embedding item not detected.");
	}

	auto embedding_item = static_cast<EmbeddingItem*>(item);

	embedding_item->show(res, tf_name + " - " + gene_name, false, "Coefficient");

};

void PandoItem::s_abstract() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Abstract Settings",
		{ "P Threshold:0.05",  "N Variable:10" },
		{soap::InputStyle::NumericLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (settings.isEmpty()) {
		return;
	}

	double p_val_threshold = settings[1].toDouble();
	int n_variable_threshold = settings[2].toInt();

	auto filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	filter *= custom::greater_equal(this->data()->mat_.get_const_integer_reference(METADATA_PANDO_N_VARIABLE_NAME), n_variable_threshold);
	if (filter.count() == 0) {
		G_NOTICE("No feature found.");
		return;
	}

	QStringList gene_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME), filter);
	QStringList tf_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME), filter);
	QVector<double> coef = custom::sliced(this->data()->mat_.get_const_double_reference(METADATA_PANDO_ESTIMATE_NAME), filter);

	std::unordered_map<QString, std::map<QString, double>> coefs;
	int n_conn = gene_names.size();
	for (int i = 0; i < n_conn; ++i) {
		coefs[tf_names[i]][gene_names[i]] += coef[i];
	}

	gene_names.clear();
	tf_names.clear();
	coef.clear();

	for (auto&& [tf_name, reg] : coefs) {
		for (auto&& [gene_name, score] : reg) {
			tf_names << tf_name;
			gene_names << gene_name;
			coef << score;
		}
	}

	auto order = custom::order(custom::abs(coef), true);

	CustomMatrix show;
	show.update("Transcription Factor", custom::reordered(tf_names, order));
	show.update("Gene Name", custom::reordered(gene_names, order));
	show.update("Coefficient", custom::reordered(coef, order));

	MatrixWindow::show_matrix(&show, "Abstract", this->signal_emitter_);
};

void PandoItem::s_show_factor_coef() {

	QStringList tfs = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME));

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Transcriptional Factor",
		{ "Transcriptional Factor Name",  "P Threshold:0.05",  "N Variable:10",
		"Max Show Number:10"},
		{soap::InputStyle::LineEditWithCompleter, 
		soap::InputStyle::NumericLineEdit, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit},
		{tfs}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString tf_name = settings[0];
	double p_val_threshold = settings[1].toDouble();
	int n_variable_threshold = settings[2].toInt();
	int n_max_show = settings[3].toInt();

	if (n_max_show <= 0) {
		G_WARN("Show Number must be positive");
		return;
	}

	auto filter = custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME), tf_name);
	filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	filter *= custom::greater_equal(this->data()->mat_.get_const_integer_reference(METADATA_PANDO_N_VARIABLE_NAME), n_variable_threshold);
	if (filter.count() == 0) {
		G_NOTICE("No significant differential feature expression found.");
		return;
	}

	QStringList gene_names = custom::sliced(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_GENE_NAME), filter);
	QVector<double> coef = custom::sliced(this->data()->mat_.get_const_double_reference(METADATA_PANDO_ESTIMATE_NAME), filter);

	QStringList unique_genes = custom::unique(gene_names);
	QVector<double> coef_use;

	for (auto&& gene_name : unique_genes) {
		coef_use << custom::sum(custom::sliced(coef, custom::equal(gene_names, gene_name)));
	}

	auto order = custom::order(coef_use);

	unique_genes = custom::reordered(unique_genes, order);
	coef_use = custom::reordered(coef_use, order);	

	int n_gene = unique_genes.size();
	int n_positive = std::ranges::count(custom::greater_than(coef_use, 0.0), true);
	int n_negative = n_gene - n_positive;

	filter = Eigen::ArrayX<bool>::Constant(n_gene, true);
	if (n_negative > n_max_show) {
		filter.segment(n_max_show, n_negative - n_max_show) = false;
		n_negative = n_max_show;
	}

	if (n_positive > n_max_show) {
		filter.segment(n_gene - n_positive, n_positive - n_max_show) = false;
		n_positive = n_max_show;
	}

	n_gene = n_positive + n_negative;

	if (n_gene < 1) {
		return;
	}

	if (custom::any(filter, false)) {
		coef_use = custom::sliced(coef_use, filter);
		unique_genes = custom::sliced(unique_genes, filter);
	}

	Eigen::ArrayXd bar_loc = Eigen::ArrayXd::LinSpaced(n_gene, 1, n_gene);
	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	custom_plot::patch::remove_left_axis(axis_rect);

	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);

	custom_plot::patch::bar_polychrome(draw_area, axis_rect, bar_loc, custom::cast<Eigen::ArrayX>(coef_use),
		QList<QColor>() << QList<QColor>(n_negative, custom_plot::color::navy) << QList<QColor>(n_positive, custom_plot::color::firebrick3),
		24, false
		);
		
	for (int i = 0; i < n_gene; ++i) {
		if (coef_use[i] > 0) {
			custom_plot::patch::add_label(
				draw_area,
				axis_rect,
				unique_genes[i] + ' ',
				0.0,
				i + 1.0,
				gs.get_left_label_font(),
				Qt::AlignVCenter | Qt::AlignRight
			);
		}
		else {
			custom_plot::patch::add_label(
				draw_area,
				axis_rect,
				' ' + unique_genes[i],
				0.0,
				i + 1.0,
				gs.get_left_label_font(),
				Qt::AlignVCenter | Qt::AlignLeft
			);
		}
	}

	QPen pen(Qt::black);
	pen.setWidth(3);
	axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
	custom_plot::set_bottom_title(axis_rect, "Coefficient", gs, true);
	custom_plot::add_title(draw_area, tf_name, gs);

	auto [min_coef, max_coef] = std::ranges::minmax(coef_use);
	if (min_coef > -0.05) {
		min_coef = -0.05;
	}

	if (max_coef < 0.05) {
		max_coef = 0.05;
	}

	custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(min_coef, max_coef), QCPRange(0.0, n_gene + 1.0));
	
	this->draw_suite_->update(draw_area);
};

void PandoItem::s_show_factor() {

	QStringList tfs = custom::unique(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME));

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Transcriptional Factor",
		{ "Transcriptional Factor Name",  "P Threshold:0.05", "N Variable:10" },
		{soap::InputStyle::LineEditWithCompleter, 
		soap::InputStyle::NumericLineEdit, soap::InputStyle::IntegerLineEdit},
		{tfs}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString tf_name = settings[0];
	double p_val_threshold = settings[1].toDouble();
	int n_variable_threshold = settings[2].toInt();

	auto filter = custom::equal(this->data()->mat_.get_const_qstring_reference(METADATA_PANDO_TF_NAME), tf_name);
	filter *= custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	filter *= custom::greater_equal(this->data()->mat_.get_const_integer_reference(METADATA_PANDO_N_VARIABLE_NAME), n_variable_threshold);
	if (filter.count() == 0) {
		G_NOTICE("No significant differential feature expression found.");
		return;
	}

	auto tmp = this->data()->mat_.row_sliced(filter);
	MatrixWindow::show_matrix(
		&tmp,
		tf_name,
		this->signal_emitter_);
};

void PandoItem::s_show_significant() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Significance Settings",
		{ "P Threshold:0.05", "N Variable:10"},
		{soap::InputStyle::NumericLineEdit, soap::InputStyle::IntegerLineEdit}
	);

	if (settings.isEmpty()) {
		return;
	}

	double p_val_threshold = settings[0].toDouble();
	int n_variable_threshold = settings[1].toInt();

	auto filter = custom::less_than(this->data()->mat_.get_const_double_reference(METADATA_PANDO_P_VAL_NAME), p_val_threshold);
	filter *= custom::greater_equal(this->data()->mat_.get_const_integer_reference(METADATA_PANDO_N_VARIABLE_NAME), n_variable_threshold);
	if (filter.count() == 0) {
		G_NOTICE("No significant differential feature expression found.");
		return;
	}

	auto tmp = this->data()->mat_.row_sliced(filter);
	MatrixWindow::show_matrix(
		&tmp,
		"Significant Result",
		this->signal_emitter_);
};

