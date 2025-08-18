#include "EmbeddingItem.h"

#include "Custom.h"
#include "CommonDialog.h"
#include "CustomPlot.h"
#include "EmbeddingViewWindow.h"
#include "ItemIOWorker.h"
#include "FileWritingWorker.h"
#include "MatrixWindow.h"

void EmbeddingItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->data_.mat_.rows()) + " | " + QString::number(this->data()->data_.mat_.cols()) + " ]");
};

void EmbeddingItem::__show_this() {

	MatrixWindow::show_matrix(this->data(), this->title_, this->signal_emitter_);
}

void EmbeddingItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Show", s_raw_plot);

	ADD_MAIN_ACTION("Feature Plot", s_feature_plot);

	ADD_MAIN_ACTION("Highlight", s_highlight);

	ADD_MAIN_ACTION("Numeric Feature Plot", s_multiple_numeric_feature_plot);

	ADD_MAIN_ACTION("View", s_view);

	if (this->data()->data_type_ == Embedding::DataType::Pca) {
		ADD_MAIN_ACTION("Deviation Plot", s_deviation_plot);
	}

	if (this->attached_to(soap::VariableType::DataField) && this->data()->data_type_ == Embedding::DataType::Pca) {

		DataField* data_field = this->trace_back<DataField>(1);
		if (data_field->data_type_ == DataField::DataType::Atac) {
			ADD_MAIN_ACTION("Show Depth Correlation", s_show_depth_correlation);
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac) && this->data()->data_type_ == Embedding::DataType::Pca) {

		ADD_MAIN_ACTION("Show Depth Correlation", s_show_depth_correlation);
	}

	if (this->attached_to(soap::VariableType::ChromVAR)) {
		ADD_MAIN_ACTION("Feature Plot (ChromVAR)", s_chromvar_plot);
	}

	if (this->attached_to(soap::VariableType::Cicero)) {
		ADD_MAIN_ACTION("Feature Plot (Cicero)", s_cicero_plot);
	}

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void EmbeddingItem::s_view() {

	if (this->stem_from(soap::VariableType::SingleCellRna)) {

		auto* single_cell_rna = this->get_root<SingleCellRna>();

		auto species = single_cell_rna->species_;
		if (species != soap::Species::Human && species != soap::Species::Mouse) {
			G_WARN("Species should be Human or Mouse for Embedding View.");
			return;
		}
		EmbeddingViewWindow::view(this->trace_back<SingleCellRna>(1), this->data(), this->signal_emitter_);
	}
	else if (this->stem_from(soap::VariableType::SingleCellAtac)) {

		auto* single_cell_atac = this->get_root<SingleCellAtac>();

		auto species = single_cell_atac->species_;
		if (species != soap::Species::Human && species != soap::Species::Mouse) {
			G_WARN("Species should be Human or Mouse for Embedding View.");
			return;
		}

		EmbeddingViewWindow::view(single_cell_atac, this->data(), this->signal_emitter_);
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {

		auto* single_cell_multiome = this->get_root<SingleCellMultiome>();

		auto species = single_cell_multiome->species_;
		if (species != soap::Species::Human && species != soap::Species::Mouse) {
			G_WARN("Species should be Human or Mouse for Embedding View.");
			return;
		}

		EmbeddingViewWindow::view(single_cell_multiome, this->data(), this->signal_emitter_);
	}
	else {
		G_WARN("Embedding is not attached to any data.");
		return;
	}
};

void EmbeddingItem::s_cicero_plot() {
	if (!this->attached_to(soap::VariableType::Cicero)) {
		G_WARN("Illegal operation.");
		return;
	}

	auto* cicero = this->trace_back<Cicero>(1);
	auto valid_features = cicero->regulation_group_normalized_.rownames_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Settings",
		{ "feature" },
		{soap::InputStyle::LineEditWithCompleter},
		{valid_features}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature = settings[0];

	int index = valid_features.indexOf(feature);

	if (index == -1) {
		G_WARN("Invalid feature.");
		return;
	}

	this->show(cicero->regulation_group_normalized_.mat_.row(index), feature, false, "Accessibility");
};

void EmbeddingItem::s_chromvar_plot() {

	if (!this->attached_to(soap::VariableType::ChromVAR)) {
		G_WARN("Illegal operation.");
		return;
	}

	auto* chrom_var = this->trace_back<ChromVAR>(1);

	auto valid_features = chrom_var->motif_names_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Settings",
		{ "feature(z value)" },
		{soap::InputStyle::LineEditWithCompleter},
		{valid_features}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature = settings[0];

	int index = valid_features.indexOf(feature);

	if (index == -1) {
		G_WARN("Invalid feature.");
		return;
	}

	this->show(chrom_var->z_.row(index), feature, false, "Z Score");
};

void EmbeddingItem::s_multiple_numeric_feature_plot() {

	Eigen::MatrixXd feature_mat;
	Eigen::ArrayXd x = this->data()->data_.mat_.col(0), y = this->data()->data_.mat_.col(1);
	int n_cell = x.size();
	bool scale{ false };
	int nrow{ 1 };
	QStringList valid_features;

	if (this->attached_to(soap::VariableType::SingleCellRna) ||
		this->attached_to(soap::VariableType::SingleCellAtac) ||
		this->attached_to(soap::VariableType::BulkRna)) {

		FeatureHandler handler;

		if (this->attached_to(soap::VariableType::SingleCellRna)) {
			SingleCellRna* single_cell_rna = this->trace_back<SingleCellRna>(1);
			handler.set(single_cell_rna);
		}
		else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
			SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);
			handler.set(single_cell_atac);
		}
		else if (this->attached_to(soap::VariableType::BulkRna)) {
			BulkRna* rna = this->trace_back<BulkRna>(1);
			handler.set(rna);
		}

		QStringList settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Feature Plot Setting",
			{ "Feature", "Normalized:yes", "Scale:no", "Number of row:1" },
			{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::SwitchButton,
			soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit },
			{ handler.get_feature_names().numeric_names }
		);

		if (settings.isEmpty()) {
			return;
		}

		QStringList features = multiple_line_edit_with_completer_to_list(settings[0]);
		if (features.isEmpty()) {
			return;
		}

		bool normalized = switch_to_bool(settings[1]);
		scale = switch_to_bool(settings[2]);
		nrow = settings[3].toInt();
		if (nrow < 0) {
			G_WARN("Number of rows can not be less than 1.");
			return;
		}

		auto data = custom::sapply(features, [&handler, normalized](auto&& t) {return handler.get_data({ t, normalized }); });

		int n_valid{ 0 };
		for (auto&& d : data) {
			if (!d.is_continuous()) {
				G_WARN("Feature " + d.name + " is not valid.");
				continue;
			}

			++n_valid;
		}
		if (n_valid == 0) {
			G_WARN("No Valid Feature.");
			return;
		}

		feature_mat.resize(n_cell, n_valid);
		int count{ 0 };
		for (auto&& d : data) {
			if (d.is_continuous()) {
				feature_mat.col(count++) = custom::cast<Eigen::ArrayX>(d.get_continuous());
				valid_features << d.name;
			}
		}

	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

		FeatureHandler handler(single_cell_multiome);

		QStringList settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Feature Plot Setting",
			{ "Feature", "Normalized:yes", "Scale:no", "Field", "Number of row:1" },
			{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton,
			soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit },
			{ handler.get_feature_names().numeric_names, {"RNA", "ATAC", "Gene Activity"} }
		);
		if (settings.isEmpty())return;

		QStringList features = multiple_line_edit_with_completer_to_list(settings[0]);

		if (features.isEmpty()) {
			return;
		}

		bool normalized = switch_to_bool(settings[1]);
		scale = switch_to_bool(settings[2]);

		bool gene_activity = settings[3] == "Gene Activity";

		nrow = settings[4].toInt();
		if (nrow < 0) {
			G_WARN("Number of rows can not be less than 1.");
			return;
		}

		auto data = custom::sapply(features,
			[&handler, normalized, gene_activity](auto&& t) {return handler.get_data({ t, normalized, gene_activity }); });

		int n_valid{ 0 };
		for (auto&& d : data) {
			if (!d.is_continuous()) {
				G_WARN("Feature " + d.name + " is not found.");
				continue;
			}

			++n_valid;
		}

		if (n_valid == 0) {
			G_WARN("No Valid Feature.");
			return;
		}

		feature_mat.resize(n_cell, n_valid);
		int count{ 0 };
		for (auto&& d : data) {
			if (d.is_continuous()) {
				feature_mat.col(count++) = custom::cast<Eigen::ArrayX>(d.get_continuous());
				valid_features << d.name;
			}
		}
	}
	else {
		return;
	}

	double min_val = feature_mat.minCoeff(), max_val = feature_mat.maxCoeff();
	int n_valid_feature = feature_mat.cols();

	if (scale) {
		for (int i = 0; i < n_valid_feature; ++i) {
			double mean = feature_mat.col(i).mean();

			feature_mat.col(i).array() -= mean;

			double sd = custom::sd(feature_mat.col(i));
			if (sd != 0.0) {
				feature_mat.col(i).array() /= sd;
			}
		}

		min_val = -1.0;
		max_val = 1.0;
	}

	auto& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, left_layout, legend_layout] = custom_plot::prepare_lg_lg(gs);

	for (int i = 0; i < n_valid_feature; ++i) {

		QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
		int col = i / nrow;
		int row = i - col * nrow;
		left_layout->addElement(row, col, sub_layout);

		QCPAxisRect* axis_rect = new QCPAxisRect(draw_area);
		sub_layout->addElement(0, 0, axis_rect);
		custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(x, y));
		custom_plot::patch::set_border_only(axis_rect, Qt::black, 2);
		custom_plot::patch::scatter_gradient(
			draw_area,
			axis_rect,
			x,
			y,
			feature_mat.col(i),
			min_val,
			max_val,
			gs.get_gradient_low_color(),
			gs.get_gradient_middle_color(),
			gs.get_gradient_high_color(),
			gs.get_scatter_point_size()
		);

		custom_plot::patch::add_title(draw_area, sub_layout, valid_features[i], gs.get_title_font());
	}

	if (scale) {
		custom_plot::add_gradient_legend(
			draw_area,
			legend_layout,
			-1.0,
			1.0,
			"Expression",
			gs,
			"Low",
			"High"
		);
	}
	else {

		custom_plot::add_gradient_legend(draw_area, legend_layout, min_val, max_val, "Expression", gs);
	}

	this->draw_suite_->update(draw_area);
};

void EmbeddingItem::s_highlight() {

	void* d = this->index_tree_->search_one(soap::VariableType::Metadata);

	if (d == nullptr) {
		G_WARN("No Metadata Found.");
		return;
	}

	Metadata* metadata = static_cast<Metadata*>(d);
	auto factor_info = metadata->mat_.get_factor_information(false);

	if (factor_info.isEmpty()) {
		G_WARN("No Suitable Metadata for Visualisation.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Highlight Settings",
		{ "Group", "Background Color:#ececec" },
		{ soap::InputStyle::FactorChoice, soap::InputStyle::ColorChoice },
		{},
		{ factor_info }
	);

	if (settings.isEmpty()) {
		return;
	}

	auto [factor_name, levels] = factor_choice_to_pair(settings[0]);
	if (levels.isEmpty()) {
		return;
	}

	QColor bg_color = QColor::fromString(settings[1]);

	auto* embedding = this->data();

	Eigen::ArrayXd x = embedding->data_.mat_.col(0);
	Eigen::ArrayXd y = embedding->data_.mat_.col(1);
	QString bottom_title = embedding->data_.colnames_[0];
	QString left_title = embedding->data_.colnames_[1];

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(x, y));
	custom_plot::set_scatter_plot_axis_style(draw_area, axis_rect, bottom_title, left_title, x, y, gs);

	const int n_level = levels.size();
	QStringList values = metadata->mat_.get_qstring(factor_name);
	auto colors = gs.palette(levels);

	Eigen::ArrayX<bool> bg_filter = Eigen::ArrayX<bool>::Constant(values.size(), false);

	for (int i = 0; i < n_level; ++i) {

		auto filter = custom::equal(values, levels[i]);
		bg_filter += filter;
	}

	bg_filter = custom::flip(bg_filter);
	if (bg_filter.count() > 0) {

		auto sub_x = custom::sliced(x, bg_filter);
		auto sub_y = custom::sliced(y, bg_filter);

		custom_plot::patch::scatter(
			draw_area,
			axis_rect,
			sub_x,
			sub_y,
			bg_color,
			gs.get_scatter_point_size()
		);
	}

	for (int i = 0; i < n_level; ++i) {

		auto filter = custom::equal(values, levels[i]);
		auto sub_x = custom::sliced(x, filter);
		auto sub_y = custom::sliced(y, filter);

		custom_plot::patch::scatter(
			draw_area,
			axis_rect,
			sub_x,
			sub_y,
			colors[i],
			gs.get_scatter_point_size()
		);
	}

	custom_plot::patch::add_round_legend(
		draw_area,
		legend_layout,
		levels,
		colors,
		gs.get_legend_title(""),
		gs.get_legend_column_width(),
		gs.get_legend_row_width(),
		gs.get_legend_title_font(),
		gs.get_legend_label_font()
	);

	custom_plot::add_title(draw_area, factor_name, gs);

	this->draw_suite_->update(draw_area);
};

void EmbeddingItem::s_feature_plot() {

	if (this->attached_to(soap::VariableType::SingleCellRna) ||
		this->attached_to(soap::VariableType::SingleCellAtac) ||
		this->attached_to(soap::VariableType::BulkRna)) {

		FeatureHandler handler;

		if (this->attached_to(soap::VariableType::SingleCellRna)) {
			SingleCellRna* single_cell_rna = this->trace_back<SingleCellRna>(1);
			handler.set(single_cell_rna);
		}
		else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
			SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);
			handler.set(single_cell_atac);
		}
		else if (this->attached_to(soap::VariableType::BulkRna)) {
			BulkRna* rna = this->trace_back<BulkRna>(1);
			handler.set(rna);
		}

		QStringList settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Feature Plot Setting",
			{ "Feature", "Normalized:yes", "Scale:no", "Number of row:1" },
			{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::SwitchButton, 
			soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit },
			{ handler.get_feature_names().numeric_names}
		);

		if (settings.isEmpty()) {
			return;
		}

		QStringList features = multiple_line_edit_with_completer_to_list(settings[0]);
		if (features.isEmpty()) {
			return;
		}

		bool normalized = switch_to_bool(settings[1]);
		bool scale = switch_to_bool(settings[2]);
		int nrow = settings[3].toInt();
		if (nrow < 0) {
			G_WARN("Number of rows can not be less than 1.");
			return;
		}

		auto data = custom::sapply(features, [&handler, normalized](auto&& t) {return handler.get_data({t, normalized}); });

		int n_valid{ 0 };
		for (auto&& d : data) {
			if (!d.is_valid()) {
				G_WARN("Feature " + d.name + " is not found.");
				continue;
			}

			++n_valid;
		}
		if (n_valid == 0) {
			G_WARN("No Valid Feature.");
			return;
		}

		auto draw_area = custom_plot::feature_plot(data, this->data(), scale, nrow, this->draw_suite_->graph_settings_);
		this->draw_suite_->update(draw_area);
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

		FeatureHandler handler(single_cell_multiome);

		QStringList settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Feature Plot Setting",
			{ "Feature", "Normalized:yes", "Scale:no", "Field", "Number of row:1"},
			{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton,
			soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit},
			{ handler.get_feature_names().numeric_names, {"RNA", "ATAC", "Gene Activity"}}
		);
		if (settings.isEmpty())return;

		QStringList features = multiple_line_edit_with_completer_to_list(settings[0]);

		if (features.isEmpty()) {
			return;
		}

		bool normalized = switch_to_bool(settings[1]);
		bool scale = switch_to_bool(settings[2]);

		bool gene_activity = settings[3] == "Gene Activity";

		int nrow = settings[4].toInt();
		if (nrow < 0) {
			G_WARN("Number of rows can not be less than 1.");
			return;
		}

		auto data = custom::sapply(features, 
			[&handler, normalized, gene_activity](auto&& t) {return handler.get_data({ t, normalized, gene_activity }); });

		int n_valid{ 0 };
		for (auto&& d : data) {
			if (!d.is_valid()) {
				G_WARN("Feature " + d.name + " is not found.");
				continue;
			}

			++n_valid;
		}
		if (n_valid == 0) {
			G_WARN("No Valid Feature.");
			return;
		}

		auto draw_area = custom_plot::feature_plot(data, this->data(), scale, nrow, this->draw_suite_->graph_settings_);
		this->draw_suite_->update(draw_area);
	}
};

void EmbeddingItem::s_raw_plot() {

	if (this->attached_to(soap::VariableType::BulkRna) && this->data()->data_type_ == Embedding::DataType::Pca) {

		QStringList embedding_names = this->data()->data_.colnames_;

		BulkRna* rna = this->trace_back<BulkRna>(1);
		if (rna->double_vectors_.contains(VARIABLE_PCA_VARIANCE_PROPORTION)) {
			auto&& vp = rna->double_vectors_[VARIABLE_PCA_VARIANCE_PROPORTION];
			if (vp.size() >= 2) {
				embedding_names[0] = embedding_names[0] + " ( " + QString::number(vp[0] * 100, 'f', 2) + "% Variation )";
				embedding_names[1] = embedding_names[1] + " ( " + QString::number(vp[1] * 100, 'f', 2) + "% Variation )";
			}
		}

		auto [draw_area, _, __] = custom_plot::embedding_single_color_plot(
			this->title_,
			this->data()->data_.mat_,
			Qt::blue,
			embedding_names,
			this->draw_suite_->graph_settings_
		);

		this->draw_suite_->update(draw_area);
	}
	else {

		auto [draw_area, _, __] = custom_plot::embedding_single_color_plot(
			this->title_,
			this->data()->data_.mat_,
			Qt::blue,
			this->data()->data_.colnames_,
			this->draw_suite_->graph_settings_
		);

		this->draw_suite_->update(draw_area);
	}
}

void EmbeddingItem::show(
	const Eigen::ArrayXd& data,
	const QString& title,
	bool scale,
	const QString& legend_title
) {
	QUERY_DATA d;
	d.name = title;
	d.type = QUERY_DATA::DataType::numeric;
	d.dd = custom::cast<QVector>(data);
	d.info["Source"] = "Others";

	auto [draw_area, _, __] = custom_plot::feature_plot(d, this->data(), scale, this->draw_suite_->graph_settings_);

	this->draw_suite_->update(draw_area);

};

void EmbeddingItem::show(
	const QStringList& data,
	const QString& name
) {

	QUERY_DATA d;
	d.name = name;
	d.type = QUERY_DATA::DataType::string;
	d.ds = custom::cast<QVector>(data);
	d.dsl = custom::unique(d.ds);
	d.info["Source"] = "Others";

	auto [draw_area, _, __] = custom_plot::feature_plot(d, this->data(), false, this->draw_suite_->graph_settings_);

	this->draw_suite_->update(draw_area);

};

void EmbeddingItem::s_deviation_plot() {

	if (this->data()->data_type_ != Embedding::DataType::Pca) {
		return;
	}

	QVector<double> y;

	if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		if (this->attached_to(soap::VariableType::DataField)) {
			SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();
			DataField* data_field = this->trace_back<DataField>(1);

			QString standard_deviation_name;
			auto type = data_field->data_type_;
			switch (type)
			{
			case DataField::DataType::Rna:
				standard_deviation_name = VARIABLE_RNA_PCA_STANDARD_DEVIATION;
				break;
			case DataField::DataType::Atac:
				standard_deviation_name = VARIABLE_ATAC_PCA_STANDARD_DEVIATION;
				break;
			case DataField::DataType::Trans:
				standard_deviation_name = VARIABLE_TRANS_PCA_STANDARD_DEVIATION;
				break;
			default:
				break;
			}

			if (!single_cell_multiome->double_vectors_.contains(standard_deviation_name)) {
				G_NOTICE("No standard deviation data.");
				return;
			}

			y = single_cell_multiome->double_vectors_[standard_deviation_name];
		}
		else {
			G_WARN("This embedding is not attached to any data field.");
			return;
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellRna)) {

		SingleCellRna* single_cell_rna = this->trace_back<SingleCellRna>(1);

		if (!single_cell_rna->double_vectors_.contains(VARIABLE_RNA_PCA_STANDARD_DEVIATION)) {
			G_WARN("No standard deviation data found.");
			return;
		}
		y = single_cell_rna->double_vectors_[VARIABLE_RNA_PCA_STANDARD_DEVIATION];
	}
	else if (this->attached_to(soap::VariableType::BulkRna)) {
		BulkRna* rna = this->trace_back<BulkRna>(1);
		if (!rna->double_vectors_.contains(VARIABLE_PCA_STANDARD_DEVIATION)) {
			G_WARN("No standard deviation data found.");
			return;
		}
		y = rna->double_vectors_[VARIABLE_PCA_STANDARD_DEVIATION];
	}
	else {
		G_WARN("No standard deviation data found.");
		return;
	}

	if (y.isEmpty()) {
		G_WARN("Invalid standard deviation data.");
		return;
	}

	QVector<double> x = custom::linspaced(y.size(), 1, y.size());
	auto& gs = this->draw_suite_->graph_settings_;
	
	auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);
	
	custom_plot::set_simple_axis(axis_rect, "PCA components", "Standard Deviation", gs);

	custom_plot::patch::scatter(draw_area, axis_rect, x, y, Qt::black, 5);

	custom_plot::patch::set_range(axis_rect, QCPRange(0, y.size() + 1), QCPRange(0, y[0] + 1));
	this->draw_suite_->update(draw_area);
};


void EmbeddingItem::s_show_depth_correlation() {

	auto& mat = this->data()->data_.mat_;
	int n_dimension = mat.cols();

	QVector<double> depth;

	if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

		depth = single_cell_multiome->metadata()->mat_.get_double(METADATA_ATAC_UMI_NUMBER);
		if (depth.isEmpty() || 
			single_cell_multiome->metadata()->mat_.data_type_[METADATA_ATAC_UMI_NUMBER] != CustomMatrix::DataType::IntegerNumeric) {
			G_NOTICE("Depth data has lost.");
			return;
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
		SingleCellAtac* single_cell_atac = this->get_root<SingleCellAtac>();

		depth = single_cell_atac->metadata()->mat_.get_double(METADATA_ATAC_UMI_NUMBER);
		if (depth.isEmpty() ||
			single_cell_atac->metadata()->mat_.data_type_[METADATA_ATAC_UMI_NUMBER] != CustomMatrix::DataType::IntegerNumeric) {
			G_NOTICE("Depth data has lost.");
			return;
		}
	}
	else {
		G_WARN("Illegal Operation");
		return;
	}

	QVector<double> y(n_dimension);

	for (int i = 0; i < n_dimension; ++i) {
		double cor = custom::correlation_pearson(depth, mat.col(i));
		if (std::isnan(cor)) {
			cor = 0.0;
		}
		y[i] = cor;
	}

	auto& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);
	
	custom_plot::set_simple_axis(axis_rect, "PCA components", "Correlation with Sequencing Depth", gs);

	custom_plot::patch::scatter(draw_area, axis_rect, custom::linspaced(n_dimension, 1, n_dimension), y, Qt::black, 5);

	custom_plot::patch::set_range(axis_rect, QCPRange(0, n_dimension + 1), custom_plot::utility::get_range(y));

	this->draw_suite_->update(draw_area);
};