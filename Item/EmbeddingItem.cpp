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
			{ handler.get_feature_names().completer_names}
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

		auto data = _Cs sapply(features, [&handler, normalized](auto&& t) {return handler.get_data({t, normalized}); });

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

		auto draw_area = _Cp feature_plot(data, this->data(), scale, nrow, this->draw_suite_->graph_settings_);
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
			{ handler.get_feature_names().completer_names, {"RNA", "ATAC", "Gene Activity"}}
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

		auto data = _Cs sapply(features, 
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

		auto draw_area = _Cp feature_plot(data, this->data(), scale, nrow, this->draw_suite_->graph_settings_);
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

		auto [draw_area, _, __] = _Cp embedding_single_color_plot(
			this->title_,
			this->data()->data_.mat_,
			Qt::blue,
			embedding_names,
			this->draw_suite_->graph_settings_
		);

		this->draw_suite_->update(draw_area);
	}
	else {

		auto [draw_area, _, __] = _Cp embedding_single_color_plot(
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
	d.dd = _Cs cast<QVector>(data);
	d.info["Source"] = "Others";

	auto [draw_area, _, __] = _Cp feature_plot(d, this->data(), scale, this->draw_suite_->graph_settings_);

	this->draw_suite_->update(draw_area);

};

void EmbeddingItem::show(
	const QStringList& data,
	const QString& name
) {

	QUERY_DATA d;
	d.name = name;
	d.type = QUERY_DATA::DataType::string;
	d.ds = _Cs cast<QVector>(data);
	d.dsl = _Cs unique(d.ds);
	d.info["Source"] = "Others";

	auto [draw_area, _, __] = _Cp feature_plot(d, this->data(), false, this->draw_suite_->graph_settings_);

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

	QVector<double> x = _Cs linspaced(y.size(), 1, y.size());
	auto& gs = this->draw_suite_->graph_settings_;
	
	auto [draw_area, axis_rect] = _Cp prepare_ar(gs);
	
	_Cp set_simple_axis(axis_rect, "PCA components", "Standard Deviation", gs);

	_CpPatch scatter(draw_area, axis_rect, x, y, Qt::black, 5);

	_CpPatch set_range(axis_rect, QCPRange(0, y.size() + 1), QCPRange(0, y[0] + 1));
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
		double cor = _Cs correlation_pearson(depth, mat.col(i));
		if (std::isnan(cor)) {
			cor = 0.0;
		}
		y[i] = cor;
	}

	auto& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect] = _Cp prepare_ar(gs);
	
	_Cp set_simple_axis(axis_rect, "PCA components", "Correlation with Sequencing Depth", gs);

	_CpPatch scatter(draw_area, axis_rect, _Cs linspaced(n_dimension, 1, n_dimension), y, Qt::black, 5);

	_CpPatch set_range(axis_rect, QCPRange(0, n_dimension + 1), _CpUtility get_range(y));

	this->draw_suite_->update(draw_area);
};