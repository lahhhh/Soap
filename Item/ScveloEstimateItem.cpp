#include "ScveloEstimateItem.h"

#include "MatrixWindow.h"
#include "EnrichWorker.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CustomPlot.h"
#include "FileWritingWorker.h"

#include "ShowEmbeddingScveloWorker.h"

#include "StreamPlot.h"
#include "FeatureHandler.h"

void ScveloEstimateItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Show Embedding Velocity");
		
	ADD_ACTION("by grid", "Show Embedding Velocity", s_show_embedding_velocity_by_grid);
	ADD_ACTION("by stream", "Show Embedding Velocity", s_show_embedding_velocity_by_stream);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
};

void ScveloEstimateItem::s_show_embedding_velocity_by_stream() {
	this->show_embedding_velocity(1);
}

void ScveloEstimateItem::s_show_embedding_velocity_by_grid() {
	this->show_embedding_velocity(0);
}

void ScveloEstimateItem::show_embedding_velocity(int mode) {

	QStringList embedding_names = this->index_tree_->search(soap::VariableType::Embedding);

	if (embedding_names.isEmpty()) {
		G_NOTICE("No embedding for velocity analysis.");
		return;
	}

	G_GETLOCK;

	QStringList valid_features, settings;

	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		SingleCellRna* single_cell_rna = this->get_root<SingleCellRna>();
		valid_features << single_cell_rna->metadata()->mat_.colnames_ << single_cell_rna->counts()->rownames_;

		settings = CommonDialog::get_response(
			this->signal_emitter_, 
			"Scvelo Analysis settings",
			{ "Embedding", "Feature", "Normalized" },
			{ soap::InputStyle::ComboBox, soap::InputStyle::LineEditWithCompleter, soap::InputStyle::SwitchButton},
			{ embedding_names, valid_features}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		if (!valid_features.contains(settings[1])) {
			G_WARN("Feature is not valid.");
			G_UNLOCK;
			return;
		}

	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();
		valid_features << single_cell_multiome->metadata()->mat_.colnames_ << single_cell_multiome->rna_counts()->rownames_;

		settings = CommonDialog::get_response(
			this->signal_emitter_, 
			"Scvelo Analysis settings",
			{ "Embedding", "Feature", "Normalized", "Field"},
			{ soap::InputStyle::ComboBox, soap::InputStyle::LineEditWithCompleter
			, soap::InputStyle::SwitchButton, soap::InputStyle::ComboBox},
			{ embedding_names, valid_features, { "RNA", "ATAC", "Gene Activity" }}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		if (!valid_features.contains(settings[1])) {
			G_WARN("Feature is not valid.");
			G_UNLOCK;
			return;
		}
	}
	else {
		G_WARN("This velocity estimate is not linked with any single cell rna or multiome data.");
		G_UNLOCK;
		return;
	}


	auto tree = this->index_tree_->search(settings[0]);

	if (tree == nullptr) {
		G_UNLOCK;
		return;
	}
	Embedding* embedding = static_cast<Embedding*>(tree->data_);

	ShowEmbeddingScveloWorker* worker = new ShowEmbeddingScveloWorker(this->data(), *embedding, settings, mode);

	if (mode == 0) {
		G_LINK_WORKER_THREAD(ShowEmbeddingScveloWorker, x_grid_graph_ready, ScveloEstimateItem, s_receive_velocity_grid);
	}
	else {
		G_LINK_WORKER_THREAD(ShowEmbeddingScveloWorker, x_stream_graph_ready, ScveloEstimateItem, s_receive_velocity_stream);
	}
};

std::pair<QCustomPlot*, QCPAxisRect*> ScveloEstimateItem::draw_feature_plot(
	const Eigen::MatrixXd& embedding, 
	const QStringList& embedding_names, 
	const QStringList& graph_settings
) {
	QString feature = graph_settings[1];
	bool normalize = switch_to_bool(graph_settings[2]);
	bool gene_activity{ false };
	QString legend_title = gene_activity ? "Activity" : "Expression";

	FeatureHandler handler;

	if (this->stem_from(soap::VariableType::SingleCellRna)) {

		if (graph_settings.size() != 3) {
			G_WARN("Unmatched graph settings!");
			return { nullptr, nullptr };
		}

		SingleCellRna* single_cell_rna = this->get_root<SingleCellRna>();
		handler.set(single_cell_rna);

	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {

		if (graph_settings.size() != 4) {
			G_WARN("Unmatched graph settings!");
			return { nullptr, nullptr };
		}

		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();
		QString field = graph_settings[3];
		gene_activity = field == "Gene Activity";

		handler.set(single_cell_multiome);
	}
	else {
		return { nullptr, nullptr };
	}

	auto feature_data = handler.get_data({ feature, normalize, gene_activity });
	if (!feature_data.is_valid()) {
		G_WARN("Feature Not Found");
		return { nullptr, nullptr };
	}

	if (feature_data.info["Source"] == "ATAC") {
		legend_title = "Activity";
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = _Cp prepare(gs);	

	Eigen::ArrayXd x = embedding.col(0), y = embedding.col(1);

	_Cp set_scatter_plot_axis_style(draw_area, axis_rect, embedding_names[0], embedding_names[1], x, y, gs);

	if (feature_data.is_continuous()) {

		Eigen::ArrayXd feat = _Cs cast<Eigen::ArrayX>(feature_data.get_continuous());

		QColor low_color = gs.get_gradient_low_color();
		low_color.setAlpha(64);
		QColor middle_color = gs.get_gradient_middle_color();
		middle_color.setAlpha(64);
		QColor high_color = gs.get_gradient_high_color();
		high_color.setAlpha(64);

		_CpPatch scatter_gradient(
			draw_area,
			axis_rect,
			x,
			y,
			feat,
			feat.minCoeff(),
			feat.maxCoeff(),
			low_color,
			middle_color,
			high_color,
			gs.get_scatter_point_size()
		);

		_Cp add_gradient_legend(draw_area, legend_layout, feat.minCoeff(), feat.maxCoeff(), legend_title, gs);
	}
	else if(feature_data.is_factor()){

		auto levels = feature_data.get_levels();
		auto colors = gs.palette(levels);
		std::ranges::for_each(colors, [](QColor& color) {color.setAlpha(64); });

		_CpPatch scatter_category(draw_area, axis_rect, x, y, feature_data.get_factor(), levels, colors, gs.get_scatter_point_size());

		_CpPatch add_round_legend(
			draw_area, 
			legend_layout, 
			levels, 
			colors, 
			feature, 
			gs.get_legend_column_width(),
			gs.get_legend_row_width(),
			gs.get_legend_title_font(),
			gs.get_legend_label_font());
	}
	else {
		return { nullptr, nullptr };
	}

	_Cp add_title(draw_area, "Scvelo Graph", gs);

	return std::make_pair(draw_area, axis_rect);
};

void ScveloEstimateItem::s_receive_velocity_stream(STREAM_PLOT_ELEMENTS res) {
	auto [draw_area, axis_rect] = this->draw_feature_plot(res.embedding, res.embedding_names, res.graph_settings);

	if (draw_area == nullptr) {
		return;
	}

	stream_plot(res.x, res.y, res.u, res.v, res.mask, draw_area, axis_rect);

	this->draw_suite_->update(draw_area);
};

void ScveloEstimateItem::s_receive_velocity_grid(VELO_GRID_PLOT_ELEMENTS res) {

	auto [draw_area, axis_rect] = this->draw_feature_plot(res.embedding, res.embedding_names, res.graph_settings);

	if (draw_area == nullptr) {
		return;
	}

	double max_x = axis_rect->axis(QCPAxis::atBottom)->range().upper, min_x = axis_rect->axis(QCPAxis::atBottom)->range().lower;
	double max_y = axis_rect->axis(QCPAxis::atLeft)->range().upper, min_y = axis_rect->axis(QCPAxis::atLeft)->range().lower;

	double graph_y_x_ratio = (max_y - min_y) / (max_x - min_x);
	double sin = 1 / std::sqrt(1 + graph_y_x_ratio * graph_y_x_ratio), cos = graph_y_x_ratio / std::sqrt(1 + graph_y_x_ratio * graph_y_x_ratio);

	QPen arrow_pen;

	arrow_pen.setColor(Qt::black);
	arrow_pen.setWidth(1);

	const int n_arrow = res.arrows_start.rows();
	for (int i = 0; i < n_arrow; ++i) {
		double start_x = res.arrows_start(i, 0), start_y = res.arrows_start(i, 1);
		double end_x = start_x + res.direction(i, 0), end_y = start_y + res.direction(i, 1);

		double middle_x = 0.8 * (end_x - start_x) + start_x;
		double middle_y = 0.8 * (end_y - start_y) + start_y;
		double arrow_x1 = 0.66667 * (end_x - start_x) + start_x + 0.2 * (end_y - start_y) * sin;
		double arrow_y1 = 0.66667 * (end_y - start_y) + start_y - 0.2 * (end_x - start_x) * cos;
		double arrow_x2 = 0.66667 * (end_x - start_x) + start_x - 0.2 * (end_y - start_y) * sin;
		double arrow_y2 = 0.66667 * (end_y - start_y) + start_y + 0.2 * (end_x - start_x) * cos;

		draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		QVector<QCPCurveData> data(5);

		data[0] = QCPCurveData(0, end_x, end_y);
		data[1] = QCPCurveData(1, arrow_x1, arrow_y1);
		data[2] = QCPCurveData(2, middle_x, middle_y);
		data[3] = QCPCurveData(3, arrow_x2, arrow_y2);
		data[4] = QCPCurveData(4, end_x, end_y);

		shape->setPen(arrow_pen);
		shape->setBrush(QBrush(Qt::black));
		shape->data()->set(data, true);

		QCPGraph* graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		graph->setData(QVector<double>({ start_x, end_x }), QVector<double>({ start_y, end_y }));
		graph->setPen(arrow_pen);

	}

	this->draw_suite_->update(draw_area);

}


