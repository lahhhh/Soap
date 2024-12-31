#include "CellChatItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"

#include "CustomPlot.h"
#include "CircosPlot.h"

void CellChatItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_MENU("Statistics");

	ADD_ACTION("Number of interactions", "Statistics", s_show_interaction_numbers);
	ADD_ACTION("Weights of interactions", "Statistics", s_show_interaction_weights);
	ADD_ACTION("Interaction probability", "Statistics", s_show_interaction_probability);
	ADD_ACTION("Interaction P value", "Statistics", s_show_interaction_p_value);
	ADD_ACTION("Pathway probability", "Statistics", s_show_pathway_probability);

	ADD_MAIN_MENU("Visualization");

	ADD_MENU("Visualization | Number of interactions", "Number of interactions", "Visualization");
	ADD_ACTION("Heatmap", "Visualization | Number of interactions", s_interaction_number_heatmap);
	ADD_ACTION("Circos Plot", "Visualization | Number of interactions", s_interaction_number_circos);
	ADD_ACTION("Source", "Visualization | Number of interactions", s_interaction_number_source);
	ADD_ACTION("Target", "Visualization | Number of interactions", s_interaction_number_target);

	ADD_MENU("Visualization | Weight of interactions", "Weight of interactions", "Visualization");
	ADD_ACTION("Heatmap", "Visualization | Weight of interactions", s_interaction_weight_heatmap);
	ADD_ACTION("Circos Plot", "Visualization | Weight of interactions", s_interaction_weight_circos);

	ADD_MENU("Visualization | Interaction probability", "Interaction probability", "Visualization");
	ADD_ACTION("Heatmap", "Visualization | Interaction probability", s_interaction_probability_heatmap);
	ADD_ACTION("Circos Plot", "Visualization | Interaction probability", s_interaction_probability_circos);

	ADD_MENU("Visualization | Interaction P Value", "Interaction P Value", "Visualization");
	ADD_ACTION("Heatmap", "Visualization | Interaction P Value", s_interaction_p_value_heatmap);
	ADD_ACTION("Circos Plot", "Visualization | Interaction P Value", s_interaction_p_value_circos);

	ADD_MENU("Visualization | Pathway probability", "Pathway probability", "Visualization");
	ADD_ACTION("Heatmap", "Visualization | Pathway probability", s_pathway_probability_heatmap);
	ADD_ACTION("Circos Plot", "Visualization | Pathway probability", s_pathway_probability_circos);

	ADD_MAIN_MENU("Summary");

	ADD_MENU("Summary | Interaction", "Interaction", "Summary");
	ADD_ACTION("Show", "Summary | Interaction", s_show_interaction_summary);
	ADD_ACTION("Extract", "Summary | Interaction", s_extract_interaction_summary);

	ADD_MENU("Summary | Pathway", "Pathway", "Summary");
	ADD_ACTION("Show", "Summary | Pathway", s_show_pathway_summary);
	ADD_ACTION("Extract", "Summary | Pathway", s_extract_pathway_summary);

	ADD_MAIN_MENU("Export");

	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void CellChatItem::s_interaction_number_source() {

	QStringList valid_sources = this->data()->levels_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose Interaction Source",
		{ "Source" },
		{ soap::InputStyle::ComboBox },
		{ valid_sources }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString level = settings[0];

	int ind = this->data()->levels_.indexOf(level);
	int n_level = valid_sources.size();
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(n_level, n_level);
	data.row(ind) = this->data()->counts_.row(ind).cast<double>();

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data,
		false
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(
		draw_area,
		"Interaction Number of " + level,
		gs
	);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_interaction_number_target() {

	QStringList valid_targets = this->data()->levels_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose Interaction Source",
		{ "Source" },
		{ soap::InputStyle::ComboBox },
		{ valid_targets }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString level = settings[0];

	int ind = this->data()->levels_.indexOf(level);
	int n_level = valid_targets.size();
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(n_level, n_level);
	data.col(ind) = this->data()->counts_.col(ind).cast<double>();

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(
		draw_area,
		"Interaction Number of " + level,
		gs
	);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_interaction_number_circos() {

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->counts_.cast<double>();

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(
		draw_area,
		"Interaction Number",
		gs
	);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_interaction_number_heatmap() {

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->counts_.cast<double>();

	custom_plot::heatmap_plot2(
		draw_area,
		axis_rect,
		10,
		10,
		0,
		this->data()->levels_,
		this->data()->levels_,
		data,
		gs
	);

	custom_plot::set_left_title(axis_rect, "Source", gs, true);
	custom_plot::set_bottom_title(axis_rect, "Target", gs, true);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		data.minCoeff(),
		data.maxCoeff(),
		"Interaction Number",
		gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_show_interaction_numbers() {
	MatrixWindow::show_matrix(
		&this->data()->counts_,
		this->data()->levels_,
		this->data()->levels_,
		"Number of interactions",
		this->signal_emitter_);
};

void CellChatItem::s_interaction_weight_circos() {

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->weights_.cast<double>();

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(
		draw_area,
		"Interaction Weight",
		gs
	);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_interaction_weight_heatmap() {

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->weights_;

	custom_plot::heatmap_plot2(
		draw_area,
		axis_rect,
		10,
		10,
		0,
		this->data()->levels_,
		this->data()->levels_,
		data,
		gs
	);

	custom_plot::set_left_title(axis_rect, "Source", gs, true);
	custom_plot::set_bottom_title(axis_rect, "Target", gs, true);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		data.minCoeff(),
		data.maxCoeff(),
		"Interaction Weight",
		gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_show_interaction_weights() {
	MatrixWindow::show_matrix(
		&this->data()->weights_,
		this->data()->levels_,
		this->data()->levels_,
		"Weights of interactions",
		this->signal_emitter_);
};

void CellChatItem::s_interaction_probability_circos() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->ligand_receptor_index_ }
	);

	if (res.isEmpty())return;

	QString interaction = res[0];

	int index = this->data()->ligand_receptor_index_.indexOf(interaction);

	if (index == -1) {
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->ligand_receptor_probability_[index];

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(
		draw_area,
		"Interaction Probability of " + interaction,
		gs
	);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_interaction_probability_heatmap() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->ligand_receptor_index_ }
	);

	if (res.isEmpty())return;

	QString interaction = res[0];

	int index = this->data()->ligand_receptor_index_.indexOf(interaction);

	if (index == -1) {
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->ligand_receptor_probability_[index];

	custom_plot::heatmap_plot2(
		draw_area,
		axis_rect,
		10,
		10,
		0,
		this->data()->levels_,
		this->data()->levels_,
		data,
		gs
	);

	custom_plot::set_left_title(axis_rect, "Source", gs, true);
	custom_plot::set_bottom_title(axis_rect, "Target", gs, true);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		data.minCoeff(),
		data.maxCoeff(),
		"Probability",
		gs);

	custom_plot::add_title(draw_area, "Interaction Probability of " + interaction, gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_show_interaction_probability() {
	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->ligand_receptor_index_ }
	);

	if (res.isEmpty())return;

	int index = this->data()->ligand_receptor_index_.indexOf(res[0]);

	if (index == -1) {
		return;
	}

	MatrixWindow::show_matrix(
		&this->data()->ligand_receptor_probability_[index],
		this->data()->levels_,
		this->data()->levels_,
		"Probability of " + res[0],
		this->signal_emitter_);
};

void CellChatItem::s_interaction_p_value_circos() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->ligand_receptor_index_ }
	);

	if (res.isEmpty())return;

	QString interaction = res[0];

	int index = this->data()->ligand_receptor_index_.indexOf(interaction);

	if (index == -1) {
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->ligand_receptor_p_value_[index];

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(draw_area, "Interaction P Value of " + interaction, gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_interaction_p_value_heatmap() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->ligand_receptor_index_ }
	);

	if (res.isEmpty())return;

	QString interaction = res[0];

	int index = this->data()->ligand_receptor_index_.indexOf(interaction);

	if (index == -1) {
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->ligand_receptor_p_value_[index];

	custom_plot::heatmap_plot2(
		draw_area,
		axis_rect,
		10,
		10,
		0,
		this->data()->levels_,
		this->data()->levels_,
		data,
		gs
	);

	custom_plot::set_left_title(axis_rect, "Source", gs, true);
	custom_plot::set_bottom_title(axis_rect, "Target", gs, true);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		data.minCoeff(),
		data.maxCoeff(),
		"P Value",
		gs);

	custom_plot::add_title(draw_area, "Interaction P Value of " + interaction, gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_show_interaction_p_value() {
	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->ligand_receptor_index_ }
	);

	if (res.isEmpty())return;

	int index = this->data()->ligand_receptor_index_.indexOf(res[0]);

	if (index == -1) {
		return;
	}

	MatrixWindow::show_matrix(
		&this->data()->ligand_receptor_p_value_[index],
		this->data()->levels_,
		this->data()->levels_,
		"P value of " + res[0],
		this->signal_emitter_);
};

void CellChatItem::s_pathway_probability_circos() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Pathway",
		{ "Pathway" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->pathway_index_ }
	);

	if (res.isEmpty())return;

	QString pathway = res[0];

	if (!this->data()->pathway_probability_.contains(pathway)) {
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->pathway_probability_[pathway];

	auto colors = gs.palette(this->data()->levels_);

	circos_plot(
		draw_area,
		axis_rect,
		this->data()->levels_,
		colors,
		data
	);

	custom_plot::add_round_legend(
		draw_area,
		legend_layout,
		this->data()->levels_,
		colors,
		this->data()->identity_,
		gs);

	custom_plot::add_title(draw_area, "Interaction Probability of " + pathway, gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_pathway_probability_heatmap() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Pathway",
		{ "Pathway" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->pathway_index_ }
	);

	if (res.isEmpty())return;

	QString pathway = res[0];

	if (!this->data()->pathway_probability_.contains(pathway)) {
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	Eigen::MatrixXd data = this->data()->pathway_probability_[pathway];

	custom_plot::heatmap_plot2(
		draw_area,
		axis_rect,
		10,
		10,
		0,
		this->data()->levels_,
		this->data()->levels_,
		data,
		gs
	);

	custom_plot::set_left_title(axis_rect, "Source", gs, true);
	custom_plot::set_bottom_title(axis_rect, "Target", gs, true);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		data.minCoeff(),
		data.maxCoeff(),
		"Probability",
		gs);

	custom_plot::add_title(draw_area, "Interaction Probability of " + pathway, gs);

	this->draw_suite_->update(draw_area);
};

void CellChatItem::s_show_pathway_probability() {

	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Pathway",
		{ "Pathway" },
		{ soap::InputStyle::LineEditWithCompleter },
		{ this->data()->pathway_index_ }
	);

	if (res.isEmpty())return;

	QString pathway = res[0];

	if (!this->data()->pathway_probability_.contains(pathway)) {
		return;
	}

	MatrixWindow::show_matrix(
		&this->data()->pathway_probability_[pathway],
		this->data()->levels_,
		this->data()->levels_,
		"Probability of " + pathway,
		this->signal_emitter_);
};

void CellChatItem::s_extract_interaction_summary() {

	DataFrame* df = new DataFrame(this->data()->interaction_summary_);

	this->signal_emitter_->x_data_create_soon(df, soap::VariableType::DataFrame, "Extracted Interaction Summary");
};

void CellChatItem::s_show_interaction_summary() {

	MatrixWindow::show_matrix(
		&this->data()->interaction_summary_,
		"Interaction Summary",
		this->signal_emitter_
	);
};

void CellChatItem::s_extract_pathway_summary() {

	DataFrame* df = new DataFrame(this->data()->pathway_summary_);

	this->signal_emitter_->x_data_create_soon(df, soap::VariableType::DataFrame, "Extracted Pathway Summary");
};

void CellChatItem::s_show_pathway_summary() {

	MatrixWindow::show_matrix(
		&this->data()->pathway_summary_,
		"Pathway Summary",
		this->signal_emitter_
	);
};
