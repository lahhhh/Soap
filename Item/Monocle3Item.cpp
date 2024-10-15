#include "Monocle3Item.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "CustomPlot.h"
#include "ItemIOWorker.h"
#include "FileWritingWorker.h"

#include "Monocle3ChooseRootDialog.h"

static
Eigen::ArrayXi find_nearest_vertex(const Eigen::MatrixXd& points, const Eigen::MatrixXd& data) {

	int n_point = points.cols(), n_data = data.cols();

	Eigen::ArrayXi nearest_vertex = Eigen::ArrayXi::Zero(n_point);

	for (int i = 0; i < n_point; ++i) {

		double best = (points.col(i).array() - data.col(0).array()).cwiseAbs2().sum();

		int v = 0;

		for (int j = 1; j < n_data; ++j) {

			double local = (points.col(i).array() - data.col(j).array()).cwiseAbs2().sum();

			if (local < best) {
				best = local;
				v = j;
			}

		}

		nearest_vertex[i] = v;

	}

	return nearest_vertex;

}

void Monocle3Item::__s_update_interface() {

	this->setText(2, "Monocle3");
}

void Monocle3Item::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Order Cells", s_order_cells);

	ADD_MAIN_ACTION("Show Pseudo Time", s_pseudo_time_graph);

	ADD_MAIN_ACTION("Pseudo Time Feature Plot", s_pseudo_time_feature_plot);

	ADD_MAIN_ACTION("Feature Plot", s_feature_plot);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void Monocle3Item::pseudo_time_feature_plot_single_cell_rna() {

	auto* single_cell_rna = this->get_root<SingleCellRna>();

	FeatureHandler handler(single_cell_rna);

	auto feature_names = handler.get_feature_names();
	if (feature_names.factor_names.isEmpty()) {
		G_WARN("No Factor for visualization.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Feature Plot Setting",
		{ "Feature", "Factor", "Normalized:yes", "Number of row:1" },
		{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton,
		soap::InputStyle::IntegerLineEdit },
		{ feature_names.numeric_names, feature_names.factor_names }
	);
	if (settings.isEmpty())return;

	QStringList features = multiple_line_edit_with_completer_to_list(settings[0]);

	if (features.isEmpty()) {
		return;
	}

	QString factor_name = settings[1];
	bool normalized = switch_to_bool(settings[2]);
	bool gene_activity = false;

	int nrow = settings[3].toInt();
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

	auto factor_data = handler.get_data({ factor_name });
	if (!factor_data.is_factor()) {
		G_WARN("Illegal Factor.");
		return;
	}

	Eigen::ArrayXd x = this->data()->pseudo_time_;

	for (auto&& d : data) {
		d.slice(this->data()->cell_included_);
	}

	factor_data.slice(this->data()->cell_included_);

	auto draw_area = custom_plot::monocle3_feature_plot(
		data, x, factor_data, "Pseudo Time", nrow, this->draw_suite_->graph_settings_
	);

	this->draw_suite_->update(draw_area);
};

void Monocle3Item::pseudo_time_feature_plot_single_cell_multiome() {

	SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

	FeatureHandler handler(single_cell_multiome);

	auto feature_names = handler.get_feature_names();
	if (feature_names.factor_names.isEmpty()) {
		G_WARN("No Factor for visualization.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Feature Plot Setting",
		{ "Feature", "Factor", "Normalized:yes", "Field", "Number of row:1" },
		{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton,
		soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit },
		{ feature_names.numeric_names, feature_names.factor_names, {"RNA", "ATAC", "Gene Activity"}}
	);
	if (settings.isEmpty())return;

	QStringList features = multiple_line_edit_with_completer_to_list(settings[0]);

	if (features.isEmpty()) {
		return;
	}

	QString factor_name = settings[1];
	bool normalized = switch_to_bool(settings[2]);
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

	auto factor_data = handler.get_data({ factor_name});
	if (!factor_data.is_factor()) {
		G_WARN("Illegal Factor.");
		return;
	}

	Eigen::ArrayXd x = this->data()->pseudo_time_;

	for (auto&& d : data) {
		d.slice(this->data()->cell_included_);
	}

	factor_data.slice(this->data()->cell_included_);

	auto draw_area = custom_plot::monocle3_feature_plot(
		data, x, factor_data, "Pseudo Time", nrow, this->draw_suite_->graph_settings_
	);
	
	this->draw_suite_->update(draw_area);

};

void Monocle3Item::s_pseudo_time_feature_plot() {

	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		this->pseudo_time_feature_plot_single_cell_rna();
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		this->pseudo_time_feature_plot_single_cell_multiome();
	}
	else {
		G_WARN("This is not attached to any object.");
	}
};

void Monocle3Item::s_pseudo_time_graph() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Pseudo Time Graph Settings",
		{ "Show Trajectory:yes", "Trajectory Width:2", "Trajectory Color",
		"Show Nodes:yes", "Nodes Size:5" , "Nodes Color", "Show All Cells:yes" },
		{soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit, soap::InputStyle::ColorChoice,
		soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit, soap::InputStyle::ColorChoice,
		soap::InputStyle::SwitchButton}
	);

	if (settings.isEmpty()) {
		return;
	}

	bool show_trajectory = switch_to_bool(settings[0]);
	int trajectory_width = settings[1].toInt();
	QColor trajectory_color = QColor::fromString(settings[2]);
	bool show_nodes = switch_to_bool(settings[3]);
	int nodes_size = settings[4].toInt();
	QColor nodes_color = QColor::fromString(settings[5]);
	bool show_all_cells = switch_to_bool(settings[6]);

	auto& gs = this->draw_suite_->graph_settings_;

	Eigen::ArrayXd x = this->data()->cell_embedding_.row(0);
	Eigen::ArrayXd y = this->data()->cell_embedding_.row(1);
	Eigen::ArrayXd values = this->data()->pseudo_time_;

	if (values.size() != x.size()) {
		G_WARN("Please Run Order Cells First.");
		return;
	}

	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	custom_plot::set_scatter_plot_axis_style(draw_area, axis_rect,
		this->data()->original_embedding_.data_.colnames_[0],
		this->data()->original_embedding_.data_.colnames_[1], x, y, gs);


	QColor low_color = gs.get_gradient_low_color(Qt::blue);
	QColor middle_color = gs.get_gradient_middle_color(Qt::yellow);
	QColor high_color = gs.get_gradient_high_color(custom_plot::color::darkorchid2);

	custom_plot::patch::scatter_gradient(
		draw_area,
		axis_rect,
		x,
		y,
		values,
		values.minCoeff(),
		values.maxCoeff(),
		low_color,
		middle_color,
		high_color,
		gs.get_scatter_point_size()
	);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		values.minCoeff(),
		values.maxCoeff(),
		"Pseudo Time",
		gs,
		low_color,
		middle_color,
		high_color);

	if (show_trajectory) {

		if (trajectory_width < 1 || trajectory_width > 10) {
			trajectory_width = 2;
			G_LOG("Trajectory Width is reset to 2.");
		}

		igraph_vector_int_t edge_list;
		igraph_vector_int_init(&edge_list, 0);
		igraph_get_edgelist(&this->data()->pr_graph_, &edge_list, 0);

		int n_edge = igraph_vector_int_size(&edge_list) / 2;

		for (int i = 0; i < n_edge; ++i) {
			int v1 = VECTOR(edge_list)[i * 2];
			int v2 = VECTOR(edge_list)[i * 2 + 1];

			custom_plot::patch::line(draw_area, axis_rect,
				this->data()->pr_embedding_(Eigen::all, { v1, v2 }).row(0),
				this->data()->pr_embedding_(Eigen::all, { v1, v2 }).row(1),
				trajectory_color, trajectory_width);
		}
		igraph_vector_int_destroy(&edge_list);
	}

	if (show_nodes) {

		if (nodes_size < 1 || nodes_size > 10) {
			nodes_size = 5;
			G_LOG("Node Size is reset to 5.");
		}

		x = this->data()->pr_embedding_.row(0);
		y = this->data()->pr_embedding_.row(1);

		custom_plot::patch::scatter(
			draw_area,
			axis_rect,
			x,
			y,
			nodes_color,
			nodes_size
		);
	}

	if (show_all_cells) {
		auto rest_index = custom::which(custom::flip(this->data()->cell_included_));
		if (!rest_index.isEmpty()) {
			x = this->data()->original_embedding_.data_.mat_.col(0)(rest_index);
			y = this->data()->original_embedding_.data_.mat_.col(1)(rest_index);
			custom_plot::patch::scatter(
				draw_area,
				axis_rect,
				x,
				y,
				custom_plot::color::gray,
				gs.get_scatter_point_size()
			);

			custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(this->data()->original_embedding_.data_.mat_.col(0),
				this->data()->original_embedding_.data_.mat_.col(1)));
		}
	}

	custom_plot::add_title(draw_area, "Pseudo Time", gs);

	this->draw_suite_->update(draw_area);
};

void Monocle3Item::s_order_cells() {

	G_GETLOCK;

	auto root_nodes = this->get_root_nodes();

	if (root_nodes.size() == 0) {
		G_UNLOCK;
		return;
	}

	this->data()->pseudo_time_ = this->extract_general_graph_ordering(root_nodes);

	G_UNLOCK;
};

Eigen::ArrayXd Monocle3Item::extract_general_graph_ordering(const Eigen::ArrayXi& root_pr_nodes) {

	int n_pr_node = this->data()->pr_embedding_.cols();

	auto closest_vertex = find_nearest_vertex(this->data()->pr_embedding_(Eigen::all, root_pr_nodes), this->data()->cell_embedding_);

	int n_closest_vertex = closest_vertex.size();

	igraph_vector_int_t _closest_vertice;
	igraph_vector_int_init(&_closest_vertice, n_closest_vertex);
	for (int i = 0; i < n_closest_vertex; ++i) {
		VECTOR(_closest_vertice)[i] = closest_vertex[i] + n_pr_node;
	}

	igraph_matrix_t distances;
	igraph_matrix_init(&distances, 0, 0);

	igraph_distances_dijkstra(&this->data()->cell_graph_, &distances, igraph_vss_all(), igraph_vss_vector(&_closest_vertice),
		&this->data()->cell_graph_weights_, IGRAPH_ALL);

	int n_vertice = igraph_matrix_nrow(&distances);
	Eigen::ArrayXd res(n_vertice);
	for (int i = 0; i < n_vertice; ++i) {
		res[i] = MATRIX(distances, i, 0);
	}

	for (int i = 0; i < n_vertice; ++i) {
		for (int j = 1; j < n_closest_vertex; ++j) {
			double dis = MATRIX(distances, i, j);
			if (dis < res[i]) {
				res[i] = dis;
			}
		}
	}

	igraph_vector_int_destroy(&_closest_vertice);
	igraph_matrix_destroy(&distances);

	return res.segment(n_pr_node, n_vertice - n_pr_node);
};

Eigen::ArrayXi Monocle3Item::get_root_nodes() {

	return custom::cast<Eigen::ArrayX>(Monocle3ChooseRootDialog::get_response(this->data()));
};

void Monocle3Item::s_feature_plot() {

	if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		auto* single_cell_multiome = this->get_root<SingleCellMultiome>();

		FeatureHandler handler(single_cell_multiome);

		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Feature Plot Setting",
			{ "Feature", "Normalized:yes", "Scale:no", "Field", "Show Trajectory:yes", "Trajectory Width:2", 
			"Trajectory Color",	"Show Nodes:yes", "Nodes Size:5" , "Nodes Color", "Lighten cell color:yes"},
			{ soap::InputStyle::LineEditWithCompleter, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton,
			soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ColorChoice, soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ColorChoice, soap::InputStyle::SwitchButton },
			{ handler.get_feature_names().numeric_names, {"RNA", "ATAC", "Gene Activity"} }
		);
		if (settings.isEmpty())return;

		QString feature = settings[0];

		bool normalized = switch_to_bool(settings[1]);
		bool scale = switch_to_bool(settings[2]);

		bool gene_activity = settings[3] == "Gene Activity";
		bool show_trajectory = switch_to_bool(settings[4]);
		int trajectory_width = settings[5].toInt();
		QColor trajectory_color = QColor::fromString(settings[6]);
		bool show_nodes = switch_to_bool(settings[7]);
		int nodes_size = settings[8].toInt();
		QColor nodes_color = QColor::fromString(settings[9]);
		bool lighten_color = switch_to_bool(settings[10]);

		auto feature_data = handler.get_data({ feature, normalized, gene_activity });

		if (!feature_data.is_valid()) {
			G_WARN("Invalid Feature");
			return;
		}

		if (feature_data.length() != this->data()->original_embedding_.data_.mat_.rows()) {
			G_WARN("Unmatched Length!");
			return;
		}

		if (lighten_color) {
			feature_data.message.insert("Lighten Color");
		}

		auto [draw_area, axis_rect, _] = custom_plot::feature_plot(feature_data, &this->data()->original_embedding_, scale, this->draw_suite_->graph_settings_);

		if (show_trajectory) {

			if (trajectory_width < 1 || trajectory_width > 10) {
				trajectory_width = 3;
				G_LOG("Trajectory Width is reset to 3.");
			}

			igraph_vector_int_t edge_list;
			igraph_vector_int_init(&edge_list, 0);
			igraph_get_edgelist(&this->data()->pr_graph_, &edge_list, 0);

			int n_edge = igraph_vector_int_size(&edge_list) / 2;

			for (int i = 0; i < n_edge; ++i) {
				int v1 = VECTOR(edge_list)[i * 2];
				int v2 = VECTOR(edge_list)[i * 2 + 1];

				custom_plot::patch::line(draw_area, axis_rect,
					this->data()->pr_embedding_(Eigen::all, { v1, v2 }).row(0),
					this->data()->pr_embedding_(Eigen::all, { v1, v2 }).row(1),
					trajectory_color, trajectory_width);
			}
			igraph_vector_int_destroy(&edge_list);
		}

		if (show_nodes) {

			if (nodes_size < 1 || nodes_size > 10) {
				nodes_size = 3;
				G_LOG("Node Size is reset to 3.");
			}

			Eigen::ArrayXd x = this->data()->pr_embedding_.row(0);
			Eigen::ArrayXd y = this->data()->pr_embedding_.row(1);

			custom_plot::patch::scatter(
				draw_area,
				axis_rect,
				x,
				y,
				nodes_color,
				nodes_size
			);
		}

		this->draw_suite_->update(draw_area);
	}
	else if (this->stem_from(soap::VariableType::SingleCellRna)) {
		auto* single_cell_rna = this->get_root<SingleCellRna>();

		FeatureHandler handler(single_cell_rna);

		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Feature Plot Setting",
			{ "Feature", "Normalized:yes", "Scale:no", "Show Trajectory:yes", "Trajectory Width:2",
			"Trajectory Color",	"Show Nodes:yes", "Nodes Size:5" , "Nodes Color", "Lighten cell color:yes" },
			{ soap::InputStyle::LineEditWithCompleter, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton,
			soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ColorChoice, soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ColorChoice, soap::InputStyle::SwitchButton },
			{ handler.get_feature_names().numeric_names }
		);
		if (settings.isEmpty())return;

		QString feature = settings[0];

		bool normalized = switch_to_bool(settings[1]);
		bool scale = switch_to_bool(settings[2]);

		bool show_trajectory = switch_to_bool(settings[3]);
		int trajectory_width = settings[4].toInt();
		QColor trajectory_color = QColor::fromString(settings[5]);
		bool show_nodes = switch_to_bool(settings[6]);
		int nodes_size = settings[7].toInt();
		QColor nodes_color = QColor::fromString(settings[8]);
		bool lighten_color = switch_to_bool(settings[9]);

		auto feature_data = handler.get_data({ feature, normalized, false });

		if (!feature_data.is_valid()) {
			G_WARN("Invalid Feature");
			return;
		}

		if (feature_data.length() != this->data()->original_embedding_.data_.mat_.rows()) {
			G_WARN("Unmatched Length!");
			return;
		}

		if (lighten_color) {
			feature_data.message.insert("Lighten Color");
		}

		auto [draw_area, axis_rect, _] = custom_plot::feature_plot(feature_data, &this->data()->original_embedding_, scale, this->draw_suite_->graph_settings_);

		if (show_trajectory) {

			if (trajectory_width < 1 || trajectory_width > 10) {
				trajectory_width = 3;
				G_LOG("Trajectory Width is reset to 3.");
			}

			igraph_vector_int_t edge_list;
			igraph_vector_int_init(&edge_list, 0);
			igraph_get_edgelist(&this->data()->pr_graph_, &edge_list, 0);

			int n_edge = igraph_vector_int_size(&edge_list) / 2;

			for (int i = 0; i < n_edge; ++i) {
				int v1 = VECTOR(edge_list)[i * 2];
				int v2 = VECTOR(edge_list)[i * 2 + 1];

				custom_plot::patch::line(draw_area, axis_rect,
					this->data()->pr_embedding_(Eigen::all, { v1, v2 }).row(0),
					this->data()->pr_embedding_(Eigen::all, { v1, v2 }).row(1),
					trajectory_color, trajectory_width);
			}
			igraph_vector_int_destroy(&edge_list);
		}

		if (show_nodes) {

			if (nodes_size < 1 || nodes_size > 10) {
				nodes_size = 3;
				G_LOG("Node Size is reset to 3.");
			}

			Eigen::ArrayXd x = this->data()->pr_embedding_.row(0);
			Eigen::ArrayXd y = this->data()->pr_embedding_.row(1);

			custom_plot::patch::scatter(
				draw_area,
				axis_rect,
				x,
				y,
				nodes_color,
				nodes_size
			);
		}

		this->draw_suite_->update(draw_area);
	}
};