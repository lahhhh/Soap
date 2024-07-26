#include "EmbeddingViewWindow.h"
#include "CommonDialog.h"
#include "PlotWindow.h"
#include "CellMarkerDatabase.h"
#include "CustomPlot.h"
#include "MatrixWindow.h"
#include "GraphSettingDialog.h"

#include "SoapGUI.h"

#include "EnrichmentUtility.h"

inline bool is_in_triangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3) {
	if ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) < 0) {
		return is_in_triangle(x, y, x1, y1, x3, y3, x2, y2);
	}
	if ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) > 0
		&& (x3 - x2) * (y - y2) - (y3 - y2) * (x - x2) > 0
		&& (x1 - x3) * (y - y3) - (y1 - y3) * (x - x3) > 0)
	{
		return true;
	}
	return false;
}

Eigen::ArrayX<bool> EmbeddingViewWindow::get_selected() {

	Eigen::ArrayXd first_dimension = this->embedding_->data_.mat_.col(0);
	Eigen::ArrayXd second_dimension = this->embedding_->data_.mat_.col(1);
	int size = this->point_x_.size() - 2;
	int n_point = first_dimension.size();
	Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(n_point, false);

	for (int i = 0; i < n_point; ++i) {

		double x = first_dimension[i], y = second_dimension[i];

		for (int j = 0; j < size; ++j) {

			double x1 = this->point_x_[0], x2 = this->point_x_[j + 1], x3 = this->point_x_[j + 2];
			double y1 = this->point_y_[0], y2 = this->point_y_[j + 1], y3 = this->point_y_[j + 2];

			if (is_in_triangle(x, y, x1, y1, x2, y2, x3, y3)) {
				filter[i] = true;
				break;
			}
		}
	}

	return filter;
};

EmbeddingViewWindow::EmbeddingViewWindow(
	SingleCellRna* single_cell_rna,
	Embedding* embedding,
	SignalEmitter* signal_emitter
) :
	QMainWindow(signal_emitter->widget_),
	handler_(single_cell_rna),
	signal_emitter_(signal_emitter),
	embedding_(embedding)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &EmbeddingViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &EmbeddingViewWindow::s_check_data);

	this->set_layout();

	this->preload();

	this->set_property();
}

void EmbeddingViewWindow::view(
	SingleCellRna* single_cell_rna,
	Embedding* embedding,
	SignalEmitter* signal_emitter
) {
	EmbeddingViewWindow* window = new EmbeddingViewWindow(single_cell_rna, embedding, signal_emitter);
};

EmbeddingViewWindow::EmbeddingViewWindow(
	SingleCellAtac* single_cell_atac,
	Embedding* embedding,
	SignalEmitter* signal_emitter
) :
	QMainWindow(signal_emitter->widget_),
	handler_(single_cell_atac),
	signal_emitter_(signal_emitter),
	embedding_(embedding)
{
	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &EmbeddingViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &EmbeddingViewWindow::s_check_data);

	this->set_layout();

	this->preload();

	this->set_property();
}

void EmbeddingViewWindow::view(
	SingleCellAtac* single_cell_atac,
	Embedding* embedding,
	SignalEmitter* signal_emitter
) {
	EmbeddingViewWindow* window = new EmbeddingViewWindow(single_cell_atac, embedding, signal_emitter);
};

EmbeddingViewWindow::EmbeddingViewWindow(
	SingleCellMultiome* single_cell_multiome,
	Embedding* embedding,
	SignalEmitter* signal_emitter
) :
	QMainWindow(signal_emitter->widget_),
	handler_(single_cell_multiome),
	signal_emitter_(signal_emitter),
	embedding_(embedding)
{

	connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &EmbeddingViewWindow::s_check_data);
	connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &EmbeddingViewWindow::s_check_data);

	this->set_layout();

	this->preload();

	this->set_property();
}

void EmbeddingViewWindow::view(
	SingleCellMultiome* single_cell_multiome,
	Embedding* embedding,
	SignalEmitter* signal_emitter
) {
	EmbeddingViewWindow* window = new EmbeddingViewWindow(single_cell_multiome, embedding, signal_emitter);
};

void EmbeddingViewWindow::set_layout() {

	this->main_interface_ = new QWidget(this);
	this->setCentralWidget(this->main_interface_);

	this->main_layout_ = new QHBoxLayout;
	this->main_interface_->setLayout(this->main_layout_);

	this->set_left_layout();

	this->set_right_layout();

	this->final_connect();

	this->draw_suite_ = new PlotsSuite();

	connect(this->graph_setting_switch_, &Switch::toggled, this->draw_suite_, &PlotsSuite::s_setting_activate);
	connect(this->draw_suite_, &PlotsSuite::x_plot_prepared, this, &EmbeddingViewWindow::s_new_plot);

	this->draw_suite_->prepare();
};

void EmbeddingViewWindow::set_property() {

	this->setAttribute(Qt::WA_DeleteOnClose);

	G_SET_ICON;

	this->resize(1200, 800);

	this->setWindowTitle("Embedding View");

	this->show();
};

void EmbeddingViewWindow::s_check_data(void* data, soap::VariableType type, void* item) {

	if (data == this->handler_.data_) {
		this->close();
	}

};

void EmbeddingViewWindow::final_connect() {

	connect(this->cell_marker_tree_widget_, &CellMarkerTreeWidget::x_gene_to_active, this->active_tree_widget_, &CellMarkerTreeWidget::s_receive_marker);
	connect(this->cell_marker_tree_widget_, &CellMarkerTreeWidget::x_gene_to_alternative, this->alternative_tree_widget_, &CellMarkerTreeWidget::s_receive_marker);
	connect(this->cell_marker_tree_widget_, &CellMarkerTreeWidget::x_type_to_active, this->active_tree_widget_, &CellMarkerTreeWidget::s_receive_type);
	connect(this->cell_marker_tree_widget_, &CellMarkerTreeWidget::x_type_to_alternative, this->alternative_tree_widget_, &CellMarkerTreeWidget::s_receive_type);

	connect(this->alternative_tree_widget_, &CellMarkerTreeWidget::x_gene_to_active, this->active_tree_widget_, &CellMarkerTreeWidget::s_receive_marker);
	connect(this->alternative_tree_widget_, &CellMarkerTreeWidget::x_type_to_active, this->active_tree_widget_, &CellMarkerTreeWidget::s_receive_type);

	connect(this->active_tree_widget_, &CellMarkerTreeWidget::x_gene_to_alternative, this->alternative_tree_widget_, &CellMarkerTreeWidget::s_receive_marker);
	connect(this->active_tree_widget_, &CellMarkerTreeWidget::x_type_to_alternative, this->alternative_tree_widget_, &CellMarkerTreeWidget::s_receive_type);

	connect(this, &EmbeddingViewWindow::x_gene_to_active, this->active_tree_widget_, &CellMarkerTreeWidget::s_receive_marker);
	connect(this, &EmbeddingViewWindow::x_gene_to_alternative, this->alternative_tree_widget_, &CellMarkerTreeWidget::s_receive_marker);
	connect(this, &EmbeddingViewWindow::x_type_to_active, this->active_tree_widget_, &CellMarkerTreeWidget::s_receive_type);
	connect(this, &EmbeddingViewWindow::x_type_to_alternative, this->alternative_tree_widget_, &CellMarkerTreeWidget::s_receive_type);
};

EmbeddingViewWindow::~EmbeddingViewWindow() {

	delete this->draw_suite_;
}
void EmbeddingViewWindow::s_update_cell_marker() {

	QString db = this->cell_marker_database_box_->currentText();

	QString tissue = this->cell_marker_tissue_box_->currentText();

	auto&& types = this->cell_marker_database_[db][tissue];

	this->cell_marker_tree_widget_->set_types(types);
};

void EmbeddingViewWindow::show_type(const QString& type_name, const QStringList& features) {

	QString legend_title = "Expression";

	bool normalize = this->normalize_switch_->value_;
	bool use_gene_activity{ false };
	if (this->gene_activity_switch_ != nullptr) {
		use_gene_activity = this->gene_activity_switch_->value_;
	}

	if (use_gene_activity) {
		legend_title = "Activity";
	}

	Eigen::ArrayXd data;

	int n_not_found{ 0 }, n_found{ 0 };
	for (auto&& feature : features) {
		auto feature_data = this->handler_.get_data({ feature, normalize, use_gene_activity });

		if (!feature_data.is_continuous()) {
			++n_not_found;
			continue;
		}

		if (data.size() == 0) {
			data = _Cs cast<Eigen::ArrayX>(feature_data.get_continuous());
		}
		else {
			data += _Cs cast<Eigen::ArrayX>(feature_data.get_continuous());
		}

		++n_found;
	}
	if (n_found == 0) {
		G_WARN("Found no feature in " + type_name);
		return;
	}
	else {
		G_LOG("Found " + QString::number(n_found) + " features in " + type_name + ", " + QString::number(n_not_found) + " not found.");
	}

	data /= n_found;

	QUERY_DATA d;
	d.name = type_name;
	d.type = QUERY_DATA::DataType::numeric;
	d.dd = _Cs cast<QVector>(data);
	d.info["Legend Title"] = legend_title;

	auto [draw_area, axis_rect, _] = _Cp feature_plot(d, this->embedding_, false, this->draw_suite_->graph_settings_);

	this->plot_to_axis_rect_[draw_area] = axis_rect;
	this->draw_suite_->update(draw_area);
};

void EmbeddingViewWindow::s_show_type(const QString& type_name, const QStringList& markers) {

	this->show_type(type_name, markers);
};

void EmbeddingViewWindow::s_show_gene(const QString& gene_name) {

	this->show_gene(gene_name);
};

void EmbeddingViewWindow::preload() {

	soap::Species species = this->handler_.get_species();

	this->cell_marker_database_ = CellMarkerDatabase::get_database(this->handler_.get_species());
	this->cell_marker_database_box_->addItems(this->cell_marker_database_.keys());
	this->cell_marker_database_box_->adjustSize();
	this->cell_marker_database_box_->setFixedHeight(30);

	QString db = this->cell_marker_database_box_->currentText();
	this->cell_marker_tissue_box_->addItems(this->cell_marker_database_[db].keys());
	this->cell_marker_tissue_box_->adjustSize();
	this->cell_marker_tissue_box_->setFixedHeight(30);

	QString tissue = this->cell_marker_tissue_box_->currentText();

	auto&& types = this->cell_marker_database_[db][tissue];

	this->cell_marker_tree_widget_->set_types(types);

	this->pathway_content_ = get_pathway_information(this->handler_.get_species());

	QStringList pathways;
	for (const auto& database : this->pathway_content_) {
		pathways << database.keys();
	}

	this->pathway_line_edit_->setCompleter(new QCompleter(pathways, this));

	this->database_box_->addItems(this->pathway_content_.keys());
	this->database_box_->adjustSize();
	this->database_box_->setFixedHeight(30);
};

void EmbeddingViewWindow::s_view_pathway() {

	QString pathway_name = this->pathway_line_edit_->text();
	QStringList gene_list;
	bool found = false;

	for (const auto& database : this->pathway_content_) {
		if (database.contains(pathway_name)) {
			gene_list = database[pathway_name];
			found = true;
			break;
		}
	}

	if (!found) {
		G_LOG(pathway_name + " is not found in database!");
		return;
	}

	this->show_type(pathway_name, gene_list);
};

void EmbeddingViewWindow::s_explore_pathway() {

	QString database_name = this->database_box_->currentText();

	auto tmp = this->pathway_content_[database_name].keys();

	MatrixWindow::show_matrix(&tmp, "Pathway of " + database_name, this->signal_emitter_);
};

void EmbeddingViewWindow::set_left_layout() {

	this->left_layout_ = new QGridLayout();

	set_left_top_left_layout();
	set_left_top_right_layout();
	set_left_bottom_right_layout();

	this->information_area_ = new InformationTextBrowser(this);
	this->left_layout_->addWidget(this->information_area_, 1, 0);

	this->main_layout_->addLayout(this->left_layout_);
};

void EmbeddingViewWindow::s_clear_active_items() {
	this->active_tree_widget_->safe_clear();
};

void EmbeddingViewWindow::s_select_cell() {
	this->plot_interacting_ = !this->plot_interacting_;

	if (this->plot_interacting_) {
		this->connections_ << connect(this->draw_area_, &QCustomPlot::mousePress, this, &EmbeddingViewWindow::s_plot_mouse_press);
		this->connections_ << connect(this->draw_area_, &QCustomPlot::mouseMove, this, &EmbeddingViewWindow::s_plot_mouse_move);
		this->connections_ << connect(this->draw_area_, &QCustomPlot::mouseDoubleClick, this, &EmbeddingViewWindow::s_plot_mouse_double_click);
	}
	else {
		for (const auto& connection : this->connections_) {
			this->disconnect(connection);
		}
		this->connections_.clear();

		for (auto p : this->plot_active_shape_.keys()) {
			p->removePlottable(this->plot_active_shape_[p]);
		}
		this->plot_active_shape_.clear();

		if (this->point_x_.size() < 2) {
			this->clear_interactive();
		}
		else {
			if (this->handler_.type_ == FeatureHandler::DataType::SingleCellRna) {
				QStringList settings = CommonDialog::get_response(
					this->signal_emitter_,
					"Metadata Assign",
					{ "Metadata", "Value" },
					{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
					{ this->handler_.get_metadata_names() }
				);

				if (!settings.isEmpty()) {
					auto filter = get_selected();
					auto* metadata = static_cast<SingleCellRna*>(this->handler_.data_)->metadata();
					this->signal_emitter_->x_data_edit_soon(metadata);
					metadata->mat_.edit(settings[0], filter, settings[1]);
				}
			}
			else if (this->handler_.type_ == FeatureHandler::DataType::SingleCellAtac) {
				QStringList settings = CommonDialog::get_response(
					this->signal_emitter_,
					"Metadata Assign",
					{ "Metadata", "Value" },
					{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
					{ this->handler_.get_metadata_names() }
				);

				if (!settings.isEmpty()) {
					auto filter = get_selected();
					auto* metadata = static_cast<SingleCellAtac*>(this->handler_.data_)->metadata();
					this->signal_emitter_->x_data_edit_soon(metadata);
					metadata->mat_.edit(settings[0], filter, settings[1]);
				}
			}
			else if (this->handler_.type_ == FeatureHandler::DataType::SingleCellMultiome) {
				QStringList settings = CommonDialog::get_response(
					this->signal_emitter_,
					"Metadata Assign",
					{ "Metadata", "Value" },
					{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
					{ this->handler_.get_metadata_names() }
				);

				if (!settings.isEmpty()) {
					auto filter = get_selected();
					auto* metadata = static_cast<SingleCellMultiome*>(this->handler_.data_)->metadata();
					this->signal_emitter_->x_data_edit_soon(metadata);
					metadata->mat_.edit(settings[0], filter, settings[1]);
				}
			}

			this->point_x_.clear();
			this->point_y_.clear();

			for (auto p : this->plot_to_graph_.keys()) {
				p->removeGraph(this->plot_to_graph_[p]);
			}
			this->plot_to_graph_.clear();

			for (auto p : this->plot_shapes_.keys()) {
				for (auto s : this->plot_shapes_[p]) {
					p->removePlottable(s);
				}
			}
			this->plot_shapes_.clear();
		}
		this->draw_area_->replot();
	}
};

void EmbeddingViewWindow::s_plot_mouse_double_click(QMouseEvent*) {

	this->plot_interacting_ = false;

	for (const auto& connection : this->connections_) {
		this->disconnect(connection);
	}

	this->connections_.clear();

	if (this->point_x_.size() < 2) {
		clear_interactive();
	}
	else {
		if (this->handler_.type_ == FeatureHandler::DataType::SingleCellRna) {
			QStringList settings = CommonDialog::get_response(
				this->signal_emitter_,
				"Metadata Assign",
				{ "Metadata", "Value" },
				{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
				{ this->handler_.get_metadata_names() }
			);

			if (!settings.isEmpty()) {
				auto filter = get_selected();
				auto* metadata = static_cast<SingleCellRna*>(this->handler_.data_)->metadata();
				this->signal_emitter_->x_data_edit_soon(metadata);
				metadata->mat_.edit(settings[0], filter, settings[1]);
			}
		}
		else if (this->handler_.type_ == FeatureHandler::DataType::SingleCellAtac) {
			QStringList settings = CommonDialog::get_response(
				this->signal_emitter_,
				"Metadata Assign",
				{ "Metadata", "Value" },
				{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
				{ this->handler_.get_metadata_names() }
			);

			if (!settings.isEmpty()) {
				auto filter = get_selected();
				auto* metadata = static_cast<SingleCellAtac*>(this->handler_.data_)->metadata();
				this->signal_emitter_->x_data_edit_soon(metadata);
				metadata->mat_.edit(settings[0], filter, settings[1]);
			}
		}
		else if (this->handler_.type_ == FeatureHandler::DataType::SingleCellMultiome) {
			QStringList settings = CommonDialog::get_response(
				this->signal_emitter_,
				"Metadata Assign",
				{ "Metadata", "Value" },
				{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
				{ this->handler_.get_metadata_names() }
			);

			if (!settings.isEmpty()) {
				auto filter = get_selected();
				auto* metadata = static_cast<SingleCellMultiome*>(this->handler_.data_)->metadata();
				this->signal_emitter_->x_data_edit_soon(metadata);
				metadata->mat_.edit(settings[0], filter, settings[1]);
			}
		}

		this->point_x_.clear();
		this->point_y_.clear();

		for (auto p : this->plot_active_shape_.keys()) {
			p->removePlottable(this->plot_active_shape_[p]);
		}
		this->plot_active_shape_.clear();

		for (auto p : this->plot_to_graph_.keys()) {
			p->removeGraph(this->plot_to_graph_[p]);
		}
		this->plot_to_graph_.clear();

		for (auto p : this->plot_shapes_.keys()) {
			for (auto s : this->plot_shapes_[p]) {
				p->removePlottable(s);
			}
		}
		this->plot_shapes_.clear();
	}
	this->draw_area_->replot();
};

void EmbeddingViewWindow::s_plot_mouse_move(QMouseEvent* event) {

	if (!this->plot_to_axis_rect_.contains(this->draw_area_) || this->point_x_.isEmpty())
		return;

	int x_pos = event->pos().x(), y_pos = event->pos().y();

	QCPAxisRect* axis_rect = this->plot_to_axis_rect_[this->draw_area_];

	double x = axis_rect->axis(QCPAxis::atBottom)->pixelToCoord(x_pos), y = axis_rect->axis(QCPAxis::atLeft)->pixelToCoord(y_pos);

	if (this->point_x_.size() == 1) {
		QVector<double> X(this->point_x_), Y(this->point_y_);
		X << x;
		Y << y;

		QCPGraph* graph = this->plot_to_graph_[this->draw_area_];

		QPen pen;
		pen.setColor(Qt::black);
		pen.setWidth(2);
		pen.setStyle(Qt::DashLine);

		graph->setPen(pen);
		graph->setLineStyle(QCPGraph::lsLine);
		graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));
		graph->setData(X, Y);

		this->draw_area_->replot();
	}
	else {
		// remove last moving shape
		if (this->plot_to_graph_.contains(this->draw_area_)) {
			this->draw_area_->removeGraph(this->plot_to_graph_[this->draw_area_]);
			this->plot_to_graph_.remove(this->draw_area_);
		}
		if (this->plot_active_shape_.contains(this->draw_area_)) {
			this->draw_area_->removePlottable(this->plot_active_shape_[this->draw_area_]);
			this->plot_active_shape_.remove(this->draw_area_);
		}

		QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		this->plot_active_shape_[this->draw_area_] = shape;

		int length = this->point_x_.size();

		QVector<QCPCurveData> data(3);
		data[0] = QCPCurveData(0, this->point_x_[0], this->point_y_[0]);
		data[1] = QCPCurveData(1, this->point_x_[length - 1], this->point_y_[length - 1]);
		data[2] = QCPCurveData(2, x, y);
		shape->setPen(Qt::NoPen);
		shape->setBrush(QBrush(QColor(221, 221, 221, 128)));
		shape->data()->set(data, true);

		this->draw_area_->replot();
	}
};

void EmbeddingViewWindow::s_plot_mouse_press(QMouseEvent* event) {
	if (!this->plot_to_axis_rect_.contains(this->draw_area_))
		return;

	int x_pos = event->pos().x(), y_pos = event->pos().y();

	QCPAxisRect* axis_rect = this->plot_to_axis_rect_[this->draw_area_];

	double x = axis_rect->axis(QCPAxis::atBottom)->pixelToCoord(x_pos), y = axis_rect->axis(QCPAxis::atLeft)->pixelToCoord(y_pos);
	this->point_x_ << x;
	this->point_y_ << y;

	if (this->point_x_.size() == 1) {
		// first dot
		QCPGraph* graph = this->draw_area_->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		this->plot_to_graph_[this->draw_area_] = graph;

		QPen pen;
		pen.setColor(Qt::black);
		pen.setWidth(2);

		graph->setPen(pen);
		graph->setLineStyle(QCPGraph::lsNone);
		graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));
		graph->setData(this->point_x_, this->point_y_);

		this->draw_area_->replot();
	}
	else if (this->point_x_.size() == 2) {
		// first border
		this->draw_area_->removeGraph(this->plot_to_graph_[this->draw_area_]);

		QCPGraph* graph = this->draw_area_->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		this->plot_to_graph_[this->draw_area_] = graph;

		QPen pen;
		pen.setColor(Qt::black);
		pen.setWidth(2);
		pen.setStyle(Qt::DashLine);

		graph->setPen(pen);
		graph->setLineStyle(QCPGraph::lsLine);
		graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 3));
		graph->setData(this->point_x_, this->point_y_);
		this->draw_area_->replot();
	}
	else {
		// fix triangle
		if (this->plot_active_shape_.contains(this->draw_area_)) {
			this->draw_area_->removePlottable(this->plot_active_shape_[this->draw_area_]);
			this->plot_active_shape_.remove(this->draw_area_);
		}

		QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		this->plot_shapes_[this->draw_area_] << shape;
		int length = this->point_x_.size();

		QVector<QCPCurveData> data(3);
		data[0] = QCPCurveData(0, this->point_x_[0], this->point_y_[0]);

		for (int j = 1; j < 3; ++j) {
			data[j] = QCPCurveData(j, this->point_x_[length - 3 + j], this->point_y_[length - 3 + j]);
		}

		shape->setPen(Qt::NoPen);
		shape->setBrush(QBrush(QColor(221, 221, 221, 128)));
		shape->data()->set(data, true);

		this->draw_area_->replot();
	}
};

void EmbeddingViewWindow::set_right_layout() {
	this->right_layout_ = new QVBoxLayout;

	G_SET_BUTTON(this->select_cell_button_, "Select Cell", soap::MiddleSize);
	connect(this->select_cell_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_select_cell);

	G_SET_BUTTON(this->graph_setting_button_, "Graph Settings", soap::MiddleSize);
	connect(this->graph_setting_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_set_graph_settings);

	G_SET_SWITCH(this->graph_setting_switch_, false, this->graph_setting_button_, soap::MiddleSize);

	G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->right_layout_, row_layout, this->select_cell_button_, this->graph_setting_button_, this->graph_setting_switch_);

	G_SET_PLOTSUITE_BUTTON;

	row_layout = new QHBoxLayout;

	row_layout->addStretch();
	row_layout->addWidget(this->previous_picture_button_);
	row_layout->addWidget(this->next_picture_button_);
	row_layout->addWidget(this->pop_picture_button_);
	row_layout->addWidget(this->clear_picture_button_);
	row_layout->addWidget(this->save_picture_button_);

	this->right_layout_->addLayout(row_layout);

	this->save_picture_menu_ = new QMenu();

	this->save_png_action_ = new QAction("PNG", this);
	this->save_jpg_action_ = new QAction("JPG", this);
	this->save_bmp_action_ = new QAction("BMP", this);
	this->save_pdf_action_ = new QAction("PDF", this);
	this->save_pdf_and_png_action_ = new QAction("PDF and PNG", this);

	this->save_picture_menu_->addAction(this->save_png_action_);
	this->save_picture_menu_->addAction(this->save_jpg_action_);
	this->save_picture_menu_->addAction(this->save_bmp_action_);
	this->save_picture_menu_->addAction(this->save_pdf_action_);
	this->save_picture_menu_->addAction(this->save_pdf_and_png_action_);

	connect(this->save_jpg_action_, &QAction::triggered, this, &EmbeddingViewWindow::s_save_jpg);
	connect(this->save_png_action_, &QAction::triggered, this, &EmbeddingViewWindow::s_save_png);
	connect(this->save_bmp_action_, &QAction::triggered, this, &EmbeddingViewWindow::s_save_bmp);
	connect(this->save_pdf_action_, &QAction::triggered, this, &EmbeddingViewWindow::s_save_pdf);
	connect(this->save_pdf_and_png_action_, &QAction::triggered, this, &EmbeddingViewWindow::s_save_pdf_and_png);

	this->save_picture_button_->setMenu(this->save_picture_menu_);

	this->main_layout_->addLayout(this->right_layout_);

	connect(this->previous_picture_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_previous_plot);
	connect(this->next_picture_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_next_plot);
	connect(this->pop_picture_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_pop_plot);
	connect(this->clear_picture_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_clear_plot);
};

void EmbeddingViewWindow::set_left_top_left_layout() {
	this->left_top_left_layout_ = new QVBoxLayout;

	G_SET_LABEL_PRECISE(this->cell_marker_area_label_, "Cell Markers", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));
	this->left_top_left_layout_->addWidget(this->cell_marker_area_label_);

	G_SET_LABEL(this->cell_marker_database_label_, "Database", soap::MiddleSize);
	this->cell_marker_database_box_ = new QComboBox(this);

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_left_layout_, row_layout, this->cell_marker_database_label_, this->cell_marker_database_box_);

	G_SET_LABEL(this->cell_marker_tissue_label_, "Tissue", soap::MiddleSize);
	this->cell_marker_tissue_box_ = new QComboBox(this);

	G_ADD_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_left_layout_, row_layout, this->cell_marker_tissue_label_, this->cell_marker_tissue_box_);

	this->cell_marker_tree_widget_ = new CellMarkerTreeWidget(CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget, this);
	this->left_top_left_layout_->addWidget(this->cell_marker_tree_widget_);

	this->left_layout_->addLayout(this->left_top_left_layout_, 0, 0);

	connect(this->cell_marker_database_box_, &QComboBox::currentIndexChanged, this, &EmbeddingViewWindow::s_update_tissue);
	connect(this->cell_marker_tissue_box_, &QComboBox::currentIndexChanged, this, &EmbeddingViewWindow::s_update_cell_marker);

	connect(this->cell_marker_tree_widget_, &CellMarkerTreeWidget::x_show_gene, this, &EmbeddingViewWindow::s_show_gene);
	connect(this->cell_marker_tree_widget_, &CellMarkerTreeWidget::x_show_type, this, &EmbeddingViewWindow::s_show_type);
};

void EmbeddingViewWindow::s_update_tissue() {
	this->cell_marker_tissue_box_->blockSignals(true);

	this->cell_marker_tissue_box_->clear();
	this->cell_marker_tissue_box_->addItems(this->cell_marker_database_[this->cell_marker_database_box_->currentText()].keys());
	this->cell_marker_tissue_box_->adjustSize();

	QString db = this->cell_marker_database_box_->currentText();

	QString tissue = this->cell_marker_tissue_box_->currentText();

	auto&& types = this->cell_marker_database_[db][tissue];

	this->cell_marker_tree_widget_->set_types(types);

	this->cell_marker_tissue_box_->blockSignals(false);
};

void EmbeddingViewWindow::show_gene(const QString& gene_name) {

	bool normalize = this->normalize_switch_->value_;
	bool use_gene_activity{ false };
	if (this->gene_activity_switch_ != nullptr) {
		use_gene_activity = this->gene_activity_switch_->value_;
	}

	auto feature_data = this->handler_.get_data({ gene_name, normalize, use_gene_activity });

	if (!feature_data.is_continuous()) {
		G_WARN("Feature not found.");
		return;
	}

	auto [draw_area, axis_rect, _] = _Cp feature_plot(feature_data, this->embedding_, false, this->draw_suite_->graph_settings_);

	this->plot_to_axis_rect_[draw_area] = axis_rect;
	this->draw_suite_->update(draw_area);
};

void EmbeddingViewWindow::s_view_gene() {

	QString gene_name = this->gene_line_edit_->text();

	this->show_gene(gene_name);
};

void EmbeddingViewWindow::s_view_metadata() {

	QString feature = this->metadata_box_->currentText();

	auto feature_data = this->handler_.get_data({ feature });

	if (!feature_data.is_valid()) {
		G_WARN("Illegal data.");
		return;
	}

	auto [draw_area, axis_rect, _] = _Cp feature_plot(feature_data, this->embedding_, false, this->draw_suite_->graph_settings_);

	this->plot_to_axis_rect_[draw_area] = axis_rect;
	this->draw_suite_->update(draw_area);
};

void EmbeddingViewWindow::s_alternative_gene() {

	QString gene_name = this->gene_line_edit_->text();

	emit x_gene_to_alternative(gene_name);
};

void EmbeddingViewWindow::s_active_gene() {

	QString gene_name = this->gene_line_edit_->text();

	emit x_gene_to_active(gene_name);
};

void EmbeddingViewWindow::set_left_top_right_layout() {

	this->left_top_right_layout_ = new QVBoxLayout;

	G_SET_LABEL_PRECISE(this->gene_label_, "Gene", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));
	this->left_top_right_layout_->addWidget(this->gene_label_);


	G_SET_LABEL(this->gene_name_label_, "Gene Name", soap::MiddleSize);
	G_SET_LINEEDIT_WITH_COMPLETER(this->gene_line_edit_, "", this->handler_.get_feature_names().numeric_names, soap::MiddleSize);
	G_SET_BUTTON_ICON(this->gene_view_button_, FILE_EYE_ICON_PNG, QSize(30, 30));

	G_SET_BUTTON(this->gene_alternative_button_, "►", QSize(30, 30));
	G_SET_BUTTON(this->gene_active_button_, "▼", QSize(30, 30));

	G_ADD_NEW_FIVE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_right_layout_, row_layout, this->gene_name_label_,
		this->gene_line_edit_, this->gene_view_button_, this->gene_alternative_button_, this->gene_active_button_);

	connect(this->gene_view_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_view_gene);
	connect(this->gene_alternative_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_alternative_gene);
	connect(this->gene_active_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_active_gene);

	G_SET_LABEL_PRECISE(this->metadata_label_, "Metadata", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));
	this->left_top_right_layout_->addWidget(this->metadata_label_);

	G_SET_LABEL(this->metadata_name_label_, "Metadata Name", soap::MiddleSize);
	this->metadata_box_ = new QComboBox(this);
	this->metadata_box_->addItems(this->handler_.get_metadata_names());
	this->metadata_box_->adjustSize();
	this->metadata_box_->setFixedHeight(30);

	G_SET_BUTTON_ICON(this->metadata_view_button_, FILE_EYE_ICON_PNG, QSize(30, 30));

	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_right_layout_, row_layout, this->metadata_name_label_, this->metadata_box_, this->metadata_view_button_);
	connect(this->metadata_view_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_view_metadata);

	G_SET_LABEL_PRECISE(this->pathway_label_, "Pathway", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));
	this->left_top_right_layout_->addWidget(this->pathway_label_);

	G_SET_LABEL(this->database_label_, "Database", soap::MiddleSize);
	this->database_box_ = new QComboBox(this);
	G_SET_BUTTON(this->database_show_button_, "View", soap::MiddleSize);

	G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_right_layout_, row_layout, this->database_label_, this->database_box_, this->database_show_button_);
	connect(this->database_show_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_explore_pathway);

	G_SET_LABEL(this->pathway_name_label_, "Pathway", soap::MiddleSize);
	G_SET_LINEEDIT(this->pathway_line_edit_, "", soap::MiddleSize);
	G_SET_BUTTON_ICON(this->pathway_view_button_, FILE_EYE_ICON_PNG, QSize(30, 30));
	G_SET_BUTTON(this->pathway_alternative_button_, "►", QSize(30, 30));
	G_SET_BUTTON(this->pathway_active_button_, "▼", QSize(30, 30));

	G_ADD_FIVE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_right_layout_, row_layout, this->pathway_name_label_, this->pathway_line_edit_,
		this->pathway_view_button_, this->pathway_alternative_button_, this->pathway_active_button_);

	connect(this->pathway_view_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_view_pathway);
	connect(this->pathway_alternative_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_pathway_to_alternative);
	connect(this->pathway_active_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_pathway_to_active);


	G_SET_LABEL_PRECISE(this->data_source_label_, "Data Source", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));
	this->left_top_right_layout_->addWidget(this->data_source_label_);

	G_SET_LABEL(this->normalize_label_, "Normalize", soap::MiddleSize);
	G_SET_SWITCH(this->normalize_switch_, true, this->normalize_label_, soap::MiddleSize);

	G_ADD_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_right_layout_, row_layout, this->normalize_label_, this->normalize_switch_);

	if (this->handler_.type_ == FeatureHandler::DataType::SingleCellMultiome || this->handler_.type_ == FeatureHandler::DataType::SingleCellAtac) {
		G_SET_LABEL(this->gene_activity_label_, "Gene Activity", soap::MiddleSize);
		G_SET_SWITCH(this->gene_activity_switch_, false, this->gene_activity_label_, soap::MiddleSize);
		G_ADD_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(this->left_top_right_layout_, row_layout, this->gene_activity_label_, this->gene_activity_switch_);
	}

	G_SET_LABEL_PRECISE(this->alternative_label_, "Alternative Items", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));
	this->left_top_right_layout_->addWidget(this->alternative_label_);

	this->alternative_tree_widget_ = new CellMarkerTreeWidget(CellMarkerTreeWidget::MarkerWidgetType::AlternativeWidget, this);
	this->left_top_right_layout_->addWidget(this->alternative_tree_widget_);

	this->left_layout_->addLayout(this->left_top_right_layout_, 0, 1);

	connect(this->alternative_tree_widget_, &CellMarkerTreeWidget::x_show_gene, this, &EmbeddingViewWindow::s_show_gene);
	connect(this->alternative_tree_widget_, &CellMarkerTreeWidget::x_show_type, this, &EmbeddingViewWindow::s_show_type);
}

void EmbeddingViewWindow::s_pathway_to_alternative() {

	QString pathway_name = this->pathway_line_edit_->text();

	QStringList gene_list;

	bool found = false;

	for (const auto& database : this->pathway_content_) {
		if (database.contains(pathway_name)) {
			gene_list = database[pathway_name];
			found = true;
			break;
		}
	}

	if (!found) {
		G_LOG(pathway_name + " is not found in database!");
	}

	emit x_type_to_alternative(pathway_name, gene_list);
};

void EmbeddingViewWindow::s_pathway_to_active() {

	QString pathway_name = this->pathway_line_edit_->text();

	QStringList gene_list;

	bool found = false;

	for (const auto& database : this->pathway_content_) {
		if (database.contains(pathway_name)) {
			gene_list = database[pathway_name];
			found = true;
			break;
		}
	}

	if (!found) {
		G_LOG(pathway_name + " is not found in database!");
	}

	emit x_type_to_active(pathway_name, gene_list);
};

void EmbeddingViewWindow::set_left_bottom_right_layout() {
	this->left_bottom_right_layout_ = new QVBoxLayout;

	QHBoxLayout* row_layout = new QHBoxLayout;

	G_SET_LABEL_PRECISE(this->active_label_, "Active Items", soap::LargeSize, QColor("#20b2aa"), QFont("Arial", 20, QFont::Bold));

	G_SET_BUTTON(this->active_clear_button_, "Clear", soap::MiddleSize);
	G_SET_BUTTON(this->active_refresh_button_, "Refresh", soap::MiddleSize);

	row_layout->addWidget(this->active_label_);
	row_layout->addStretch();
	row_layout->addWidget(this->active_clear_button_);
	row_layout->addWidget(this->active_refresh_button_);
	this->left_bottom_right_layout_->addLayout(row_layout);

	this->active_tree_widget_ = new CellMarkerTreeWidget(CellMarkerTreeWidget::MarkerWidgetType::ActiveWidget, this);
	this->left_bottom_right_layout_->addWidget(this->active_tree_widget_);

	this->left_layout_->addLayout(this->left_bottom_right_layout_, 1, 1);

	connect(this->active_tree_widget_, &CellMarkerTreeWidget::x_show_gene, this, &EmbeddingViewWindow::s_show_gene);
	connect(this->active_tree_widget_, &CellMarkerTreeWidget::x_show_type, this, &EmbeddingViewWindow::s_show_type);

	connect(this->active_clear_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_clear_active_items);
	connect(this->active_refresh_button_, &QPushButton::clicked, this, &EmbeddingViewWindow::s_refresh_active_items);
}

void EmbeddingViewWindow::s_set_graph_settings() {

	GraphSettingDialog::set_graph_setting(this->draw_suite_->graph_settings_);

	this->graph_setting_switch_->set_status(this->draw_suite_->graph_settings_.active());
};

void EmbeddingViewWindow::s_refresh_active_items() {

	if (this->active_tree_widget_->types_.size() == 0)return;

	if (this->active_tree_widget_->types_.size() == 1) {
		show_type(this->active_tree_widget_->types_[0], this->active_tree_widget_->markers_[0]);
	}
	else {
		QStringList genelist;
		for (const auto& markers : this->active_tree_widget_->markers_) {
			genelist << markers;
		}
		show_type("Active Items", genelist);
	}
};

void EmbeddingViewWindow::s_save_png() {

	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit }
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "PNG(*.png)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->savePng(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void EmbeddingViewWindow::s_save_bmp() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit }
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "BMP(*.bmp)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->saveBmp(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void EmbeddingViewWindow::s_save_jpg() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit }
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "JPG(*.jpg)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->saveJpg(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void EmbeddingViewWindow::s_save_pdf_and_png() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit }
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set PNG Name", "", "PNG(*.PNG)");
	if (picture_name.isEmpty())return;

	QString pdf_name = QFileDialog::getSaveFileName(this, "Set PDF Name", "", "PDF(*.pdf)");
	if (pdf_name.isEmpty())return;

	bool success = this->draw_area_->savePdf(pdf_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (!success) {
		G_LOG("PDF saving failed");
		return;
	}
	else {
		G_LOG("PDF saved.");
	}
	std::string pdf_file_name = pdf_name.toStdString();
	std::string picture_file_name = picture_name.toStdString();

	if (!_Cs save_pdf_page_as_png(pdf_file_name, 0, picture_file_name)) {
		G_WARN("PNG Saving Failed.");
	}
	else {
		G_LOG("PNG saved.");
	}

};

void EmbeddingViewWindow::s_save_pdf() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"PDF Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit }
	);

	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set PDF Name", "", "PDF(*.pdf)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->savePdf(picture_name, picture_size[1].toInt(), picture_size[0].toInt());

	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void EmbeddingViewWindow::s_previous_plot() {
	if (this->plot_interacting_) {
		clear_interactive();
	}
	if (this->draw_suite_->current_plot_id_ <= 1)return;
	this->right_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[-- this->draw_suite_->current_plot_id_];
	this->right_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
};

void EmbeddingViewWindow::s_next_plot() {
	if (this->plot_interacting_) {
		clear_interactive();
	}
	if (this->draw_suite_->current_plot_id_ == this->draw_suite_->maximum_plot_id_)return;
	this->right_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[++ this->draw_suite_->current_plot_id_];
	this->right_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
};

void EmbeddingViewWindow::s_clear_plot() {
	if (this->plot_interacting_) {
		clear_interactive();
	}
	this->main_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);
	this->draw_area_ = nullptr;
	this->draw_suite_->clear();
};

void EmbeddingViewWindow::clear_interactive() {
	this->plot_interacting_ = false;

	this->point_x_.clear();
	this->point_y_.clear();
	for (auto p : this->plot_active_shape_.keys()) {
		p->removePlottable(this->plot_active_shape_[p]);
	}
	this->plot_active_shape_.clear();
	for (auto p : this->plot_to_graph_.keys()) {
		p->removeGraph(this->plot_to_graph_[p]);
	}
	this->plot_to_graph_.clear();
	for (auto p : this->plot_shapes_.keys()) {
		for (auto s : this->plot_shapes_[p]) {
			p->removePlottable(s);
		}
	}
	this->plot_shapes_.clear();
};

void EmbeddingViewWindow::s_pop_plot() {
	if (this->plot_interacting_) {
		clear_interactive();
	}
	PlotWindow::show_plot(this->draw_suite_, "Figure", this->draw_area_->width(), this->draw_area_->height(), this);
};

void EmbeddingViewWindow::s_new_plot() {
	if (this->draw_area_ == nullptr) {
		this->draw_area_ = this->draw_suite_->current_plot_;
		this->right_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();
	}
	else {
		this->right_layout_->removeWidget(this->draw_area_);
		this->draw_area_->setVisible(false);
		this->draw_area_ = this->draw_suite_->current_plot_;
		this->right_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();
	}

};

void EmbeddingViewWindow::s_refresh_plot() {
	if (this->active_tree_widget_->types_.isEmpty())return;
	if (this->active_tree_widget_->types_.size() == 1) {
		show_type(this->active_tree_widget_->types_[0], this->active_tree_widget_->markers_[0]);
	}
	else {
		show_type("Active Items", _Cs unroll(this->active_tree_widget_->markers_));
	}
};
