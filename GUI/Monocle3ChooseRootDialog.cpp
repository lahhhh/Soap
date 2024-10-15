#include "Monocle3ChooseRootDialog.h"

#include "CustomPlot.h"

#include "SoapGUI.h"

Monocle3ChooseRootDialog::Monocle3ChooseRootDialog(Monocle3* data):
	data_(data)
{

	QVBoxLayout* main_layout = new QVBoxLayout;

	this->draw_area_ = new QCustomPlot();
	this->draw_area_->plotLayout()->clear();

	this->axis_rect_ = new QCPAxisRect(this->draw_area_);
	this->draw_area_->plotLayout()->addElement(0, 0, this->axis_rect_);

	Eigen::ArrayXd x = this->data_->cell_embedding_.row(0);
	Eigen::ArrayXd y = this->data_->cell_embedding_.row(1);

	custom_plot::patch::set_range(this->axis_rect_, custom_plot::utility::get_range(x, y));

	custom_plot::patch::scatter(
		this->draw_area_,
		this->axis_rect_,
		x,
		y,
		QColor(221, 221, 221, 128),
		5
	);

	x = this->data_->pr_embedding_.row(0);
	y = this->data_->pr_embedding_.row(1);

	custom_plot::patch::scatter(
		this->draw_area_,
		this->axis_rect_,
		x,
		y,
		Qt::black,
		4
	);

	custom_plot::patch::remove_left_bottom_axis(this->axis_rect_);

	igraph_vector_int_t edge_list;
	igraph_vector_int_init(&edge_list, 0);
	igraph_get_edgelist(&this->data_->pr_graph_, &edge_list, 0);

	int n_edge = igraph_vector_int_size(&edge_list) / 2;

	for (int i = 0; i < n_edge; ++i) {
		int v1 = VECTOR(edge_list)[i * 2];
		int v2 = VECTOR(edge_list)[i * 2 + 1];

		custom_plot::patch::line(this->draw_area_, this->axis_rect_,
			this->data_->pr_embedding_(Eigen::all, { v1, v2 }).row(0),
			this->data_->pr_embedding_(Eigen::all, { v1, v2 }).row(1),
			Qt::black, 2);
	}
	igraph_vector_int_destroy(&edge_list);

	connect(this->draw_area_, &QCustomPlot::mousePress, this, &Monocle3ChooseRootDialog::s_plot_mouse_press);
	connect(this->draw_area_, &QCustomPlot::mouseMove, this, &Monocle3ChooseRootDialog::s_plot_mouse_move);

	main_layout->addWidget(this->draw_area_);	
	this->draw_area_->replot();

	G_SET_FINISH_BUTTON;
	G_SET_CANCEL_BUTTON;

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT(main_layout, row_layout, this->finish_button_, this->cancel_button_);

	this->setLayout(main_layout);
	this->setFixedSize(800, 800);

	connect(this->finish_button_, &QPushButton::clicked, this, &Monocle3ChooseRootDialog::accept);
	connect(this->cancel_button_, &QPushButton::clicked, this, &Monocle3ChooseRootDialog::reject);

	G_SET_ICON;
	this->setWindowTitle("Choose Root Nodes");
	this->exec();
};

void Monocle3ChooseRootDialog::s_plot_mouse_press(QMouseEvent* event) {

	int x_pos = event->pos().x(), y_pos = event->pos().y();

	double x = this->axis_rect_->axis(QCPAxis::atBottom)->pixelToCoord(x_pos), y = this->axis_rect_->axis(QCPAxis::atLeft)->pixelToCoord(y_pos);

	double min_dist{ 0.0 }, xx{ 0.0 }, yy{ 0.0 };
	int index{ -1 };
	int n_pr_node = this->data_->pr_embedding_.cols();
	auto& emb = this->data_->pr_embedding_;

	for (int i = 0; i < n_pr_node; ++i) {
		double dist = (emb(0, i) - x) * (emb(0, i) - x) + (emb(1, i) - y) * (emb(1, i) - y);

		if (min_dist == 0.0 || dist < min_dist) {
			index = i;
			min_dist = dist;
			xx = emb(0, i);
			yy = emb(1, i);
		}
	}

	if (this->choosed_graph_.contains(index)) {

		this->draw_area_->removeGraph(this->choosed_graph_[index]);
	}
	else {
		this->choosed_graph_[index] = custom_plot::patch::scatter(
			this->draw_area_,
			this->axis_rect_,
			QVector<double>{xx},
			QVector<double>{yy},
			Qt::red,
			6
		);
	}

	this->draw_area_->replot();
};

void Monocle3ChooseRootDialog::s_plot_mouse_move(QMouseEvent* event) {

	int x_pos = event->pos().x(), y_pos = event->pos().y();

	double x = this->axis_rect_->axis(QCPAxis::atBottom)->pixelToCoord(x_pos), y = this->axis_rect_->axis(QCPAxis::atLeft)->pixelToCoord(y_pos);

	double min_dist{ 0.0 }, xx{ 0.0 }, yy{ 0.0 };
	int n_pr_node = this->data_->pr_embedding_.cols();
	auto& emb = this->data_->pr_embedding_;

	for (int i = 0; i < n_pr_node; ++i) {
		double dist = (emb(0, i) - x) * (emb(0, i) - x) + (emb(1, i) - y) * (emb(1, i) - y);

		if (min_dist == 0.0 || dist < min_dist) {
			min_dist = dist;
			xx = emb(0, i);
			yy = emb(1, i);
		}
	}

	if (this->move_graph_ != nullptr) {
		this->draw_area_->removeGraph(this->move_graph_);
	}

	this->move_graph_ = custom_plot::patch::scatter(
		this->draw_area_,
		this->axis_rect_,
		QVector<double>{xx},
		QVector<double>{yy},
		Qt::red,
		6
	);

	this->draw_area_->replot();

};

void Monocle3ChooseRootDialog::accept() {
	
	QDialog::accept();

	this->is_accepted_ = true;
}

void Monocle3ChooseRootDialog::reject() {

	QDialog::reject();
}

QVector<int> Monocle3ChooseRootDialog::get_response(Monocle3* data) {

	Monocle3ChooseRootDialog dlg(data);

	if (dlg.is_accepted_) {
		return custom::keys(dlg.choosed_graph_);
	}
	else return {};
};