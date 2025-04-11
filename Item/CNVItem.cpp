#include "CNVItem.h"

#include "MatrixWindow.h"
#include "CustomPlot.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"

void CNVItem::__s_update_interface() {

	this->setText(1, "CNV Detection");

	if (this->data()->data_type_ == CNV::DataType::InferCnv) {
		this->setText(2, "Infer CNV");
	}

	if (this->data()->data_type_ == CNV::DataType::SciCnv) {
		this->setText(2, "SciCNV");
	}
};

void CNVItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("View", s_view);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void CNVItem::draw_cnv(
	const Eigen::MatrixXd& cnv_mat,
	const std::vector<std::tuple<QString, int, int>>& chromosome_location,
	const std::vector<std::tuple<QString, int, int>>& metadata_location
) {
	auto&& gs = this->draw_suite_->graph_settings_;

	int nrow = cnv_mat.rows(), ncol = cnv_mat.cols();

	QCustomPlot* draw_area = custom_plot::initialize_plot(gs);
	SoapTextElement* title = new SoapTextElement(draw_area, gs.get_title("CNV"), gs.get_title_font());
	draw_area->plotLayout()->addElement(0, 0, title);

	QCPLayoutGrid* main_layout = new QCPLayoutGrid;
	draw_area->plotLayout()->addElement(1, 0, main_layout);

	QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	axis_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	axis_rect->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);
	QCPAxisRect* left_bottom_legend = new QCPAxisRect(draw_area, true);

	main_layout->addElement(0, 0, left_bottom_legend);
	main_layout->addElement(0, 1, axis_rect);

	auto cluster_names = custom::sapply(metadata_location, [](auto&& loc) {return std::get<0>(loc); });

	auto cluster_colors = gs.palette(cluster_names);
	int index = 0;
	for (const auto& [cluster, start, n] : metadata_location) {

		double end = start + n;

		custom_plot::patch::rectangle_borderless(
			draw_area, left_bottom_legend, 0, start, 1, end - start, cluster_colors[index++]
		);

		QCPItemText* label_x = new QCPItemText(draw_area);
		label_x->setClipToAxisRect(false);
		label_x->position->setAxisRect(left_bottom_legend);
		label_x->position->setAxes(left_bottom_legend->axis(QCPAxis::atBottom), left_bottom_legend->axis(QCPAxis::atLeft));
		label_x->position->setType(QCPItemPosition::ptPlotCoords);
		label_x->setPositionAlignment(Qt::AlignRight | Qt::AlignVCenter);
		label_x->position->setCoords(-1, (start + end) / 2);
		label_x->setText(cluster);
		label_x->setFont(gs.get_left_label_font());
	}

	custom_plot::patch::remove_left_bottom_axis(left_bottom_legend);
	left_bottom_legend->setMinimumSize(custom_plot::utility::get_max_text_width(custom::sapply(metadata_location, [](auto&& t) {return std::get<0>(t); }), gs.get_left_label_font()) * 1.2, 200);

	left_bottom_legend->axis(QCPAxis::atBottom)->setRange(-11, 1);
	left_bottom_legend->axis(QCPAxis::atLeft)->setRange(0, ncol);


	QCPColorMap* heatmap = new QCPColorMap(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
	heatmap->data()->setSize(nrow, ncol);
	heatmap->data()->setRange(QCPRange(0, nrow - 1), QCPRange(0, ncol - 1));
	custom_plot::patch::remove_left_bottom_axis(axis_rect);
	custom_plot::patch::set_range(axis_rect, QCPRange(-0.5, nrow - 0.5), QCPRange(-0.5, ncol - 0.5));

	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j) {
			heatmap->data()->setCell(i, j, cnv_mat(i, j));
		}
	}
	heatmap->setInterpolate(false);
	heatmap->setTightBoundary(false);

	double min_val = cnv_mat.minCoeff(), max_val = cnv_mat.maxCoeff();
	double span = std::max(std::abs(1.0 - min_val), std::abs(1.0 - max_val));

	QCPColorGradient gradient;
	gradient.setColorStopAt(0.5, QColor(255, 255, 255));
	gradient.setColorStopAt(1.0, QColor(205, 38, 38));
	gradient.setColorStopAt(0.0, QColor(0, 0, 128));

	heatmap->setGradient(gradient);
	if (this->data()->data_type_ == CNV::DataType::InferCnv) {

		if (span > 0.5) {
			span = 0.5;
		}

		heatmap->setDataRange({ 1.0 - 0.5 * span, 1.0 + 0.5 * span });
	}
	else {
		heatmap->setDataRange({ -1.0, 1.0});
	}

	QCPColorScale* color_scale = new QCPColorScale(draw_area);
	color_scale->setType(QCPAxis::atRight);

	if (this->data()->data_type_ == CNV::DataType::InferCnv) {
		color_scale->setDataRange({ 1.0 - 0.5 * span, 1.0 + 0.6 * span });
	}
	else {
		color_scale->setDataRange({ -1.0, 1.0 });
	}
	color_scale->setGradient(gradient);

	color_scale->setMargins(QMargins(0, 0, 20, 0));

	color_scale->axis()->setTickPen(Qt::NoPen);
	color_scale->axis()->setSubTickPen(Qt::NoPen);
	color_scale->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);

	color_scale->axis()->setTickLabelFont(QFont("Arial", 15, QFont::Bold));

	color_scale->mAxisRect->axis(QCPAxis::atBottom)->setBasePen(Qt::NoPen);
	color_scale->mAxisRect->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
	color_scale->mAxisRect->axis(QCPAxis::atTop)->setBasePen(Qt::NoPen);
	color_scale->mAxisRect->axis(QCPAxis::atRight)->setBasePen(Qt::NoPen);
	main_layout->addElement(0, 2, color_scale);

	QCPAxisRect* bottom_legend = new QCPAxisRect(draw_area, true);
	QCPAxisRect* bottom_left = new QCPAxisRect(draw_area, false);
	QCPAxisRect* bottom_right = new QCPAxisRect(draw_area, false);
	main_layout->addElement(1, 0, bottom_left);
	main_layout->addElement(1, 1, bottom_legend);
	main_layout->addElement(1, 2, bottom_right);

	auto chr_colors = custom_plot::utility::kmeans_palette(chromosome_location.size());
	index = 0;
	for (const auto& [chr, start, n] : chromosome_location) {

		double end = start + n;
		draw_area->addGraph(bottom_legend->axis(QCPAxis::atBottom), bottom_legend->axis(QCPAxis::atLeft));
		QCPCurve* legend = new QCPCurve(bottom_legend->axis(QCPAxis::atBottom), bottom_legend->axis(QCPAxis::atLeft));
		QVector<QCPCurveData> data(5);
		data[0] = QCPCurveData(0, start, 0);
		data[1] = QCPCurveData(1, start, 1);
		data[2] = QCPCurveData(2, end, 1);
		data[3] = QCPCurveData(3, end, 0);
		data[4] = QCPCurveData(0, start, 0);
		legend->setPen(Qt::NoPen);
		legend->setBrush(QBrush(chr_colors[index++]));
		legend->data()->set(data, true);

		QCPItemText* label_x = new QCPItemText(draw_area);
		label_x->setClipToAxisRect(false);
		label_x->position->setAxisRect(bottom_legend);
		label_x->position->setAxes(bottom_legend->axis(QCPAxis::atBottom), bottom_legend->axis(QCPAxis::atLeft));
		label_x->position->setType(QCPItemPosition::ptPlotCoords);
		label_x->setPositionAlignment(Qt::AlignTop | Qt::AlignHCenter);
		label_x->position->setCoords((start + end) / 2, -(index % 2) * 1.5);
		label_x->setText(chr);
		label_x->setFont(gs.get_bottom_label_font());
	}

	bottom_legend->axis(QCPAxis::atBottom)->grid()->setVisible(false);
	bottom_legend->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	bottom_legend->axis(QCPAxis::atBottom)->setBasePen(Qt::NoPen);
	bottom_legend->axis(QCPAxis::atBottom)->setTicks(false);
	bottom_legend->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
	bottom_legend->axis(QCPAxis::atLeft)->setTicks(false);

	bottom_legend->axis(QCPAxis::atBottom)->setRange(0, nrow);
	bottom_legend->axis(QCPAxis::atLeft)->setRange(-5, 1);

	bottom_legend->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	bottom_legend->setMargins({ 0, 0, 0, 0 });

	main_layout->setColumnStretchFactor(0, 1);
	main_layout->setColumnStretchFactor(1, 6);
	main_layout->setColumnStretchFactor(2, 1);

	main_layout->setRowStretchFactor(0, 7);
	main_layout->setRowStretchFactor(1, 1);

	draw_area->plotLayout()->setRowStretchFactor(0, 1);
	draw_area->plotLayout()->setRowStretchFactor(1, 5);

	this->draw_suite_->update(draw_area);
};

void CNVItem::s_view() {

	this->draw_cnv(this->data()->mat_, this->data()->chromosome_info_, this->data()->cluster_info_);

};

