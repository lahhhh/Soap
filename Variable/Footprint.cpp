#include "Footprint.h"

#include "CustomPlot.h"

Footprint::Footprint(
	const PatternWeightMatrix& motif,
	const Eigen::MatrixXi& insertion_matrix,
	const QStringList& cell_names,
	const QStringList& base_position_names,
	const QVector<double>& expected_insertions,
	const GenomicRange& motif_location
) :
	motif_(motif),
	insertion_matrix_(DenseInt::DataType::Plain, insertion_matrix, cell_names, base_position_names),
	expected_insertions_(expected_insertions),
	motif_location_(motif_location)
{};

bool Footprint::is_empty() const {
	return this->expected_insertions_.isEmpty();
};

bool Footprint::draw(
	const QString& file_name,
	const QString& factor_name,
	const QStringList& levels,
	const QStringList& factors,
	const GraphSettings& gs,
	const int height,
	const int width
) const {
	const int nrow = this->insertion_matrix_.mat_.rows(), ncol = this->insertion_matrix_.mat_.cols();
	if (nrow != factors.size() || ncol <= 495) {
		return false;
	}
	Eigen::ArrayXd flank_mean = this->insertion_matrix_.mat_.block(0, 0, nrow, 50).rowwise().sum().cast<double>() + this->insertion_matrix_.mat_.block(0, ncol - 50, nrow, 50).rowwise().sum().cast<double>();
	flank_mean /= 100;
	double all_mean = flank_mean.mean();
	if (all_mean == 0) {
		return false;
	}
	for (int i = 0; i < nrow; ++i) {
		if (flank_mean[i] == 0) {
			flank_mean[i] = all_mean;
		}
	}
	Eigen::MatrixXd normalized = this->insertion_matrix_.mat_.cast<double>();
	normalized.array().colwise() /= flank_mean;

	QList<Eigen::ArrayXd> locations;
	for (const auto& factor : levels) {
		auto index = custom::match(factors, factor);

		Eigen::MatrixXd sub_matrix = normalized(index, Eigen::all);

		Eigen::ArrayXd means = sub_matrix.colwise().mean();

		locations << means - custom::cast<Eigen::ArrayX>(this->expected_insertions_);

	}

	auto [draw_area, left_layout, legend_layout] = custom_plot::prepare_lg_lg(gs);

	QCPAxisRect* observation_axisrect = new QCPAxisRect(draw_area), * expect_axisrect = new QCPAxisRect(draw_area);
	left_layout->addElement(0, 0, observation_axisrect);
	left_layout->addElement(1, 0, expect_axisrect);
	left_layout->setRowStretchFactors({4.0, 1.0});

	int ncolor = levels.size();
	auto colors = gs.palette(levels);

	QVector<double> x_axis = custom::cast<double>(this->insertion_matrix_.colnames_);

	for (int i = 0; i < ncolor; ++i) {
		custom_plot::patch::line(draw_area, observation_axisrect, x_axis, custom::cast<QVector>(locations[i]), colors[i], 2);
	}
	double min_ob = std::ranges::min(custom::sapply(locations, [](auto&& arr) {return arr.minCoeff(); }));
	double max_ob = std::ranges::max(custom::sapply(locations, [](auto&& arr) {return arr.maxCoeff(); }));

	auto [min_x, max_x] = std::ranges::minmax(x_axis);

	custom_plot::patch::set_range(observation_axisrect, QCPRange(min_x, max_x), QCPRange(1.1 * min_ob - 0.1 * max_ob, 1.1 * max_ob - 0.1 * min_ob));
	custom_plot::set_left_title(observation_axisrect, "Insertion Frequency", gs);

	QPen pen(Qt::black);
	pen.setWidth(3);

	custom_plot::patch::remove_bottom_axis(observation_axisrect);

	observation_axisrect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	observation_axisrect->axis(QCPAxis::atLeft)->setSubTickPen(Qt::NoPen);
	observation_axisrect->axis(QCPAxis::atLeft)->setBasePen(pen);

	auto [min_expect, max_expect] = std::ranges::minmax(this->expected_insertions_);
	custom_plot::patch::set_range(expect_axisrect, QCPRange(min_x, max_x), QCPRange(1.1 * min_expect - 0.1 * max_expect, 1.1 * max_expect - 0.1 * min_expect));
	custom_plot::patch::line(draw_area, expect_axisrect, x_axis, this->expected_insertions_, Qt::black, 2);

	QCPMarginGroup* edge_margin_group = new QCPMarginGroup(draw_area);
	observation_axisrect->setMarginGroup(QCP::msLeft | QCP::msRight, edge_margin_group);
	expect_axisrect->setMarginGroup(QCP::msLeft | QCP::msRight, edge_margin_group);

	custom_plot::patch::remove_grid(expect_axisrect);
	expect_axisrect->axis(QCPAxis::atBottom)->setBasePen(pen);
	expect_axisrect->axis(QCPAxis::atBottom)->setSubTickPen(Qt::NoPen);
	expect_axisrect->axis(QCPAxis::atLeft)->setBasePen(pen);
	expect_axisrect->axis(QCPAxis::atLeft)->setSubTickPen(Qt::NoPen);
	expect_axisrect->axis(QCPAxis::atLeft)->ticker()->setTickCount(3);

	custom_plot::set_left_title(expect_axisrect, "Expected", gs);
	custom_plot::set_bottom_title(expect_axisrect, "Relative Position (bp)", gs, true);

	custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, factor_name, gs);

	custom_plot::add_title(draw_area, this->motif_.motif_name_ + " Footprint", gs);

	bool success = draw_area->savePng(file_name, width, height);

	delete draw_area;

	return success;
};