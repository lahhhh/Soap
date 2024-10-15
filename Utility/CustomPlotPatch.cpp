#include "CustomPlotPatch.h"

#include "Custom.h"
#include "CustomPlotUtility.h"

namespace custom_plot {
	
	namespace patch {

		QCPGraph* shape(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor border_color,
			QColor fill_color,
			int width
		) {

			auto graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			
			QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			
			int length = x.size();
			
			QVector<QCPCurveData> data(length + 1);
			for (int j = 0; j < length; ++j) {
				data[j] = QCPCurveData(j, x[j], y[j]);
			}
			data[length] = QCPCurveData(length, x[0], y[0]);

			QPen pen(border_color);
			pen.setWidth(width);
			
			shape->setPen(pen);
			shape->setBrush(QBrush(fill_color));
			shape->data()->set(data);

			return graph;
		};

		QCPGraph* shape(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			QColor border_color,
			QColor fill_color,
			int width
		) {
			return shape(
				draw_area,
				axis_rect,
				custom::cast<QVector>(x),
				custom::cast<QVector>(y),
				border_color,
				fill_color,
				width
			);
		};

		QCPGraph* shape_borderless(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor fill_color
		) {
		
			auto graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));

			QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));

			int length = x.size();

			QVector<QCPCurveData> data(length + 1);
			for (int j = 0; j < length; ++j) {
				data[j] = QCPCurveData(j, x[j], y[j]);
			}
			data[length] = QCPCurveData(length, x[0], y[0]);

			shape->setPen(Qt::NoPen);
			shape->setBrush(QBrush(fill_color));
			shape->data()->set(data);

			return graph;
		};

		QCPGraph* shape_borderless(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			QColor fill_color
		) {
			return shape_borderless(
				draw_area,
				axis_rect,
				custom::cast<QVector>(x),
				custom::cast<QVector>(y),
				fill_color
			);
		};

		QCPGraph* rectangle(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			double start_left_x,
			double start_bottom_y,
			double width,
			double height,
			QColor color,
			int border_width,
			QColor border_color
		) {
			return shape(
				draw_area,
				axis_rect,
				QVector<double>{start_left_x, start_left_x, start_left_x + width, start_left_x + width},
				QVector<double>{start_bottom_y, start_bottom_y + height, start_bottom_y + height, start_bottom_y},
				border_color,
				color,
				border_width
			);
		};

		QCPGraph* rectangle_borderless(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			double start_left_x,
			double start_bottom_y,
			double width,
			double height,
			QColor color
		) {
			return shape_borderless(
				draw_area,
				axis_rect,
				QVector<double>{start_left_x, start_left_x, start_left_x + width, start_left_x + width},
				QVector<double>{start_bottom_y, start_bottom_y + height, start_bottom_y + height, start_bottom_y},
				color
			);
		};

		QCPGraph* scatter(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			const QColor& color,
			int scatter_size) 
		{
			QPen pen(color);

			QCPGraph* graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			graph->setLineStyle(QCPGraph::lsNone);
			graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, scatter_size));
			graph->setPen(pen);
			graph->setData(x, y);

			return graph;
		};

		QCPGraph* scatter(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			const QColor& color,
			int scatter_size
		) {
		
			return scatter(
				draw_area,
				axis_rect,
				custom::cast<QVector>(x),
				custom::cast<QVector>(y),
				color,
				scatter_size
			);
		};

		void scatter_category(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			const QStringList& values,
			const QStringList& levels,
			const QList<QColor>& colors,
			int scatter_size
		) {
			const int n_level = levels.size();

			for (int i = 0; i < n_level; ++i) {

				auto filter = custom::equal(values, levels[i]);
				QVector<double> sub_x = custom::sliced(x, filter);
				QVector<double> sub_y = custom::sliced(y, filter);

				scatter(
					draw_area,
					axis_rect,
					sub_x,
					sub_y,
					colors[i],
					scatter_size
				);
			}
		};

		void scatter_category(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			const QStringList& values,
			const QStringList& levels,
			const QList<QColor>& colors,
			int scatter_size
		) {
			
			scatter_category(
				draw_area,
				axis_rect,
				custom::cast<QVector>(x),
				custom::cast<QVector>(y),
				values,
				levels,
				colors,
				scatter_size
			);
		};

		void scatter_gradient(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			const Eigen::ArrayXd& value,
			double min,
			double max,
			QColor low_color,
			QColor middle_color,
			QColor high_color,
			int scatter_size
		) {
			if (min >= max) {
				
				scatter(
					draw_area,
					axis_rect,
					x,
					y,
					low_color,
					scatter_size
				);

				return;
			}

			int low_red = low_color.red(), low_green = low_color.green(), low_blue = low_color.blue();
			int high_red = high_color.red(), high_green = high_color.green(), high_blue = high_color.blue();
			int middle_red = middle_color.red(), middle_green = middle_color.green(), middle_blue = middle_color.blue();
			int r_space_lower = middle_red - low_red, g_space_lower = middle_green - low_green, b_space_lower = middle_blue - low_blue;
			int r_space_upper = high_red - middle_red, g_space_upper = high_green - middle_green, b_space_upper = high_blue - middle_blue;

			double span = max - min;
			double middle = (max + min) / 2;
			double space = span / 99.0;
			min -= space / 2.0;
			max += space / 2.0;

			QVector<QVector<double>> xs(100), ys(100);
			int size = x.size();
			for (int i = 0; i < size; ++i) {

				int loc = std::floor((value[i] - min) / space);

				if (loc < 0) {
					xs[0] << x[i];
					ys[0] << y[i];
				}
				else if (loc > 99) {
					xs[99] << x[i];
					ys[99] << y[i];
				}
				else {
					xs[loc] << x[i];
					ys[loc] << y[i];
				}
			}

			for (int i = 0; i < 100; ++i) {

				if (xs[i].isEmpty())
				{
					continue;
				}

				double val = min + (i + 0.5) * space;

				int r, g, b;
				if (val < middle) {
					r = low_red + (int)(((double)(val - min) / span) * r_space_lower);
					b = low_blue + (int)(((double)(val - min) / span) * b_space_lower);
					g = low_green + (int)(((double)(val - min) / span) * g_space_lower);
				}
				else {
					r = middle_red + (int)(((double)(val - middle) / span) * r_space_upper);
					b = middle_blue + (int)(((double)(val - middle) / span) * b_space_upper);
					g = middle_green + (int)(((double)(val - middle) / span) * g_space_upper);
				}

				scatter(
					draw_area,
					axis_rect,
					xs[i],
					ys[i],
					QColor(r, g, b),
					scatter_size
				);
			}
		}

		QCPGraph* curve(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor color,
			int width,
			Qt::PenStyle style
		) {
			QPen pen(color);
			pen.setWidth(width);
			pen.setStyle(style);

			auto graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			QCPCurve* shape = new QCPCurve(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));

			const int n_point = x.size();
			QVector<QCPCurveData> data(n_point);

			for (int i = 0; i < n_point; ++i) {
				data[i] = QCPCurveData(0, x[i], y[i]);
			}

			shape->setPen(pen);
			shape->data()->set(data);

			return graph;
		};

		QCPGraph* line(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& x,
			const QVector<double>& y,
			QColor color,
			int width,
			Qt::PenStyle style
		) {

			QPen pen(color);
			pen.setWidth(width);
			pen.setStyle(style);

			auto graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			graph->setPen(pen);
			graph->setData(x, y);

			return graph;
		};

		QCPGraph* line(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			QColor color,
			int width,
			Qt::PenStyle style
		) {
		
			return line(
				draw_area,
				axis_rect,
				custom::cast<QVector>(x),
				custom::cast<QVector>(y),
				color,
				width,
				style
			);
		};

		QCPGraph* bar(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const Eigen::ArrayXd& bar_location,
			const Eigen::ArrayXd& bar_length,
			int bar_width,
			bool vertical
		) {
			QPen pen(color);
			pen.setWidth(bar_width);

			QCPGraph* graph;

			if (vertical) {
				graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			}
			else {
				graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atLeft), axis_rect->axis(QCPAxis::atBottom));
			}

			graph->setLineStyle(QCPGraph::lsImpulse);
			graph->setPen(pen);
			graph->setData(custom::cast<QVector>(bar_location), custom::cast<QVector>(bar_length));

			return graph;
		};

		QCPGraph* bar(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const QVector<double>& bar_location,
			const QVector<double>& bar_length,
			int bar_width,
			bool vertical
		) {
			QPen pen(color);
			pen.setWidth(bar_width);

			QCPGraph* graph;

			if (vertical) {
				graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			}
			else {
				graph = draw_area->addGraph(axis_rect->axis(QCPAxis::atLeft), axis_rect->axis(QCPAxis::atBottom));
			}

			graph->setLineStyle(QCPGraph::lsImpulse);
			graph->setPen(pen);
			graph->setData(bar_location, bar_length);

			return graph;
		};

		void bar_polychrome(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& bar_location,
			const Eigen::ArrayXd& bar_length,
			const QList<QColor>& colors,
			int bar_width,
			bool vertical)
		{
			const int n_bar = bar_location.size();

			for (int i = 0; i < n_bar; ++i) {

				QPen pen(colors[i]);
				pen.setWidth(bar_width);

				if (vertical) {
					draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
				}
				else {
					draw_area->addGraph(axis_rect->axis(QCPAxis::atLeft), axis_rect->axis(QCPAxis::atBottom));
				}

				draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
				draw_area->graph()->setPen(pen);
				draw_area->graph()->setData(QVector<double>{bar_location[i]}, QVector<double>{bar_length[i]});
			}
		}

		void bar_gradient(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& bar_location,
			const Eigen::ArrayXd& bar_length,
			const Eigen::ArrayXd& values,
			QColor low_color,
			QColor middle_color,
			QColor high_color,
			int bar_width,
			bool vertical
		) {

			double min = values.minCoeff();
			double max = values.maxCoeff();
			double span = (max - min) / 2;
			double middle = (min + max) / 2;

			int bar_number = values.size();

			for (int i = 0; i < bar_number; ++i) {

				QColor color;
				double val = values[i];
				if (val < middle) {
					color = custom_plot::utility::gradient_color((val - min) / span, low_color, middle_color);
				}
				else {
					color = custom_plot::utility::gradient_color((val - middle) / span, middle_color, high_color);
				}

				bar(
					draw_area,
					axis_rect,
					color,
					QVector<double>{ bar_location[i] },
					QVector<double>{ bar_length[i] },
					bar_width,
					vertical
				);
			}
		}

		void bar_single_stack(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QMap<QString, int>& distribution,
			const QStringList& levels,
			const QVector<QColor>& colors,
			double start_point_x,
			double width
		) {

			double sum = custom::sum(distribution.values());
			int size = levels.size();
			double height{ 0.0 };

			for (int i = 0; i < size; ++i) {

				double sub_height = distribution[levels[i]] / sum;				

				custom_plot::patch::shape_borderless(
					draw_area,
					axis_rect,
					QVector<double>{start_point_x - width, start_point_x - width, start_point_x + width, start_point_x + width},
					QVector<double>{height, height + sub_height, height + sub_height, height },
					colors[i]
				);

				height += sub_height;
			}
		};

		void bar_stack(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QStringList& data,
			const QStringList& data_levels,
			const QVector<QColor>& colors,
			double start_x,
			double width,
			double margin
		) {
			int n_level = levels.size();

			for (int i = 0; i < n_level; ++i) {

				auto filter = custom::equal(group, levels[i]);
				auto subdata = custom::sliced(data, filter);
				if (subdata.size() == 0) {
					continue;
				}

				custom_plot::patch::bar_single_stack(draw_area, axis_rect, custom::table(subdata), data_levels, colors, start_x + i * margin, width);
			}
		};

		void add_title(QCustomPlot* draw_area, QCPLayoutGrid* layout, const QString& title, QFont font) {

			layout->insertRow(0);

			SoapTextElement* text = new SoapTextElement(draw_area, title, font);
			layout->addElement(0, 0, text);
		};

		void add_label(
			QCustomPlot* draw_area,
			QCPAxisRect* rect,
			const QString& label,
			double x,
			double y,
			QFont font,
			Qt::Alignment alignment,
			double degree
		) {

			QCPItemText* item = new QCPItemText(draw_area);
			
			item->position->setAxisRect(rect);
			item->position->setAxes(rect->axis(QCPAxis::atBottom), rect->axis(QCPAxis::atLeft));
			item->position->setType(QCPItemPosition::ptPlotCoords);
			item->position->setCoords(x, y);

			item->setPositionAlignment(alignment);
			item->setClipToAxisRect(false);
			item->setText(label);
			item->setFont(font);
			item->setRotation(degree);
		};

		void dot(QCustomPlot* draw_area, QCPAxisRect* axis_rect, double x, double y, int size, QColor color) {

			draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
			draw_area->graph()->setLineStyle(QCPGraph::lsNone);
			draw_area->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, size));

			QVector<double> sub_x{ x }, sub_y{ y };
			QPen pen(color);
			draw_area->graph()->setPen(pen);
			draw_area->graph()->setData(sub_x, sub_y);
		};

		void set_proportion_left_axis(QCPAxisRect* axis_rect, QFont tick_label_font) {

			QSharedPointer<QCPAxisTickerText> ticker(new QCPAxisTickerText);
			ticker->setTicks(
				QVector<double>{ 0.0, 0.25, 0.5, 0.75, 1}, 
				{ "0%", "25%", "50%", "75%", "100%" });

			QPen pen(Qt::black);
			pen.setWidth(3);
			
			axis_rect->axis(QCPAxis::atLeft)->setTicks(true);
			axis_rect->axis(QCPAxis::atLeft)->setTickLabels(true);
			axis_rect->axis(QCPAxis::atLeft)->setTicker(ticker);
			axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 5);
			axis_rect->axis(QCPAxis::atLeft)->setBasePen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setTickPen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setSubTickPen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setTickLabelFont(tick_label_font);
		};

		void remove_axis(QCPAxisRect* axis_rect, QCPAxis::AxisType type) {

			axis_rect->axis(type)->setBasePen(Qt::NoPen);
			axis_rect->axis(type)->setTicks(false);
			axis_rect->axis(type)->grid()->setVisible(false);
		};

		void clear_axis(QCPAxisRect* axis_rect, QCPAxis::AxisType type) {

			axis_rect->axis(type)->setTicks(false);
			axis_rect->axis(type)->grid()->setVisible(false);
		};

		void remove_left_axis(QCPAxisRect* axis_rect) {
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atLeft);
		};

		void remove_bottom_axis(QCPAxisRect* axis_rect) {
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atBottom);
		};

		void clear_left_axis(QCPAxisRect* axis_rect) {
			custom_plot::patch::clear_axis(axis_rect, QCPAxis::atLeft);
		};

		void clear_bottom_axis(QCPAxisRect* axis_rect) {
			custom_plot::patch::clear_axis(axis_rect, QCPAxis::atBottom);
		};

		void remove_all_axis(QCPAxisRect* axis_rect) {

			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atBottom);
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atTop);
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atRight);
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atLeft);
		};

		void remove_left_bottom_axis(QCPAxisRect* axis_rect) {

			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atBottom);
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atLeft);
		};

		void remove_left_ticks(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atLeft)->setTicks(false);
		};

		void remove_bottom_ticks(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atBottom)->setTicks(false);
		};

		void remove_left_subticks(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atLeft)->setSubTicks(false);
		};

		void remove_bottom_subticks(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atBottom)->setSubTicks(false);
		};

		void remove_left_grid(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
		};

		void remove_bottom_grid(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
		};

		void remove_grid(QCPAxisRect* axis_rect) {

			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
		};

		void set_range(QCPAxisRect* axis_rect, QCPRange bottom_range, QCPRange left_range) {

			axis_rect->axis(QCPAxis::atBottom)->setRange(bottom_range);
			axis_rect->axis(QCPAxis::atLeft)->setRange(left_range);
		};

		void set_range(QCPAxisRect* axis_rect, std::pair<QCPRange, QCPRange> range) {

			axis_rect->axis(QCPAxis::atBottom)->setRange(range.first);
			axis_rect->axis(QCPAxis::atLeft)->setRange(range.second);
		};

		void set_fixed_size(QCPAxisRect* axis_rect, int width, int height) {

			axis_rect->setMaximumSize(width, height);
			axis_rect->setMinimumSize(width, height);
		};

		void set_single_square_legend(
			QCustomPlot* draw_area,
			QCPAxisRect* legend_rect,
			const QString& label,
			QColor color,
			QFont font,
			double start_point_x,
			double start_point_y,
			int legend_size
		) {

			patch::shape_borderless(
				draw_area,
				legend_rect,
				QVector<double>() << -legend_size << -legend_size << legend_size << legend_size,
				QVector<double>() << start_point_y - legend_size << start_point_y + legend_size << start_point_y + legend_size << start_point_y - legend_size,
				color
			);

			custom_plot::patch::add_label(
				draw_area,
				legend_rect,
				label,
				start_point_x,
				start_point_y,
				font,
				Qt::AlignVCenter | Qt::AlignLeft
			);
		};

		void set_single_round_legend(
			QCustomPlot* draw_area,
			QCPAxisRect* legend_rect,
			const QString& label,
			QColor color,
			QFont font,
			double start_point_x,
			double start_point_y,
			int legend_size
		) {
			QVector<double> x{ 0 }, y{ start_point_y };
			QPen pen(color);

			draw_area->addGraph(legend_rect->axis(QCPAxis::atBottom), legend_rect->axis(QCPAxis::atLeft));
			draw_area->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, legend_size));
			draw_area->graph()->setPen(pen);
			draw_area->graph()->setData(x, y);

			custom_plot::patch::add_label(
				draw_area,
				legend_rect,
				label,
				start_point_x,
				start_point_y,
				font,
				Qt::AlignVCenter | Qt::AlignLeft
			);
		};

		void set_border_only(
			QCPAxisRect* axis_rect,
			QColor color,
			int width
		) {
			QPen pen(color);
			pen.setWidth(width);

			axis_rect->setupFullAxesBox();

			axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
			axis_rect->axis(QCPAxis::atBottom)->setTicks(false);
			axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);

			axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			axis_rect->axis(QCPAxis::atLeft)->setTicks(false);
			axis_rect->axis(QCPAxis::atLeft)->setBasePen(pen);

			axis_rect->axis(QCPAxis::atTop)->grid()->setVisible(false);
			axis_rect->axis(QCPAxis::atTop)->setTicks(false);
			axis_rect->axis(QCPAxis::atTop)->setBasePen(pen);

			axis_rect->axis(QCPAxis::atRight)->grid()->setVisible(false);
			axis_rect->axis(QCPAxis::atRight)->setTicks(false);
			axis_rect->axis(QCPAxis::atRight)->setBasePen(pen);
		};

		QCPGraph* sector(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			double x,
			double y,
			double start_angle,
			double stop_angle,
			double radius,
			int segments
		)
		{
			QVector<double> point_x(segments + 3), point_y(segments + 3);
			point_x[0] = point_x[segments + 2] = x;
			point_y[0] = point_y[segments + 2] = y;
			double angle_segment = (stop_angle - start_angle) / segments;

			for (int i = 0; i < segments + 1; ++i) {
				point_x[i + 1] = x + std::cos(start_angle + angle_segment * i) * radius;
				point_y[i + 1] = y + std::sin(start_angle + angle_segment * i) * radius;
			}

			return custom_plot::patch::shape_borderless(
				draw_area,
				axis_rect,
				point_x,
				point_y,
				color
			);
		};

		void pie(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QMap<QString, int>& distribution,
			const QStringList& levels,
			const QVector<QColor>& colors,
			double start_point_x,
			double start_point_y,
			double radius,
			bool white_border
		) {
			double sum = custom::sum(distribution.values());
			constexpr int control_point_number = 360;
			int size = levels.size();
			double angle = 0.0;
			std::vector<std::pair<double, double>> end_points;
			
			for (int i = 0; i < size; ++i) {
				QString level = levels[i];

				double start_angle = angle;
				angle += (distribution[level] / sum * 2 * M_PI);
				double end_angle = angle;

				if (white_border) {
					end_points.emplace_back(start_point_x + std::cos(start_angle) * radius,
						start_point_y + std::sin(start_angle) * radius);
				}

				 custom_plot::patch::sector(
					 draw_area,
					 axis_rect,
					 colors[i],
					 start_point_x, 
					 start_point_y, 
					 start_angle, 
					 end_angle, 
					 radius, 
					 control_point_number
				 );
			}

			if (white_border) {

				for (int i = 0; i < size; ++i) {

					custom_plot::patch::line(
						draw_area,
						axis_rect,
						QVector<double>{start_point_x, end_points[i].first},
						QVector<double>{start_point_y, end_points[i].second},
						Qt::white,
						2
					);
				}
			}
		};

		static QVector<double> fivenum(const QVector<double>& value) {

			QVector<double> x = custom::sorted(value);

			auto n = x.size();

			QVector<double> res(5);
			res[0] = x[0];
			double n4 = std::floor(double(n + 3) / 2) / 2;
			res[1] = (x[qsizetype(std::floor(n4) - 1)] + x[qsizetype(std::ceil(n4) - 1)]) / 2;
			res[2] = (x[qsizetype(std::floor(double(n + 1) / 2) - 1)] + x[qsizetype(std::ceil(double(n + 1) / 2) - 1)]) / 2;
			res[3] = (x[qsizetype(std::floor(n + 1 - n4) - 1)] + x[qsizetype(std::ceil(n + 1 - n4) - 1)]) / 2;
			res[4] = x[n - 1];

			return res;
		};

		void box(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const QVector<double>& value,
			double center,
			bool draw_outlier,
			int outlier_scatter_size
		) {
			QVector<double> stats = fivenum(value);
			constexpr double coef{ 1.5 }; // should be a parameter
			double iqr = stats[3] - stats[1];
			auto outs = custom::sliced(value, custom::greater_than(value, stats[3] + coef * iqr) +
				custom::less_than(value, stats[1] - coef * iqr));
			QVector<double> filtered = value;

			if (!outs.isEmpty()) {
				filtered = custom::sliced(value, custom::less_equal(value, stats[3] + coef * iqr) *
					custom::greater_equal(value, stats[1] - coef * iqr));
				auto [min, max] = std::ranges::minmax(filtered);
				stats[0] = min;
				stats[4] = max;
			}

			constexpr double box_width{ 0.8 };
			constexpr double staple_width{ 0.4 };
			custom_plot::patch::line(draw_area, axis_rect,
				QVector<double>{ center, center},
				QVector<double>{ stats[0], stats[1] },
				Qt::black, 1, Qt::DashLine);
			custom_plot::patch::line(draw_area, axis_rect,
				QVector<double>{ center, center},
				QVector<double>{ stats[3], stats[4] },
				Qt::black, 1, Qt::DashLine);
			custom_plot::patch::line(draw_area, axis_rect, 
				QVector<double>{ center - staple_width, center + staple_width }, 
				QVector<double>{ stats[0], stats[0] }, 
				Qt::black);
			custom_plot::patch::line(draw_area, axis_rect,
				QVector<double>{ center - staple_width, center + staple_width },
				QVector<double>{ stats[4], stats[4] },
				Qt::black);
			custom_plot::patch::rectangle(draw_area, axis_rect, center - box_width, stats[1], 2 * box_width, iqr, color, 1, Qt::black);
			custom_plot::patch::line(draw_area, axis_rect,
				QVector<double>{ center - box_width, center + box_width },
				QVector<double>{ stats[2], stats[2] },
				Qt::black, 3);

			if (draw_outlier && outs.size() > 0) {
				custom_plot::patch::scatter(
					draw_area, axis_rect, QVector<double>(outs.size(), center), outs, color, outlier_scatter_size
				);
			}
		};

		void box_batch(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QVector<QColor>& colors,
			const Eigen::ArrayXd& data,
			double start,
			double margin,
			bool draw_outlier,
			int outlier_scatter_size
		) {
			const int n_level = levels.size();

			for (int i = 0; i < n_level; ++i) {

				auto filter = custom::equal(group, levels[i]);
				auto subdata = custom::sliced(data, filter);
				if (subdata.size() == 0) {
					continue;
				}

				custom_plot::patch::box(
					draw_area,
					axis_rect,
					colors[i],
					custom::cast<QVector>(subdata),
					start + i * margin,
					draw_outlier,
					outlier_scatter_size);
			}
		};

		std::pair<double, double> violin(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			QColor color,
			const QVector<double>& value,
			double center,
			int unit,
			double zero_space
		) {
			auto [x, y] = custom_plot::utility::violin_curve(custom::cast<Eigen::ArrayX>(value), center, unit, zero_space);
			custom_plot::patch::shape_borderless(draw_area, axis_rect, x, y, color);
			auto [min, max] = std::ranges::minmax(y);

			return { min, max };
		}

		std::pair<double, double> violin_facet(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QVector<QColor>& colors,
			const Eigen::ArrayXd& data,
			double center,
			int unit
		) {
			auto filter = custom::equal(group, levels[0]);
			auto subdata = custom::sliced(data, filter);
			double min = data.minCoeff(), max = data.maxCoeff(), space = (max - min) / unit;
			double retmin = min, retmax = max;

			if (subdata.size() != 0) {
				auto [x, y] = custom_plot::utility::left_violin_curve(subdata, center);
				custom_plot::patch::shape_borderless(draw_area, axis_rect, x, y, colors[0]);

				custom::extend_minmax(y, retmin, retmax);
			}

			filter = custom::equal(group, levels[1]);
			subdata = custom::sliced(data, filter);

			if (subdata.size() != 0) {
				auto [x, y] = custom_plot::utility::right_violin_curve(subdata, center);
				custom_plot::patch::shape_borderless(draw_area, axis_rect, x, y, colors[1]);

				custom::extend_minmax(y, retmin, retmax);
			}
			return { retmin, retmax };
		};

		std::pair<double, double> violin_batch(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QStringList& group,
			const QStringList& levels,
			const QVector<QColor>& colors,
			const Eigen::ArrayXd& data,
			double start,
			double margin,
			int unit
		) {

			const int n_level = levels.size();

			auto [retmin, retmax] = std::ranges::minmax(data);

			for (int i = 0; i < n_level; ++i) {

				auto filter = custom::equal(group, levels[i]);
				auto subdata = custom::sliced(data, filter);
				if (subdata.size() == 0) {
					continue;
				}

				auto [min, max] = custom_plot::patch::violin(
					draw_area, 
					axis_rect, 
					colors[i], 
					custom::cast<QVector>(subdata), 
					start + margin * i, 
					unit);

				if (min < retmin) {
					retmin = min;
				}

				if (max > retmax) {
					retmax = max;
				}
			}

			return { retmin, retmax };
		};

		QCPAxisRect* new_axis_rect(QCustomPlot* draw_area) {
		
			auto axis_rect = new QCPAxisRect(draw_area, true);

			custom_plot::patch::remove_grid(axis_rect);

			return axis_rect;
		}

		QCPLayoutGrid* set_legend_layout(
			QCustomPlot* draw_area,
			QCPLayoutGrid* legend_layout) {

			QCPLayoutGrid* layout = new QCPLayoutGrid;
			QCPAxisRect* legend_top = new QCPAxisRect(draw_area, false);
			QCPAxisRect* legend_bottom = new QCPAxisRect(draw_area, false);

			legend_top->setMinimumSize(QSize(0, 0));
			legend_top->setMargins(QMargins(0, 0, 0, 0));

			legend_bottom->setMinimumSize(QSize(0, 0));
			legend_bottom->setMargins(QMargins(0, 0, 0, 0));

			legend_layout->addElement(0, 0, legend_top);
			legend_layout->addElement(1, 0, layout);
			legend_layout->addElement(2, 0, legend_bottom);

			layout->setRowSpacing(15);

			legend_layout->setMargins(QMargins(0, 0, 20, 0));

			return layout;
		}

		void add_gradient_legend(
			QCustomPlot* draw_area,
			QCPLayoutGrid* legend_layout,
			double minimum_value,
			double maximum_value,
			const QString& title,
			const QString& lower_label,
			const QString& upper_label,
			const QFont& legend_title_font,
			const QFont& legend_label_font,
			QColor low_color,
			QColor middle_color,
			QColor high_color
		) {
			QCPLayoutGrid* sub_legend_layout = new QCPLayoutGrid;
			QCPLayoutGrid* color_scale_layout = new QCPLayoutGrid;
			QCPAxisRect* color_scale_left_rect = new QCPAxisRect(draw_area, false);
			QCPAxisRect* color_scale_right_rect = new QCPAxisRect(draw_area, false);
			sub_legend_layout->setRowSpacing(20);
			auto [row, col] = custom_plot::patch::find_next_empty_position(legend_layout);
			legend_layout->addElement(row, col, sub_legend_layout);

			QCPColorGradient gradient;
			gradient.setColorStopAt(0.0, low_color);
			gradient.setColorStopAt(0.5, middle_color);
			gradient.setColorStopAt(1.0, high_color);
			QCPColorScale* color_scale = new QCPColorScale(draw_area);
			color_scale->setType(QCPAxis::atRight);
			color_scale->setDataRange(QCPRange(minimum_value, maximum_value));
			color_scale->setGradient(gradient);
			QSharedPointer<QCPAxisTickerText> ticker(new QCPAxisTickerText);
			ticker->setTicks(QVector<double>{minimum_value, maximum_value}, { lower_label, upper_label });

			color_scale->axis()->setTicks(true);
			color_scale->axis()->setTickLength(0);
			color_scale->axis()->setTickPen(Qt::NoPen);
			color_scale->axis()->setTickLabels(true);
			color_scale->axis()->setTicker(ticker);
			color_scale->axis()->setTickLabelFont(legend_label_font);

			color_scale->mAxisRect->axis(QCPAxis::atBottom)->setBasePen(Qt::NoPen);
			color_scale->mAxisRect->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
			color_scale->mAxisRect->axis(QCPAxis::atTop)->setBasePen(Qt::NoPen);
			color_scale->mAxisRect->axis(QCPAxis::atRight)->setBasePen(Qt::NoPen);

			color_scale_layout->addElement(0, 0, color_scale_left_rect);
			color_scale_layout->addElement(0, 1, color_scale);
			color_scale_layout->addElement(0, 2, color_scale_right_rect);
			color_scale->setBarWidth(28);
			color_scale_layout->setMinimumSize(200, 200);
			color_scale_layout->setMaximumSize(200, 200);
			sub_legend_layout->addElement(0, 0, color_scale_layout);
			custom_plot::patch::add_title(draw_area, sub_legend_layout, title, legend_title_font);
		}

		void set_fixed_width(QCPAxisRect* axis_rect, int width) {

			axis_rect->setMinimumSize(width, 0);
			axis_rect->setMaximumSize(width, 8000);
		};

		void set_fixed_height(QCPAxisRect* axis_rect, int height) {

			axis_rect->setMinimumSize(0, height);
			axis_rect->setMaximumSize(8000, height);
		};

		std::pair<int, int> find_next_empty_position(QCPLayoutGrid* layout) {
			int nrow = layout->rowCount();
			int ncol = layout->columnCount();

			for (int row = 0; row < nrow; ++row) {
				for (int col = 0; col < ncol; ++col) {

					QCPLayoutElement* element = layout->element(row, col);

					if (element == nullptr) {
						return { row, col };
					}
				}
			}

			if (nrow > ncol + 1) {
				return { 0, ncol };
			}
			else {
				return { nrow, 0 };
			}
		};

		void set_axis_label(
			QCPAxis::AxisType type,
			QCPAxisRect* axis_rect,
			const Eigen::ArrayXd& location,
			const QStringList& labels,
			int tick_length,
			QFont label_font,
			int label_angle
		) {					
			
			QSharedPointer<QCPAxisTickerText> ticker(new QCPAxisTickerText);
			ticker->setTicks(QVector<double>(location.begin(), location.end()), labels);
			axis_rect->axis(type)->setTicker(ticker);

			axis_rect->axis(type)->setTicks(true);
			axis_rect->axis(type)->setTickLength(0, tick_length);
			axis_rect->axis(type)->setTickLabels(true);
			axis_rect->axis(type)->setTickLabelFont(label_font);
			axis_rect->axis(type)->setTickLabelRotation(label_angle);

			if (tick_length == 0) {
				
				axis_rect->axis(type)->setTickPen(Qt::NoPen);
			}
			else {
								
				QPen pen(Qt::black);
				pen.setWidth(3);
				axis_rect->axis(type)->setTickPen(pen);
			}

		};

		void add_round_legend(
			QCustomPlot* draw_area,
			QCPLayoutGrid* legend_layout,
			const QStringList& levels,
			const QList<QColor>& colors,
			const QString& legend_title,
			int legend_column_width,
			int legend_row_width,
			QFont legend_title_font,
			QFont legend_label_font
		)
		{
			QCPLayoutGrid* sub_legend_layout = new QCPLayoutGrid, * item_layout_ = new QCPLayoutGrid;
			sub_legend_layout->setMargins(QMargins(0, 0, 0, 0));
			item_layout_->setMargins(QMargins(0, 0, 0, 0));

			sub_legend_layout->setRowSpacing(0);

			auto [row, col] = custom_plot::patch::find_next_empty_position(legend_layout);
			legend_layout->addElement(row, col, sub_legend_layout);

			sub_legend_layout->addElement(0, 0, item_layout_);

			double size = levels.size();
			int ncol = ceil(size / 20), nrow = ceil(size / ncol), index = 0;
			QVector<int> rows(ncol, nrow);
			rows[ncol - 1] = size - (ncol - 1) * nrow;

			QTextDocument td;
			td.setDefaultFont(legend_label_font);

			for (int i = 0; i < ncol; ++i) {
				int column_width = 0;
				int row = rows[i];
				QCPAxisRect* legend_rect = new QCPAxisRect(draw_area, true);
				item_layout_->addElement(0, i, legend_rect);
				int j;
				for (j = 0; j < row; ++j, ++index) {
					custom_plot::patch::set_single_round_legend(draw_area, legend_rect, levels[index], colors[index], legend_label_font, 10, row - j - 1);
					td.setHtml(levels[index]);
					auto [width, _] = td.size().toSize();
					column_width = column_width > width ? column_width : width;
				}
				if (legend_column_width == 0) {
					legend_column_width = column_width + 20;
				}
				custom_plot::patch::remove_left_bottom_axis(legend_rect);
				custom_plot::patch::set_fixed_size(legend_rect, legend_column_width, legend_row_width * j);
				custom_plot::patch::set_range(legend_rect, QCPRange(-10, legend_column_width - 10), QCPRange(-0.5, j - 0.5));
			}
			custom_plot::patch::add_title(draw_area, sub_legend_layout, legend_title, legend_title_font);
		};

		void single_violin_plot(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& values,
			int unit,
			QColor color,
			const QString& bottom_title,
			const QString& left_title)
		{
			
			QPen pen(Qt::black);
			pen.setWidth(3);
			axis_rect->axis(QCPAxis::atLeft)->setBasePen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setTickPen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setSubTicks(false);
			axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 6);

			axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
			axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
			axis_rect->axis(QCPAxis::atBottom)->setSubTicks(false);
			axis_rect->axis(QCPAxis::atBottom)->setTickLength(0, 6);

			axis_rect->axis(QCPAxis::atBottom)->setUpperEnding(QCPLineEnding::esSpikeArrow);
			axis_rect->axis(QCPAxis::atLeft)->setUpperEnding(QCPLineEnding::esSpikeArrow);

			axis_rect->axis(QCPAxis::atBottom)->setLabel(bottom_title);
			axis_rect->axis(QCPAxis::atBottom)->setLabelFont(QFont("Arial", 18, QFont::Bold));

			axis_rect->axis(QCPAxis::atLeft)->setLabel(left_title);
			axis_rect->axis(QCPAxis::atLeft)->setLabelFont(QFont("Arial", 18, QFont::Bold));

			auto [min, max] = custom_plot::patch::violin(
				draw_area,
				axis_rect,
				color,
				values,
				1.0,
				unit
			);

			custom_plot::patch::set_range(axis_rect, { 0.0, 2.0 }, custom_plot::utility::get_range(min, max));
		};

		void single_box_plot(
			QCustomPlot* draw_area,
			QCPAxisRect* axis_rect,
			const QVector<double>& values,
			QColor color,
			const QString& bottom_title,
			const QString& left_title,
			bool draw_outlier,
			int outlier_scatter_size)
		{
			QPen pen(Qt::black);
			pen.setWidth(3);
			axis_rect->axis(QCPAxis::atLeft)->setBasePen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setTickPen(pen);
			axis_rect->axis(QCPAxis::atLeft)->setSubTicks(false);
			axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 6);

			axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
			axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
			axis_rect->axis(QCPAxis::atBottom)->setSubTicks(false);
			axis_rect->axis(QCPAxis::atBottom)->setTickLength(0, 6);

			axis_rect->axis(QCPAxis::atBottom)->setUpperEnding(QCPLineEnding::esSpikeArrow);
			axis_rect->axis(QCPAxis::atLeft)->setUpperEnding(QCPLineEnding::esSpikeArrow);

			axis_rect->axis(QCPAxis::atBottom)->setLabel(bottom_title);
			axis_rect->axis(QCPAxis::atBottom)->setLabelFont(QFont("Arial", 18, QFont::Bold));

			axis_rect->axis(QCPAxis::atLeft)->setLabel(left_title);
			axis_rect->axis(QCPAxis::atLeft)->setLabelFont(QFont("Arial", 18, QFont::Bold));

			custom_plot::patch::box(
				draw_area,
				axis_rect,
				color,
				values,
				1.0,
				draw_outlier,
				outlier_scatter_size
			);

			custom_plot::patch::set_range(axis_rect, { 0.0, 2.0 }, custom_plot::utility::get_range(values));
		};
	};
};