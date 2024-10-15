#include "CustomPlot.h"

#include "SoapTextElement.h"
#include "Custom.h"

auto operator<=>(const QColor& color1, const QColor& color2) {
	if (auto cmp = color1.red() <=> color2.red(); cmp != 0) {
		return cmp;
	}
	if (auto cmp = color1.green() <=> color2.green(); cmp != 0) {
		return cmp;
	}
	if (auto cmp = color1.blue() <=> color2.blue(); cmp != 0) {
		return cmp;
	}

	return color1.alpha() <=> color2.alpha();
}


namespace custom_plot {

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> feature_plot(
		const QUERY_DATA& data,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		bool scale,
		const GraphSettings& gs) {

		auto [draw_area, layout] = custom_plot::prepare_lg(gs);
		auto [axis_rect, legend_layout] = custom_plot::__feature_plot(
			draw_area,
			layout,
			data,
			x,
			y,
			bottom_title,
			left_title,
			scale,
			gs
		);

		return { draw_area, axis_rect, legend_layout };
	}

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> feature_plot(
		const QUERY_DATA& data,
		Embedding* embedding,
		bool scale,
		const GraphSettings& gs) {

		Eigen::ArrayXd x = embedding->data_.mat_.col(0);
		Eigen::ArrayXd y = embedding->data_.mat_.col(1);
		QString bottom_title = embedding->data_.colnames_[0];
		QString left_title = embedding->data_.colnames_[1];

		auto [draw_area, layout] = custom_plot::prepare_lg(gs);
		auto [axis_rect, legend_layout] = custom_plot::__feature_plot(
			draw_area,
			layout,
			data,
			x,
			y,
			bottom_title,
			left_title,
			scale,
			gs
		);

		return { draw_area, axis_rect, legend_layout };
	}

	std::pair<QCPAxisRect*, QCPLayoutGrid*> __feature_plot(
		QCustomPlot* draw_area,
		QCPLayoutGrid* sub_layout,
		const QUERY_DATA& data,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		bool scale,
		const GraphSettings& gs) {

		QString legend_title{ "Expression" };

		if (data.info.contains("Source") && (data.info["Source"] == "ATAC" || data.info["Source"] == "Gene Activity")) {
			legend_title = "Activity";
		}
		if (data.info.contains("Source") && (data.info["Source"] == "Metadata" || data.info["Source"] == "Others")) {
			legend_title = data.name;
		}
		if (data.info.contains("Legend Title")) {
			legend_title = data.info["Legend Title"];
		}

		bool lighten_color = data.message.contains("Lighten Color");

		QCPLayoutGrid* right_layout = new QCPLayoutGrid;

		QCPAxisRect* axis_rect = custom_plot::patch::new_axis_rect(draw_area);
		sub_layout->addElement(0, 0, axis_rect);
		sub_layout->addElement(0, 1, right_layout);
		auto legend_layout = custom_plot::patch::set_legend_layout(draw_area, right_layout);

		custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(x, y));
		custom_plot::set_scatter_plot_axis_style(draw_area, axis_rect, bottom_title, left_title, x, y, gs);
		if (data.is_continuous()) {

			QColor low_color = gs.get_gradient_low_color();
			QColor middle_color = gs.get_gradient_middle_color();
			QColor high_color = gs.get_gradient_high_color();

			if (lighten_color) {
				low_color.setAlpha(64);
				middle_color.setAlpha(64);
				high_color.setAlpha(64);
			}

			Eigen::ArrayXd d = custom::cast<Eigen::ArrayX>(data.get_continuous());
			auto [min, max] = std::ranges::minmax(d);

			if (scale) {

				double sd = custom::sd(d);

				if (sd != 0.0) {

					double mean = d.mean();
					d -= mean;
					d /= sd;

					custom_plot::patch::scatter_gradient(
						draw_area,
						axis_rect,
						x,
						y,
						d,
						-1.0,
						1.0,
						low_color,
						middle_color,
						high_color,
						gs.get_scatter_point_size()
					);

					custom_plot::patch::add_gradient_legend(
						draw_area,
						legend_layout,
						-1.0,
						1.0,
						legend_title,
						"Low",
						"High",
						gs.get_legend_title_font(),
						gs.get_legend_label_font(),
						low_color,
						middle_color,
						high_color
					);

					custom_plot::patch::add_title(draw_area, sub_layout, data.name, gs.get_title_font());
				}
				else {

					custom_plot::patch::scatter_gradient(
						draw_area,
						axis_rect,
						x,
						y,
						d,
						min,
						max,
						low_color,
						middle_color,
						high_color,
						gs.get_scatter_point_size()
					);

					custom_plot::add_gradient_legend(draw_area, legend_layout, min, max, legend_title, gs, low_color, middle_color, high_color);

					custom_plot::patch::add_title(draw_area, sub_layout, data.name, gs.get_title_font());
				}

			}
			else {

				custom_plot::patch::scatter_gradient(
					draw_area,
					axis_rect,
					x,
					y,
					d,
					min,
					max,
					low_color,
					middle_color,
					high_color,
					gs.get_scatter_point_size()
				);
				custom_plot::add_gradient_legend(draw_area, legend_layout, min, max, legend_title, gs, low_color, middle_color, high_color);
				custom_plot::patch::add_title(draw_area, sub_layout, data.name, gs.get_title_font());
			}
		}
		else if (data.is_factor()) {
			auto factors = data.get_factor();
			auto levels = data.get_levels();

			auto colors = gs.palette(levels);

			if (lighten_color) {
				for (auto&& c : colors) {
					c.setAlpha(64);
				}
			}

			custom_plot::patch::scatter_category(draw_area, axis_rect, x, y, factors, levels, colors, gs.get_scatter_point_size());

			custom_plot::add_round_legend(draw_area, legend_layout, levels, colors, data.name, gs);

			if (gs.is_scatter_labeled()) {

				for (auto&& level : levels) {

					QCPItemText* scatter_label = new QCPItemText(draw_area);
					scatter_label->position->setAxisRect(axis_rect);
					scatter_label->position->setAxes(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
					scatter_label->position->setType(QCPItemPosition::ptPlotCoords);
					scatter_label->setPositionAlignment(Qt::AlignCenter);
					scatter_label->position->setCoords(custom::sliced(x, custom::equal(factors, level)).mean(), custom::sliced(y, custom::equal(factors, level)).mean());
					scatter_label->setText(level);
					scatter_label->setFont(gs.get_scatter_label_font());
				}
			}

			custom_plot::patch::add_title(draw_area, sub_layout, gs.get_title(data.name), gs.get_title_font());
		}

		return { axis_rect, legend_layout };
	};
	
	QCustomPlot* monocle3_feature_plot(
		const QList<QUERY_DATA>& data,
		const Eigen::ArrayXd& x,
		const QUERY_DATA& f,
		const QString& bottom_title,
		int nrow,
		const GraphSettings& gs) {

		int n_cell = x.size();
		int n_block = 100;

		if (n_block > n_cell) {
			n_block = n_cell;
		}

		int times{ 0 };

		if (n_cell % n_block == 0) {
			times = n_cell / n_block;
		}
		else {
			times = n_cell / n_block + 1;
		}

		auto [min, max] = std::ranges::minmax(x);
		double span = (max - min) / (n_block - 1);
		double limit = min - 0.5 * span;
		Eigen::ArrayXi block_loc = Eigen::ArrayXi::Zero(n_cell);
		Eigen::ArrayXi block_count = Eigen::ArrayXi::Zero(n_block);
		for (int i = 0; i < n_cell; ++i) {
			int loc = (x[i] - limit) / span;
			++block_count[loc];
			block_loc[i] = loc;
		}

		Eigen::ArrayXd loc_multiplier(n_cell);
		for (int i = 0; i < n_cell; ++i) {
			loc_multiplier[i] = 1.0 / block_count[block_loc[i]];
		}

		auto [draw_area, layout] = custom_plot::prepare_lg(gs);
		int ind{ 0 };
		for (auto&& d : data) {
			if (!d.is_continuous()) {
				continue;
			}

			Eigen::ArrayXd y = custom::cast<Eigen::ArrayX>(d.get_continuous());

			QString left_title{ "Expression of " + d.name };

			if (d.info.contains("Source")) {
				if (d.info["Source"] == "Gene Activity" || d.info["Source"] == "ATAC") {
					left_title = "Activity of " + d.name;
				}
				if (d.info["Source"] == "Metadata" || d.info["Source"] == "Others") {
					left_title = d.name;
				}
			}

			int col = ind / nrow;
			int row = ind - col * nrow;
			QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
			layout->addElement(row, col, sub_layout);
			auto [axis_rect, _] = custom_plot::__feature_plot(draw_area, sub_layout, f, x, y, bottom_title, left_title, false, gs);

			Eigen::ArrayXd x1 = Eigen::ArrayXd::LinSpaced(n_block, min, max);
			Eigen::ArrayXd y1 = Eigen::ArrayXd::Zero(n_block);

			for (int i = 0; i < n_cell; ++i) {
				y1[block_loc[i]] += y[i] * loc_multiplier[i];
			}

			custom_plot::patch::line(
				draw_area,
				axis_rect,
				x1,
				y1,
				Qt::black,
				2
			);

			++ind;
		}

		return draw_area;
	};

	QCustomPlot* feature_plot(
		const QList<QUERY_DATA>& data,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title,
		const QString& left_title,
		bool scale,
		int nrow,
		const GraphSettings& gs) {

		auto [draw_area, layout] = custom_plot::prepare_lg(gs);
		int ind{ 0 };
		for (auto&& d : data) {
			if (!d.is_valid()) {
				continue;
			}
			int col = ind / nrow;
			int row = ind - col * nrow;
			QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
			layout->addElement(row, col, sub_layout);
			custom_plot::__feature_plot(draw_area, sub_layout, d, x, y, bottom_title, left_title, scale, gs);
			++ind;
		}

		return draw_area;
	};

	QCustomPlot* feature_plot(
		const QList<QUERY_DATA>& data, 
		Embedding* embedding,
		bool scale, 
		int nrow, 
		const GraphSettings& gs)
	{

		Eigen::ArrayXd x = embedding->data_.mat_.col(0);
		Eigen::ArrayXd y = embedding->data_.mat_.col(1);
		QString bottom_title = embedding->data_.colnames_[0];
		QString left_title = embedding->data_.colnames_[1];

		auto [draw_area, layout] = custom_plot::prepare_lg(gs);
		int ind{ 0 };
		for (auto&& d : data) {
			if (!d.is_valid()) {
				continue;
			}
			int col = ind / nrow;
			int row = ind - col * nrow;
			QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
			layout->addElement(row, col, sub_layout);
			custom_plot::__feature_plot(draw_area, sub_layout, d, x, y, bottom_title, left_title, scale, gs);
			++ind;
		}

		return draw_area;
	};

	QCustomPlot* coverage_plot(COVERAGE_PLOT_ELEMENTS res, const GraphSettings& gs) {

		bool draw_ccan = !res.ccan_locations.isEmpty();

		const int width = res.normalized_matrix.cols();

		auto& [sequence_name, region_start, region_end] = res.region;

		QVector<double> bottom_axis_value = custom::linspaced(width, region_start, region_end);

		auto [draw_area, main_layout] = custom_plot::prepare_lg(gs);
		draw_area->setNoAntialiasingOnDrag(true);

		QPen annotation_axis_pen(Qt::black), coverage_bottom_axis_pen(Qt::black);
		annotation_axis_pen.setWidth(3);
		coverage_bottom_axis_pen.setWidth(1);

		QCPMarginGroup* first_column_margin = new QCPMarginGroup(draw_area);
		QCPMarginGroup* second_column_margin = new QCPMarginGroup(draw_area);

		int anno_width = std::ceil(custom_plot::utility::get_max_text_height({ "Coverage", "Gene", "Peak" }, gs.get_left_label_font()) * 1.2);

		QCPAxisRect* coverage_annotation_rect = new QCPAxisRect(draw_area);
		main_layout->addElement(0, 0, coverage_annotation_rect);

		custom_plot::patch::remove_left_bottom_axis(coverage_annotation_rect);
		custom_plot::patch::rectangle_borderless(
			draw_area, coverage_annotation_rect, 0, 0, 1, 1, Qt::black
		);
		custom_plot::patch::add_label(draw_area, coverage_annotation_rect, "Coverage", -1.0, 0.5, gs.get_left_label_font(), Qt::AlignBottom | Qt::AlignHCenter, -90);
		custom_plot::patch::set_range(coverage_annotation_rect, { -11.0, 1.0 }, { 0.0, 1.0 });
		custom_plot::patch::set_fixed_width(coverage_annotation_rect, anno_width);
		coverage_annotation_rect->setMarginGroup(QCP::msLeft | QCP::msRight, first_column_margin);

		QCPLayoutGrid* coverage_layout = new QCPLayoutGrid;
		coverage_layout->setRowSpacing(0);
		main_layout->addElement(0, 1, coverage_layout);
		main_layout->setColumnSpacing(0);

		const int group_size = res.group_factors.size();

		auto colors = gs.palette(res.group_factors);

		Eigen::ArrayXd group_factor_name_height(1);
		group_factor_name_height[0] = 0.5;

		for (int i = 0; i < group_size; ++i) {
			QCPAxisRect* group_rect = new QCPAxisRect(draw_area);
			group_rect->setMinimumSize(20, 30);
			group_rect->setMargins({ 0, 0, 0, 0 });
			group_rect->setMinimumMargins({ 0, 0, 0, 0 });
			coverage_layout->addElement(i, 0, group_rect);

			custom_plot::set_left_axis_label(
				group_rect,
				group_factor_name_height,
				{ res.group_factors[i] },
				0,
				gs
			);
			group_rect->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
			custom_plot::patch::remove_left_grid(group_rect);
			custom_plot::patch::clear_bottom_axis(group_rect);
			group_rect->axis(QCPAxis::atBottom)->setBasePen(coverage_bottom_axis_pen);
			group_rect->setMarginGroup(QCP::msLeft | QCP::msRight, second_column_margin);

			QVector<double> normalized_value = custom::cast<QVector>(res.normalized_matrix.row(i));
			const Eigen::Index length = normalized_value.size();

			Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(length, true);

			for (Eigen::Index i = 1; i < length - 1; ++i) {
				if (normalized_value[i] < 0.01) {
					if (normalized_value[i - 1] < 0.01 && normalized_value[i + 1] < 0.01) {
						filter[i] = false;
					}
					else {
						normalized_value[i] = 0;
					}
				}
			}

			QCPGraph* graph = draw_area->addGraph(group_rect->axis(QCPAxis::atBottom), group_rect->axis(QCPAxis::atLeft));

			graph->setPen(Qt::NoPen);
			graph->setBrush(QBrush(colors[i]));
			graph->setData(custom::sliced(bottom_axis_value, filter), custom::sliced(normalized_value, filter), true);

			if (draw_ccan) {
				for (auto&& p : res.ccan_locations) {
					custom_plot::patch::rectangle_borderless(draw_area, group_rect, p.first, 0, p.second - p.first, 1, QColor(221, 221, 221, 32));
				}
			}

			custom_plot::patch::set_range(group_rect, QCPRange(region_start, region_end + 1), QCPRange(0, 1));

		}

		//-----------------------------------//

		if (!res.draw_bac) {
			QStringList gene_names = res.gene_structure.keys();

			for (auto&& gene_name : gene_names) {
				if (gene_name.startsWith("RP") && gene_name.contains('.')) {
					res.gene_structure.remove(gene_name);
				}
			}
		}

		if (!res.draw_ac) {
			QStringList gene_names = res.gene_structure.keys();

			for (auto&& gene_name : gene_names) {
				if (gene_name.startsWith("AC") && gene_name.contains('.')) {
					res.gene_structure.remove(gene_name);
				}
			}
		}

		bool draw_gene = res.gene_structure.size() > 0;

		if (draw_gene) {

			//--------------------------------//
			QCPAxisRect* gene_annotation_rect = new QCPAxisRect(draw_area);
			main_layout->addElement(1, 0, gene_annotation_rect);

			custom_plot::patch::remove_left_bottom_axis(gene_annotation_rect);
			custom_plot::patch::rectangle_borderless(
				draw_area, gene_annotation_rect, 0, 0, 1, 1, Qt::black
			);
			custom_plot::patch::add_label(draw_area, gene_annotation_rect, "Gene", -1, 0.5, gs.get_left_label_font(), Qt::AlignBottom | Qt::AlignHCenter, -90);
			custom_plot::patch::set_range(gene_annotation_rect, { -11.0, 1.0 }, { 0.0, 1.0 });
			custom_plot::patch::set_fixed_width(gene_annotation_rect, anno_width);

			gene_annotation_rect->setMarginGroup(QCP::msLeft | QCP::msRight, first_column_margin);
			//--------------------------------//

			QCPLayoutGrid* gene_layout = new QCPLayoutGrid;
			gene_layout->setRowSpacing(0);
			main_layout->addElement(1, 1, gene_layout);

			const int gene_size = res.gene_structure.size();
			int gene_count = 0;

			QList<std::tuple<QString, int, int>> locs;

			for (const auto& [gene_name, info] : res.gene_structure.asKeyValueRange()) {
				int s = std::ranges::min(custom::sapply(info.second, [](auto&& p) {return std::get<0>(p); }));
				int e = std::ranges::max(custom::sapply(info.second, [](auto&& p) {return std::get<1>(p); }));

				int span = region_end - region_start;
				s -= span / 20;
				e += span / 20;

				if (e - s < (span / 6)) {
					s -= span / 30;
					e += span / 30;
				}

				if (s < region_start) {
					s = region_start;
				}

				if (e > region_end) {
					e = region_end;
				}

				locs << std::make_tuple(gene_name, s, e);
			}

			while (locs.size() > 0) {
				QCPAxisRect* gene_rect = new QCPAxisRect(draw_area);
				gene_rect->setMarginGroup(QCP::msLeft | QCP::msRight, second_column_margin);
				gene_layout->addElement(gene_count++, 0, gene_rect);
				custom_plot::patch::remove_left_bottom_axis(gene_rect);
				custom_plot::patch::set_fixed_height(gene_rect, 25);

				QVector<int> current_loc;
				current_loc << region_start << region_end;

				bool update{ true };

				while (update) {
					update = false;

					int n_loc = locs.size();

					if (n_loc == 0) {
						break;
					}

					for (int i = 0; i < n_loc; ++i) {
						auto&& [name, s, e] = locs[i];

						bool allow{ false };

						int current_loc_size = current_loc.size() / 2;
						for (int j = 0; j < current_loc_size; ++j) {
							if (s >= current_loc[j * 2] && e <= current_loc[j * 2 + 1]) {
								allow = true;
								break;
							}
						}

						if (!allow) {
							continue;
						}

						const auto& location = res.gene_structure[name];
						char strand = location.first;
						const auto& location_list = location.second;

						double gene_start = -1.0, gene_end = -1.0;
						for (const auto& segment : location_list) {

							auto&& [start, end, type] = segment;

							if (type != "gap") {
								continue;
							}

							double segment_start = start < region_start ? region_start : start;
							double segment_end = end > region_end ? region_end : end;

							if (gene_start < 0) {
								gene_start = segment_start;
								gene_end = segment_end;
							}
							else {
								if (start < gene_start) {
									gene_start = segment_start;
								}
								if (end > gene_end) {
									gene_end = segment_end;
								}
							}

							custom_plot::patch::shape_borderless(
								draw_area,
								gene_rect,
								QVector<double>{ segment_start, segment_start, segment_end, segment_end },
								QVector<double>{9.0, 11.0, 11.0, 9.0},
								Qt::black);
						}

						for (const auto& segment : location_list) {

							auto&& [start, end, type] = segment;

							if (type != "exon") {
								continue;
							}

							double segment_start = start < region_start ? region_start : start;
							double segment_end = end > region_end ? region_end : end;

							if (gene_start < 0) {
								gene_start = segment_start;
								gene_end = segment_end;
							}
							else {
								if (start < gene_start) {
									gene_start = segment_start;
								}
								if (end > gene_end) {
									gene_end = segment_end;
								}
							}

							custom_plot::patch::shape(draw_area, gene_rect,
								QVector<double>{ segment_start, segment_start, segment_end, segment_end },
								QVector<double>{2.0, 18.0, 18.0, 2.0},
								Qt::black, Qt::black, 2);
						}

						for (const auto& segment : location_list) {

							auto&& [start, end, type] = segment;

							if (type != "utr") {
								continue;
							}

							double segment_start = start < region_start ? region_start : start;
							double segment_end = end > region_end ? region_end : end;

							if (gene_start < 0) {
								gene_start = segment_start;
								gene_end = segment_end;
							}
							else {
								if (start < gene_start) {
									gene_start = segment_start;
								}
								if (end > gene_end) {
									gene_end = segment_end;
								}
							}

							custom_plot::patch::shape(draw_area, gene_rect,
								QVector<double>{ segment_start, segment_start, segment_end, segment_end },
								QVector<double>{2.0, 18.0, 18.0, 2.0},
								Qt::black, Qt::white, 2);
						}

						for (const auto& segment : location_list) {

							auto&& [start, end, type] = segment;

							if (type != "cds") {
								continue;
							}

							double segment_start = start < region_start ? region_start : start;
							double segment_end = end > region_end ? region_end : end;

							if (gene_start < 0) {
								gene_start = segment_start;
								gene_end = segment_end;
							}
							else {
								if (start < gene_start) {
									gene_start = segment_start;
								}
								if (end > gene_end) {
									gene_end = segment_end;
								}
							}

							custom_plot::patch::shape(draw_area, gene_rect,
								QVector<double>{ segment_start, segment_start, segment_end, segment_end },
								QVector<double>{2.0, 18.0, 18.0, 2.0},
								custom_plot::color::firebrick3, custom_plot::color::firebrick3, 2);
						}

						double span = region_end - region_start;
						if (gene_end < region_end) {
							double ar_start = span / 60 + gene_end;

							double ar_span = span / 30;
							double ar_end = ar_start + ar_span;
							double ar_mid = ar_start + ar_span / 2;
							double ar_mid1 = ar_start + ar_span / 4;
							double ar_mid2 = ar_end - ar_span / 4;

							if (strand == '+') {
								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_start, ar_mid1, ar_mid, ar_mid1, ar_start, ar_mid1},
									QVector<double>{2.0, 2.0, 10.0, 18.0, 18.0, 10.0},
									custom_plot::color::navy);

								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_mid, ar_mid2, ar_end, ar_mid2, ar_mid, ar_mid2},
									QVector<double>{2.0, 2.0, 10.0, 18.0, 18.0, 10.0},
									custom_plot::color::navy);
							}
							else if (strand == '-') {
								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_mid1, ar_start, ar_mid1, ar_mid, ar_mid1, ar_mid},
									QVector<double>{2.0, 10.0, 18.0, 18.0, 10.0, 2.0},
									custom_plot::color::navy);

								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_mid2, ar_mid, ar_mid2, ar_end, ar_mid2, ar_end},
									QVector<double>{2.0, 10.0, 18.0, 18.0, 10.0, 2.0},
									custom_plot::color::navy);
							}
						}

						if (gene_start > region_start) {
							double ar_end = gene_start - span / 60;

							double ar_span = span / 30;
							double ar_start = ar_end - ar_span;
							double ar_mid = ar_start + ar_span / 2;
							double ar_mid1 = ar_start + ar_span / 4;
							double ar_mid2 = ar_end - ar_span / 4;

							if (strand == '+') {
								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_start, ar_mid1, ar_mid, ar_mid1, ar_start, ar_mid1},
									QVector<double>{2.0, 2.0, 10.0, 18.0, 18.0, 10.0},
									custom_plot::color::navy);

								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_mid, ar_mid2, ar_end, ar_mid2, ar_mid, ar_mid2},
									QVector<double>{2.0, 2.0, 10.0, 18.0, 18.0, 10.0},
									custom_plot::color::navy);
							}
							else if (strand == '-') {
								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_mid1, ar_start, ar_mid1, ar_mid, ar_mid1, ar_mid},
									QVector<double>{2.0, 10.0, 18.0, 18.0, 10.0, 2.0},
									custom_plot::color::navy);

								custom_plot::patch::shape_borderless(draw_area, gene_rect,
									QVector<double>{ar_mid2, ar_mid, ar_mid2, ar_end, ar_mid2, ar_end},
									QVector<double>{2.0, 10.0, 18.0, 18.0, 10.0, 2.0},
									custom_plot::color::navy);
							}
						}
						custom_plot::patch::set_range(gene_rect, QCPRange(region_start, region_end + 1), QCPRange(0, 20));
						custom_plot::patch::add_label(draw_area, gene_rect, name, (gene_start + gene_end) / 2, 0, gs.get_bottom_label_font(), Qt::AlignTop | Qt::AlignHCenter);

						update = true;

						current_loc << s << e;

						std::ranges::sort(current_loc);

						locs.removeAt(i);

						break;
					}
				}
			}

		}

		int layer = (int)draw_gene + 1;

		//------------------------------------//
		bool draw_link = res.peak_links.size() > 0;

		if (draw_link) {
			QCPAxisRect* link_annotation_rect = new QCPAxisRect(draw_area);
			QCPAxisRect* link_rect = new QCPAxisRect(draw_area);
			link_rect->setMarginGroup(QCP::msLeft | QCP::msRight, second_column_margin);

			main_layout->addElement(layer, 0, link_annotation_rect);
			main_layout->addElement(layer, 1, link_rect);

			custom_plot::patch::remove_left_bottom_axis(link_annotation_rect);
			custom_plot::patch::rectangle_borderless(
				draw_area, link_annotation_rect, 0, 0, 1, 1, Qt::black
			);
			custom_plot::patch::add_label(draw_area, link_annotation_rect, "Link", -1, 0.5, gs.get_left_label_font(), Qt::AlignBottom | Qt::AlignHCenter, -90);
			custom_plot::patch::set_range(link_annotation_rect, { -11.0, 1.0 }, { 0.0, 1.0 });

			constexpr int link_height = 60;
			custom_plot::patch::set_fixed_size(link_annotation_rect, anno_width, link_height);
			custom_plot::patch::set_fixed_height(link_rect, link_height);

			link_annotation_rect->setMarginGroup(QCP::msLeft | QCP::msRight, first_column_margin);

			auto max_link_val = std::ranges::max(custom::sapply(res.peak_links, [](auto&& l) {return std::get<2>(l); }));

			for (auto&& [p1, p2, link_val] : res.peak_links) {

				auto [start1, end1] = res.peak_locations[p1];
				double s = (start1 + end1) / 2;
				auto [start2, end2] = res.peak_locations[p2];
				double e = (start2 + end2) / 2;

				double mid = (s + e) / 2;
				double a = -link_val / ((s - mid) * (s - mid));

				auto link_x = custom::linspaced(100, s, e);
				auto link_y = custom::sapply(link_x, [a, mid, link_val](auto x) {return a * (x - mid) * (x - mid) + link_val; });

				custom_plot::patch::line(draw_area, link_rect, link_x, link_y, custom_plot::color::gray, 2);
			}

			custom_plot::patch::set_range(link_rect, QCPRange(region_start, region_end + 1), { 0.0, 1.0 });

			custom_plot::patch::remove_bottom_axis(link_rect);
			link_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
			custom_plot::set_left_axis_label(link_rect, custom::cast<Eigen::ArrayX>(QVector<double>{0.0, 1.0}), { "0.0", "1.0" }, 2, gs);

			QPen p(Qt::black);
			p.setWidth(2);
			link_rect->axis(QCPAxis::atLeft)->setBasePen(p);
		}

		layer += (int)draw_link;
		//------------------------------//
		QCPAxisRect* peak_annotation_rect = new QCPAxisRect(draw_area);
		QCPAxisRect* peak_rect = new QCPAxisRect(draw_area);

		main_layout->addElement(layer, 0, peak_annotation_rect);
		main_layout->addElement(layer, 1, peak_rect);

		custom_plot::patch::remove_left_bottom_axis(peak_annotation_rect);
		custom_plot::patch::rectangle_borderless(
			draw_area, peak_annotation_rect, 0, 0, 1, 1, Qt::black
		);
		custom_plot::patch::add_label(draw_area, peak_annotation_rect, "Peak", -1, 0.5, gs.get_left_label_font(), Qt::AlignBottom | Qt::AlignHCenter, -90);
		custom_plot::patch::set_range(peak_annotation_rect, { -11.0, 1.0 }, { 0.0, 1.0 });

		constexpr int peak_height = 80;
		custom_plot::patch::set_fixed_size(peak_annotation_rect, anno_width, peak_height);
		custom_plot::patch::set_fixed_height(peak_rect, peak_height);

		peak_annotation_rect->setMarginGroup(QCP::msLeft | QCP::msRight, first_column_margin);
		//------------------------------//

		peak_rect->setMarginGroup(QCP::msLeft | QCP::msRight, second_column_margin);
		custom_plot::patch::remove_left_axis(peak_rect);
		custom_plot::patch::remove_bottom_subticks(peak_rect);
		custom_plot::patch::remove_bottom_grid(peak_rect);
		peak_rect->axis(QCPAxis::atBottom)->setNumberFormat("f");
		peak_rect->axis(QCPAxis::atBottom)->setNumberPrecision(0);
		peak_rect->axis(QCPAxis::atBottom)->setBasePen(annotation_axis_pen);
		peak_rect->axis(QCPAxis::atBottom)->setTickPen(annotation_axis_pen);

		int peak_rect_height = std::ceil(custom_plot::utility::get_max_text_height({ sequence_name }, gs.get_bottom_title_font()) * 2);
		custom_plot::patch::set_fixed_height(peak_rect, peak_rect_height);

		custom_plot::patch::set_range(peak_rect, QCPRange(region_start, region_end + 1), QCPRange(0, 10));
		for (auto&& p : res.peak_locations) {

			custom_plot::patch::rectangle_borderless(draw_area, peak_rect, p.first, 2.5, p.second - p.first, 5, Qt::black);

		}

		if (draw_ccan) {
			for (auto&& p : res.ccan_locations) {
				custom_plot::patch::rectangle_borderless(draw_area, peak_rect, p.first, 2.5, p.second - p.first, 5, custom_plot::color::firebrick3);
			}
		}

		custom_plot::set_bottom_title(peak_rect, sequence_name, gs);

		++layer;

		//-------------------//
		if (res.draw_legend) {
			QCPLayoutGrid* legend_layout = new QCPLayoutGrid;

			QCPAxisRect* t1 = new QCPAxisRect(draw_area);
			custom_plot::patch::remove_left_bottom_axis(t1);
			QCPAxisRect* t2 = new QCPAxisRect(draw_area);
			custom_plot::patch::remove_left_bottom_axis(t2);

			legend_layout->addElement(0, 0, t1);

			QCPAxisRect* r1 = new QCPAxisRect(draw_area);
			custom_plot::patch::set_range(r1, QCPRange(-1, 1), QCPRange(0, 1));
			custom_plot::patch::add_label(draw_area, r1, "CDS", 0, 0.5, QFont("Arial", 16), Qt::AlignVCenter | Qt::AlignLeft);
			custom_plot::patch::shape(draw_area, r1, QVector<double>{-1, -1, 0, 0}, QVector<double>{0, 1, 1, 0}, custom_plot::color::firebrick3, custom_plot::color::firebrick3, 2);
			custom_plot::patch::set_fixed_size(r1, 120, 40);
			custom_plot::patch::remove_left_bottom_axis(r1);
			legend_layout->addElement(0, 1, r1);

			QCPAxisRect* r2 = new QCPAxisRect(draw_area);
			custom_plot::patch::set_range(r2, QCPRange(-1, 1), QCPRange(0, 1));
			custom_plot::patch::add_label(draw_area, r2, "UTR", 0, 0.5, QFont("Arial", 16), Qt::AlignVCenter | Qt::AlignLeft);
			custom_plot::patch::shape(draw_area, r2, QVector<double>{-1, -1, 0, 0}, QVector<double>{0, 1, 1, 0}, Qt::black, Qt::white, 2);
			custom_plot::patch::set_fixed_size(r2, 120, 40);
			custom_plot::patch::remove_left_bottom_axis(r2);
			legend_layout->addElement(0, 2, r2);

			QCPAxisRect* r3 = new QCPAxisRect(draw_area);
			custom_plot::patch::set_range(r3, QCPRange(-1, 1), QCPRange(0, 1));
			custom_plot::patch::add_label(draw_area, r3, "Exon", 0, 0.5, QFont("Arial", 16), Qt::AlignVCenter | Qt::AlignLeft);
			custom_plot::patch::shape(draw_area, r3, QVector<double>{-1, -1, 0, 0}, QVector<double>{0, 1, 1, 0}, Qt::black, Qt::black, 2);
			custom_plot::patch::set_fixed_size(r3, 120, 40);
			custom_plot::patch::remove_left_bottom_axis(r3);
			legend_layout->addElement(0, 3, r3);

			QCPAxisRect* r4 = new QCPAxisRect(draw_area);
			custom_plot::patch::set_range(r4, QCPRange(-1, 1), QCPRange(0, 1));
			custom_plot::patch::add_label(draw_area, r4, "Gap", 0, 0.5, QFont("Arial", 16), Qt::AlignVCenter | Qt::AlignLeft);
			custom_plot::patch::shape_borderless(draw_area, r4, QVector<double>{-1, -1, 0, 0}, QVector<double>{0.4, 0.6, 0.6, 0.4}, Qt::black);
			custom_plot::patch::set_fixed_size(r4, 120, 40);
			custom_plot::patch::remove_left_bottom_axis(r4);
			legend_layout->addElement(0, 4, r4);

			legend_layout->addElement(0, 5, t2);

			legend_layout->setMarginGroup(QCP::msLeft | QCP::msRight, second_column_margin);

			main_layout->addElement(layer, 1, legend_layout);
		}

		//--------------------------------//

		custom_plot::add_title(draw_area, res.plot_title, gs);

		return draw_area;
	};

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> embedding_single_color_plot(
		const QString& title,
		const Eigen::MatrixXd& embedding_matrix,
		const QColor& color,
		const QStringList& embedding_names,
		const GraphSettings& gs
	) {		
		Eigen::ArrayXd x = embedding_matrix.col(0), y = embedding_matrix.col(1);
		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::set_scatter_plot_axis_style(
			draw_area,
			axis_rect,
			embedding_names[0],
			embedding_names[1],
			x,
			y,
			gs
		);

		custom_plot::patch::scatter(
			draw_area,
			axis_rect,
			x,
			y,
			color,
			gs.get_scatter_point_size()
		);

		custom_plot::add_title(draw_area, title, gs);

		return std::make_tuple(draw_area, axis_rect, legend_layout);
	};

	QCustomPlot* initialize_plot(const GraphSettings& gs) {
		QCustomPlot* draw_area = new QCustomPlot();
		draw_area->plotLayout()->clear();

		if (gs.is_transparent_background()) {
			draw_area->setBackground(Qt::transparent);
		}

		return draw_area;
	};

	std::pair<QCustomPlot*, QCPAxisRect*> prepare_ar(const GraphSettings& gs) {

		QCustomPlot* draw_area = new QCustomPlot();
		draw_area->plotLayout()->clear();

		QCPAxisRect* axis_rect = custom_plot::patch::new_axis_rect(draw_area);

		draw_area->plotLayout()->addElement(0, 0, axis_rect);

		if (gs.is_transparent_background()) {
			draw_area->setBackground(Qt::transparent);
		}

		return std::make_pair(draw_area, axis_rect);
	};

	std::pair<QCustomPlot*, QCPLayoutGrid*> prepare_lg(const GraphSettings& gs) {

		QCustomPlot* draw_area = new QCustomPlot();
		draw_area->plotLayout()->clear();

		QCPLayoutGrid* layout = new QCPLayoutGrid;
		draw_area->plotLayout()->addElement(0, 0, layout);

		if (gs.is_transparent_background()) {
			draw_area->setBackground(Qt::transparent);
		}

		return std::make_pair(draw_area, layout);
	};

	std::tuple<QCustomPlot*, QCPLayoutGrid*, QCPLayoutGrid*> prepare_lg_lg(const GraphSettings& gs) {
		QCustomPlot* draw_area = new QCustomPlot();
		draw_area->plotLayout()->clear();

		QCPLayoutGrid* layout = new QCPLayoutGrid;

		QCPLayoutGrid* left_layout = new QCPLayoutGrid;
		QCPLayoutGrid* right_layout = new QCPLayoutGrid;

		draw_area->plotLayout()->addElement(0, 0, layout);
		layout->addElement(0, 0, left_layout);
		layout->addElement(0, 1, right_layout);

		auto legend_layout = custom_plot::patch::set_legend_layout(draw_area, right_layout);

		if (gs.is_transparent_background()) {
			draw_area->setBackground(Qt::transparent);
		}

		return std::make_tuple(draw_area, left_layout, legend_layout);
	};

	std::tuple<QCustomPlot*, QCPAxisRect*, QCPLayoutGrid*> prepare(const GraphSettings& gs) {
		QCustomPlot* draw_area = new QCustomPlot();
		draw_area->plotLayout()->clear();

		QCPAxisRect* axis_rect = custom_plot::patch::new_axis_rect(draw_area);

		QCPLayoutGrid* layout = new QCPLayoutGrid;
		QCPLayoutGrid* right_layout = new QCPLayoutGrid;

		draw_area->plotLayout()->addElement(0, 0, layout);
		layout->addElement(0, 0, axis_rect);
		layout->addElement(0, 1, right_layout);

		auto legend_layout = custom_plot::patch::set_legend_layout(draw_area, right_layout);

		if (gs.is_transparent_background()) {
			draw_area->setBackground(Qt::transparent);
		}

		return std::make_tuple(draw_area, axis_rect, legend_layout);
	};

	void set_axis_label(
		QCPAxis::AxisType type,
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	) {
		custom_plot::patch::set_axis_label(
			type, 
			axis_rect, 
			location, 
			labels, 
			tick_length, 
			gs.get_label_font(type),
			gs.get_label_angle(type));		
	};

	void set_left_axis_label(
		QCPAxisRect* axis_rect, 
		const Eigen::ArrayXd& location, 
		const QStringList& labels, 
		int tick_length,
		const GraphSettings& gs
	) {
		set_axis_label(QCPAxis::atLeft, axis_rect, location, labels, tick_length, gs);
	};

	void set_top_axis_label(
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	) {
		set_axis_label(QCPAxis::atTop, axis_rect, location, labels, tick_length, gs);
	};

	void set_bottom_axis_label(
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& location,
		const QStringList& labels,
		int tick_length,
		const GraphSettings& gs
	) {
		set_axis_label(QCPAxis::atBottom, axis_rect, location, labels, tick_length, gs);
	};

	void set_left_title(QCPAxisRect* axis_rect, const QString& label, const GraphSettings& gs, bool use_setting_title) {
		axis_rect->axis(QCPAxis::atLeft)->setLabel(use_setting_title ? gs.get_left_title(label) : label);
		axis_rect->axis(QCPAxis::atLeft)->setLabelFont(gs.get_left_title_font());
	};

	void set_bottom_title(QCPAxisRect* axis_rect, const QString& label, const GraphSettings& gs, bool use_setting_title) {
		axis_rect->axis(QCPAxis::atBottom)->setLabel(use_setting_title ? gs.get_bottom_title(label) : label);
		axis_rect->axis(QCPAxis::atBottom)->setLabelFont(gs.get_bottom_title_font());
	};

	void add_title(QCustomPlot* draw_area, const QString& title, const GraphSettings& gs) {
		draw_area->plotLayout()->insertRow(0);
		auto* text = new SoapTextElement(draw_area, gs.get_title(title), gs.get_title_font());
		text->setMargins(QMargins(0, 10, 0, 0));
		draw_area->plotLayout()->addElement(0, 0, text);
	};


	void set_arrow_axis(
		QCustomPlot* draw_area, 
		QCPAxisRect* axis_rect, 
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const QString& bottom_title, 
		const QString& left_title, 
		double margin, 
		int type,
		const GraphSettings& gs
	) {

		custom_plot::patch::remove_left_bottom_axis(axis_rect);

		double maximum_x = x.maxCoeff(), minimum_x = x.minCoeff(), maximum_y = y.maxCoeff(), minimum_y = y.minCoeff();
		double x_span = maximum_x - minimum_x, y_span = maximum_y - minimum_y;

		QCPItemLine* arrow_x = new QCPItemLine(draw_area);
		arrow_x->start->setCoords(minimum_x - margin / 2 * x_span, minimum_y - margin / 2 * y_span);
		arrow_x->end->setCoords(minimum_x + margin * 1.5 * x_span, minimum_y - margin / 2 * y_span);
		arrow_x->setHead(QCPLineEnding::esSpikeArrow);

		QCPItemLine* arrow_y = new QCPItemLine(draw_area);
		arrow_y->start->setCoords(minimum_x - margin / 2 * x_span, minimum_y - margin / 2 * y_span);
		arrow_y->end->setCoords(minimum_x - margin / 2 * x_span, minimum_y + margin * 1.5 * y_span);
		arrow_y->setHead(QCPLineEnding::esSpikeArrow);

		QCPItemText* label_x = new QCPItemText(draw_area);
		label_x->setClipToAxisRect(false);
		label_x->position->setAxisRect(axis_rect);
		label_x->position->setAxes(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		label_x->position->setType(QCPItemPosition::ptPlotCoords);
		if (type == 0) {
			label_x->setPositionAlignment(Qt::AlignCenter);
			label_x->position->setCoords(minimum_x + margin / 2 * x_span, minimum_y - margin * 0.75 * y_span);
		}
		else {
			label_x->setPositionAlignment(Qt::AlignLeft | Qt::AlignVCenter);
			label_x->position->setCoords(minimum_x + margin * 1.5 * x_span, minimum_y - margin / 2 * y_span);
		}
		label_x->setText(gs.get_bottom_title(bottom_title));
		label_x->setFont(gs.get_bottom_title_font());

		QCPItemText* label_y = new QCPItemText(draw_area);
		label_y->setClipToAxisRect(false);
		label_y->position->setAxisRect(axis_rect);
		label_y->position->setAxes(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
		label_y->position->setType(QCPItemPosition::ptPlotCoords);
		if (type == 0) {
			label_y->setPositionAlignment(Qt::AlignCenter);
			label_y->position->setCoords(minimum_x - margin * 0.75 * x_span, minimum_y + margin / 2 * y_span);
		}
		else {
			label_y->setPositionAlignment(Qt::AlignLeft | Qt::AlignVCenter);
			label_y->position->setCoords(minimum_x - margin / 2 * x_span, minimum_y + margin * 1.5 * y_span);
		}
		label_y->setRotation(-90);
		label_y->setText(gs.get_left_title(left_title));
		label_y->setFont(gs.get_left_title_font());
	};

	void set_scatter_plot_axis_style(
		QCustomPlot* draw_area, 
		QCPAxisRect* axis_rect, 
		const QString& bottom_title,
		const QString& left_title,
		const Eigen::ArrayXd& x,
		const Eigen::ArrayXd& y,
		const GraphSettings& gs) {

		custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(x, y));

		auto axis_style = gs.get_axis_style();

		if (axis_style == GraphSettings::AxisStyle::Arrow) {
			custom_plot::set_arrow_axis(draw_area, axis_rect, x, y, bottom_title, left_title, 0.1, 0, gs);
		}
		else if (axis_style == GraphSettings::AxisStyle::Arrow2) {
			custom_plot::set_arrow_axis(draw_area, axis_rect, x, y, bottom_title, left_title, 0.1, 1, gs);
		}
		else if (axis_style == GraphSettings::AxisStyle::NoAxis) {
			custom_plot::patch::remove_left_bottom_axis(axis_rect);
		}
		else {
			custom_plot::set_simple_axis(axis_rect, bottom_title, left_title, gs);
		}
	};
	
	void bar_plot(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		QColor color,
		const Eigen::ArrayXd& bar_location,
		const Eigen::ArrayXd& bar_length,
		int bar_width,
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs,
		bool vertical
	) {
		custom_plot::set_simple_axis(axis_rect, bottom_title, left_title, gs);

		custom_plot::patch::bar(draw_area, axis_rect, color, bar_location, bar_length, bar_width, vertical);

		int ymax = std::max(bar_length.maxCoeff() * 1.1, 0.0);
		int ymin = std::min(bar_length.minCoeff() * 1.1, 0.0);

		custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(bar_location), custom_plot::utility::get_range(ymin, ymax));
	};

	void bar_plot_enrichment(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		const Eigen::ArrayXd& bar_location,
		const Eigen::ArrayXd& bar_length,
		const Eigen::ArrayXd& values,
		int bar_width,
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs,
		bool vertical 
	) {

		custom_plot::patch::bar_gradient(
			draw_area,
			axis_rect,
			bar_location,
			bar_length,
			values,
			gs.get_gradient_low_color(),
			gs.get_gradient_middle_color(),
			gs.get_gradient_high_color(),
			bar_width,
			vertical
		);

		custom_plot::patch::remove_left_axis(axis_rect);

		QPen pen(Qt::black);
		pen.setWidth(3);
		axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
		axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
		axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
		custom_plot::set_bottom_title(axis_rect, "Gene Number", gs, true);

		int ymax = std::max(bar_length.maxCoeff() * 1.1, 0.0);
		int ymin = std::min(bar_length.minCoeff() * 1.1, 0.0);

		custom_plot::patch::set_range(axis_rect, custom_plot::utility::get_range(ymin, ymax), custom_plot::utility::get_range(bar_location));
	}

	void histogram_plot(
		QCustomPlot* draw_area, 
		QCPAxisRect* axis_rect, 
		const QVector<double>& x, 
		int unit, 
		QColor color, 
		int width,
		const QString& bottom_title,
		const QString& left_title,
		const GraphSettings& gs) 
	{
		custom_plot::set_simple_axis(axis_rect, bottom_title, left_title, gs);

		auto [counts, _, centers] = custom::histogram(custom::cast<Eigen::ArrayX>(x), unit);

		double minimum_x = 3 * centers[0] - 2 * centers[1], maximum_x = 3 * centers[unit - 1] - 2 * centers[unit - 2], maximum_y = counts.maxCoeff() * 1.1;

		custom_plot::patch::bar(
			draw_area,
			axis_rect,
			color,
			centers,
			counts.cast<double>(),
			width,
			true
		);

		custom_plot::patch::set_range(axis_rect, { minimum_x, maximum_x }, { 0, maximum_y });
	};

	void heatmap_plot2(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		int cell_width,
		int cell_height,
		int border_width,
		const QStringList& left_labels,
		const QStringList& bottom_labels,
		const Eigen::MatrixXd& values,
		const GraphSettings& gs
	) {
		double minimum_value = values.minCoeff(), 
			maximum_value = values.maxCoeff(), 
			span = (maximum_value - minimum_value) / 2, 
			middle = (minimum_value + maximum_value) / 2;
		int nrow = left_labels.size(), ncol = bottom_labels.size();

		custom_plot::patch::remove_left_bottom_axis(axis_rect);

		int all_width = ncol * cell_width + (ncol - 1) * border_width;
		int all_height = nrow * cell_height + (nrow - 1) * border_width;

		custom_plot::patch::set_range(axis_rect, QCPRange(0, all_width), QCPRange(0, all_height + 10));

		QColor low_color = gs.get_gradient_low_color(), 
			high_color = gs.get_gradient_high_color(),
			middle_color = gs.get_gradient_middle_color();
		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {

				double val = values(i, j);
				QColor this_color;
				if (val < middle) {
					this_color = custom_plot::utility::gradient_color((val - minimum_value) / span, low_color, middle_color);
				}
				else {
					this_color = custom_plot::utility::gradient_color((val - middle) / span, middle_color, high_color);
				}

				custom_plot::patch::shape_borderless(
					draw_area,
					axis_rect,
					QVector<double>{ double(j* (cell_width + border_width)),
					double(j* (cell_width + border_width) + cell_width),
					double(j* (cell_width + border_width) + cell_width),
					double(j* (cell_width + border_width)) },
					QVector<double>{ double(i* (cell_height + border_width)),
					double(i* (cell_height + border_width)),
					double(i* (cell_height + border_width) + cell_height),
					double(i* (cell_height + border_width) + cell_height)	},
					this_color
				);
			}
		}

		custom_plot::set_left_axis_label(
			axis_rect,
			Eigen::ArrayXd::LinSpaced(nrow, 0.5 * cell_height, 0.5 * cell_height + (nrow - 1) * (cell_height + border_width)),
			left_labels,
			0,
			gs
		);

		custom_plot::set_bottom_axis_label(
			axis_rect,
			Eigen::ArrayXd::LinSpaced(ncol, 0.5 * cell_width, 0.5 * cell_width + (ncol - 1) * (cell_width + border_width)),
			bottom_labels,
			0,
			gs
		);
	};

	void set_heatmap_color(
		QCustomPlot* draw_area,
		const GraphSettings& gs, 
		QCPColorMap* heatmap, 
		double min_val, 
		double max_val) {

		QCPColorGradient gradient;
		gradient.setColorStopAt(0.5, gs.get_gradient_middle_color());
		gradient.setColorStopAt(1.0, gs.get_gradient_high_color());
		gradient.setColorStopAt(0.0, gs.get_gradient_low_color());

		heatmap->setGradient(gradient);
		heatmap->setDataRange({ min_val, max_val });
	};

	QCustomPlot* heatmap_plot(
		const HEATMAP_PLOT_ELEMENTS& ele, 
		const GraphSettings& gs) {

		auto [draw_area, main_layout, legend_layout] = custom_plot::prepare_lg_lg(gs);

		QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

		QCPAxisRect* heatmap_rect = new QCPAxisRect(draw_area);
		main_layout->addElement(0, 0, heatmap_rect);

		heatmap_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
		heatmap_rect->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);

		double minimum_value = ele.mat.minCoeff();
		double maximum_value = ele.mat.maxCoeff();
		double span = (maximum_value - minimum_value) / 2;
		double middle = (minimum_value + maximum_value) / 2;
		int nrow = ele.mat.rows(), ncol = ele.mat.cols();

		custom_plot::patch::remove_all_axis(heatmap_rect);

		QCPColorMap* heatmap = new QCPColorMap(heatmap_rect->axis(QCPAxis::atBottom), heatmap_rect->axis(QCPAxis::atLeft));
		heatmap->data()->setSize(ncol, nrow);
		if (nrow > 1) {
			heatmap->data()->setRange(QCPRange(0, ncol - 1), QCPRange(0, nrow - 1));
		}
		else {
			heatmap->data()->setRange(QCPRange(0, ncol - 1), QCPRange(-1, 1));
		}
		custom_plot::patch::set_range(heatmap_rect, QCPRange(-0.5, ncol - 0.5), QCPRange(-0.5, nrow - 0.5));

		for (int i = 0; i < ncol; ++i) {
			for (int j = 0; j < nrow; ++j) {
				heatmap->data()->setCell(i, j, ele.mat(j, i));
			}
		}
		heatmap->setInterpolate(false);
		heatmap->setTightBoundary(false);
		custom_plot::set_heatmap_color(draw_area, gs, heatmap, minimum_value, maximum_value);

		heatmap_rect->setupFullAxesBox();

		if (ele.show_column_names) {

			QCPAxis::AxisType type = ele.show_column_names_at_bottom ? QCPAxis::atBottom : QCPAxis::atTop;

			custom_plot::set_axis_label(
				type,
				heatmap_rect,
				Eigen::ArrayXd::LinSpaced(ncol, 0, ncol - 1),
				ele.column_names,
				0,
				gs
			);
		}

		if (ele.show_row_names) {

			QCPAxis::AxisType type = ele.show_row_names_at_left ? QCPAxis::atLeft : QCPAxis::atRight;

			custom_plot::set_axis_label(
				type,
				heatmap_rect,
				Eigen::ArrayXd::LinSpaced(nrow, 0, nrow - 1),
				ele.row_names,
				0,
				gs
			);

		}

		if (ele.info.contains("Legend Low Label") && ele.info.contains("Legend High Label")) {

			custom_plot::patch::add_gradient_legend(
				draw_area,
				legend_layout,
				minimum_value,
				maximum_value,
				gs.get_legend_title(ele.legend_title),
				ele.info["Legend Low Label"],
				ele.info["Legend High Label"],
				gs.get_legend_title_font(),
				gs.get_legend_label_font(),
				gs.get_gradient_low_color(),
				gs.get_gradient_middle_color(),
				gs.get_gradient_high_color()
			);
		}
		else {
			custom_plot::add_gradient_legend(
				draw_area,
				legend_layout,
				minimum_value,
				maximum_value,
				ele.legend_title,
				gs);
		}

		if (!ele.column_annotations.isEmpty()) {
			int n_valid_annotation = std::ranges::count(custom::sapply(ele.column_annotations, [ncol](auto&& p) {
				return p.second.size() == ncol;
			}), true);

			if (n_valid_annotation > 0) {

				QCPAxisRect* column_annotation_rect = new QCPAxisRect(draw_area);

				column_annotation_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
				
				custom_plot::patch::remove_all_axis(column_annotation_rect);

				if (!ele.annotate_column_at_top) {
					main_layout->addElement(1, 0, column_annotation_rect);

					main_layout->setRowStretchFactor(0, 15);
					main_layout->setRowStretchFactor(1, n_valid_annotation);
				}
				else {
					main_layout->insertRow(0);
					main_layout->addElement(0, 0, column_annotation_rect);

					main_layout->setRowStretchFactor(1, 15);
					main_layout->setRowStretchFactor(0, n_valid_annotation);
				}													

				int count{ 0 };
				QStringList column_annotation_names;
				QList<QColor> all_colors;
				QList<QMap<QString, QColor>> color_maps;
				std::ranges::for_each(ele.column_annotations, [&](auto&& p) {
					if (p.second.size() != ncol) {
						return;
					}
				auto levels = custom::unique(p.second);
				auto color_map = gs.palette_map(levels);
				auto colors = color_map.values();

				all_colors << colors;
				color_maps << color_map;
				column_annotation_names << p.first;

				if (ele.show_annotation_legend) {

					custom_plot::add_round_legend(
						draw_area,
						legend_layout,
						levels,
						colors,
						p.first,
						gs
					);
				}
				++count;
				});

				all_colors = custom::unique(all_colors);
				QMap<QColor, double> color_value_map;
				double inc = 1.0 / all_colors.size();
				int n_color = all_colors.size();
				for (int i = 0; i < n_color; ++i) {
					color_value_map[all_colors[i]] = i * inc;
				}

				QCPColorMap* annotation_heatmap = new QCPColorMap(column_annotation_rect->axis(QCPAxis::atBottom), column_annotation_rect->axis(QCPAxis::atLeft));
				annotation_heatmap->data()->setSize(ncol, n_valid_annotation);
				if (n_valid_annotation > 1) {
					annotation_heatmap->data()->setRange(QCPRange(0, ncol - 1), QCPRange(0, n_valid_annotation - 1));
				}
				else {
					annotation_heatmap->data()->setRange(QCPRange(0, ncol - 1), QCPRange(-1, 1));
				}

				int column_annotation_height = n_valid_annotation;
				custom_plot::patch::set_range(column_annotation_rect, QCPRange(-0.5, ncol - 0.5), QCPRange(-0.5, column_annotation_height - 0.5));

				count = 0;
				std::ranges::for_each(ele.column_annotations,
					[&](auto&& p) {
					if (p.second.size() != ncol) {
						return;
					}
				for (int j = 0; j < ncol; ++j) {
					annotation_heatmap->data()->setCell(j, count, color_value_map[color_maps[count][p.second[j]]]);
				}
				++count;
				});
				annotation_heatmap->setInterpolate(false);
				annotation_heatmap->setTightBoundary(false);
				QCPColorGradient gradient;
				for (auto&& [c, val] : color_value_map.asKeyValueRange()) {
					gradient.setColorStopAt(val, c);
				}

				annotation_heatmap->setGradient(gradient);
				annotation_heatmap->setDataRange({ 0.0, 1.0 });
				
				if (ele.show_column_annotation_name) {

					column_annotation_rect->setupFullAxesBox();

					QCPAxis::AxisType type = ele.show_column_annotation_name_at_left ? QCPAxis::atLeft : QCPAxis::atRight;

					custom_plot::patch::set_axis_label(
						type,
						column_annotation_rect,
						Eigen::ArrayXd::LinSpaced(n_valid_annotation, 0, n_valid_annotation - 1),
						column_annotation_names,
						0,
						gs.get_left_label_font(),
						0
					);

				}
			}
		}

		if (!ele.row_annotations.isEmpty()) {
			int n_valid_annotation = std::ranges::count(custom::sapply(ele.row_annotations, [nrow](auto&& p) {
				return p.second.size() == nrow;
			}), true);

			if (n_valid_annotation > 0) {

				QCPAxisRect* row_annotation_rect = new QCPAxisRect(draw_area);

				row_annotation_rect->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);

				custom_plot::patch::remove_all_axis(row_annotation_rect);

				int row_anno_row = ele.annotate_column_at_top ? 1 : 0;
				if (ele.annotate_row_at_left) {
					main_layout->insertColumn(0);
					main_layout->addElement(1, row_anno_row, row_annotation_rect);

					main_layout->setColumnStretchFactor(1, 15);
					main_layout->setColumnStretchFactor(0, n_valid_annotation);
				}
				else {

					main_layout->addElement(0, row_anno_row, row_annotation_rect);

					main_layout->setColumnStretchFactor(0, 15);
					main_layout->setColumnStretchFactor(1, n_valid_annotation);
				}

				int count{ 0 };
				QStringList row_annotation_names;
				QList<QColor> all_colors;
				QList<QMap<QString, QColor>> color_maps;
				std::ranges::for_each(ele.row_annotations, [&](auto&& p) {
					if (p.second.size() != nrow) {
						return;
					}
				auto levels = custom::unique(p.second);
				auto color_map = gs.palette_map(levels);
				auto colors = color_map.values();

				all_colors << colors;
				color_maps << color_map;
				row_annotation_names << p.first;

				if (ele.show_annotation_legend) {

					custom_plot::add_round_legend(
						draw_area,
						legend_layout,
						levels,
						colors,
						p.first,
						gs
					);
				}
				++count;
				});

				all_colors = custom::unique(all_colors);
				QMap<QColor, double> color_value_map;
				double inc = 1.0 / all_colors.size();
				int n_color = all_colors.size();
				for (int i = 0; i < n_color; ++i) {
					color_value_map[all_colors[i]] = i * inc;
				}

				QCPColorMap* annotation_heatmap = new QCPColorMap(row_annotation_rect->axis(QCPAxis::atBottom), row_annotation_rect->axis(QCPAxis::atLeft));
				annotation_heatmap->data()->setSize(n_valid_annotation, nrow);
				if (n_valid_annotation > 1) {
					annotation_heatmap->data()->setRange(QCPRange(0, n_valid_annotation - 1), QCPRange(0, nrow - 1));
				}
				else {
					annotation_heatmap->data()->setRange(QCPRange(-1, 1), QCPRange(0, nrow - 1));
				}

				int row_annotation_height = nrow;
				custom_plot::patch::set_range(row_annotation_rect, QCPRange(-0.5, n_valid_annotation - 0.5), QCPRange(-0.5, row_annotation_height - 0.5));

				count = 0;
				std::ranges::for_each(ele.row_annotations,
					[&](auto&& p) {
					if (p.second.size() != nrow) {
						return;
					}
				for (int j = 0; j < nrow; ++j) {
					annotation_heatmap->data()->setCell(count, j, color_value_map[color_maps[count][p.second[j]]]);
				}
				++count;
				});
				annotation_heatmap->setInterpolate(false);
				annotation_heatmap->setTightBoundary(false);
				QCPColorGradient gradient;
				for (auto&& [c, val] : color_value_map.asKeyValueRange()) {
					gradient.setColorStopAt(val, c);
				}

				annotation_heatmap->setGradient(gradient);
				annotation_heatmap->setDataRange({ 0.0, 1.0 });

				if (ele.show_row_annotation_name) {

					row_annotation_rect->setupFullAxesBox();

					QCPAxis::AxisType type = ele.show_row_annotation_name_at_bottom ? QCPAxis::atBottom : QCPAxis::atTop;

					custom_plot::patch::set_axis_label(
						type,
						row_annotation_rect,
						Eigen::ArrayXd::LinSpaced(n_valid_annotation, 0, n_valid_annotation - 1),
						row_annotation_names,
						0,
						gs.get_left_label_font(),
						0
					);

				}
			}
		}

		return draw_area;
	};

	void heatmap_plot(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		QCPLayoutGrid* legend_layout,
		int cell_width,
		int cell_height,
		int border_width,
		const QStringList& left_labels,
		const QStringList& bottom_labels,
		const Eigen::MatrixXd& values,
		bool normalize_by_row,
		const QString& group_legend_title,
		bool annotation_at_top,
		const GraphSettings& gs
	) {
		double minimum_value = values.minCoeff(), maximum_value = values.maxCoeff(), span = (maximum_value - minimum_value) / 2, middle = (minimum_value + maximum_value) / 2;
		int nrow = left_labels.size(), ncol = bottom_labels.size();

		int annotation_height = annotation_at_top ? 0 : 10;

		custom_plot::patch::remove_left_bottom_axis(axis_rect);

		int all_width = ncol * cell_width + (ncol - 1) * border_width;
		int all_height = nrow * cell_height + (nrow - 1) * border_width;

		custom_plot::patch::set_range(axis_rect, QCPRange(0, all_width), QCPRange(0, all_height + 10));

		QColor low_color = gs.get_gradient_low_color(), high_color = gs.get_gradient_high_color(), middle_color = gs.get_gradient_middle_color();
		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {

				double val = values(i, j);
				QColor this_color;
				if (val < middle) {
					this_color = custom_plot::utility::gradient_color((val - minimum_value) / span, low_color, middle_color);
				}
				else {
					this_color = custom_plot::utility::gradient_color((val - middle) / span, middle_color, high_color);
				}

				custom_plot::patch::shape_borderless(
					draw_area,
					axis_rect,
					QVector<double>{ double(j* (cell_width + border_width)),
					double(j* (cell_width + border_width) + cell_width),
					double(j* (cell_width + border_width) + cell_width),
					double(j* (cell_width + border_width)) },
					QVector<double>{ double(i* (cell_height + border_width) + annotation_height),
					double(i* (cell_height + border_width) + annotation_height),
					double(i* (cell_height + border_width) + cell_height + annotation_height),
					double(i* (cell_height + border_width) + cell_height + annotation_height)	},
					this_color
				);
			}
		}

		auto annotation_colors = gs.palette(bottom_labels);

		double start_height = annotation_at_top ? nrow * (cell_height + border_width) - border_width + 5 : 0;

		for (int i = 0; i < ncol; ++i) {

			custom_plot::patch::shape_borderless(
				draw_area,
				axis_rect,
				QVector<double>{ double(i* (cell_width + border_width)),
				double(i* (cell_width + border_width) + cell_width),
				double(i* (cell_width + border_width) + cell_width),
				double(i* (cell_width + border_width)) },
				QVector<double>{ start_height,
				start_height,
				start_height + 5,
				start_height + 5 },
				annotation_colors[i]
			);
		}

		custom_plot::set_left_axis_label(
			axis_rect, 
			Eigen::ArrayXd::LinSpaced(nrow, annotation_height + 0.5 * cell_height , annotation_height + 0.5 * cell_height + (nrow - 1) * (cell_height + border_width)),
			left_labels, 
			0,
			gs
		);

		if (annotation_at_top) {

			axis_rect->setupFullAxesBox();
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atRight);
			custom_plot::patch::remove_axis(axis_rect, QCPAxis::atTop);
			
			axis_rect->axis(QCPAxis::atTop)->setRange(QCPRange(0, all_width));

			custom_plot::set_top_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(ncol, 0.5 * cell_width, 0.5 * cell_width + (ncol - 1) * (cell_width + border_width)),
				bottom_labels,
				0,
				gs
			);
		}
		else {
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(ncol, 0.5 * cell_width, 0.5 * cell_width + (ncol - 1) * (cell_width + border_width)),
				bottom_labels,
				0,
				gs
			);
		}

		if (normalize_by_row) {
			custom_plot::patch::add_gradient_legend(
				draw_area,
				legend_layout,
				minimum_value,
				maximum_value,
				gs.get_legend_title("Normalized Expression"),
				"Low",
				"High",
				gs.get_legend_title_font(),
				gs.get_legend_label_font(),
				gs.get_gradient_low_color(),
				gs.get_gradient_middle_color(),
				gs.get_gradient_high_color()
			);
		}
		else {
			custom_plot::add_gradient_legend(
				draw_area,
				legend_layout,
				minimum_value,
				maximum_value,
				"Expression",
				gs);
		}
		custom_plot::add_round_legend(
			draw_area,
			legend_layout,
			bottom_labels,
			annotation_colors,
			group_legend_title,
			gs
		);
	};

	void bubble_plot(
		QCustomPlot* draw_area,
		QCPAxisRect* axis_rect,
		QCPLayoutGrid* legend_layout,
		const QStringList& left_labels,
		const QStringList& bottom_labels,
		const Eigen::MatrixXd& values,
		const Eigen::MatrixXi& dot_size,
		const QString& value_legend_title,
		const QString& size_legend_title,
		const GraphSettings& gs
	) {
		double minimum_value = values.minCoeff(), maximum_value = values.maxCoeff(), span = (maximum_value - minimum_value) / 2, middle = (minimum_value + maximum_value) / 2;
		int nrow = left_labels.size(), ncol = bottom_labels.size();

		custom_plot::patch::set_border_only(axis_rect, Qt::black, 3);

		custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * ncol), QCPRange(0, 2 * nrow));
		QColor low_color = gs.get_gradient_low_color(), high_color = gs.get_gradient_high_color(), middle_color = gs.get_gradient_middle_color();
		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {
				double val = values(i, j);
				if (val < middle) {
					custom_plot::patch::dot(draw_area, axis_rect, j * 2 + 1, i * 2 + 1, dot_size(i, j), custom_plot::utility::gradient_color((val - minimum_value) / span, low_color, middle_color));
				}
				else {
					custom_plot::patch::dot(draw_area, axis_rect, j * 2 + 1, i * 2 + 1, dot_size(i, j), custom_plot::utility::gradient_color((val - middle) / span, middle_color, high_color));
				}
			}
		}

		custom_plot::set_left_axis_label(axis_rect, Eigen::ArrayXd::LinSpaced(nrow, 1, 2 * nrow - 1), left_labels, 6, gs);
		custom_plot::set_bottom_axis_label(axis_rect, Eigen::ArrayXd::LinSpaced(ncol, 1, 2 * ncol - 1), bottom_labels, 6, gs);
		custom_plot::proportion_legend(draw_area, legend_layout, size_legend_title, gs);
		custom_plot::add_gradient_legend(draw_area, legend_layout, minimum_value, maximum_value, value_legend_title, gs);
	};

	void set_simple_axis(
		QCPAxisRect* axis_rect, 
		const QString& bottom_title, 
		const QString& left_title, 
		const GraphSettings& gs,
		bool use_setting_title
	) {
		QPen pen(Qt::black);
		pen.setWidth(3);
		axis_rect->axis(QCPAxis::atLeft)->setBasePen(pen);
		axis_rect->axis(QCPAxis::atLeft)->setTickPen(pen);
		axis_rect->axis(QCPAxis::atLeft)->setSubTicks(false);

		axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 6);
		axis_rect->axis(QCPAxis::atLeft)->setTickLabelFont(gs.get_left_label_font());

		auto left_label_rotation = gs.get_left_label_angle();
		if (left_label_rotation != 0) {
			axis_rect->axis(QCPAxis::atLeft)->setTickLabelRotation(left_label_rotation);
		}

		axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
		axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
		axis_rect->axis(QCPAxis::atBottom)->setSubTicks(false);

		axis_rect->axis(QCPAxis::atBottom)->setTickLength(0, 6);
		axis_rect->axis(QCPAxis::atBottom)->setTickLabelFont(gs.get_bottom_label_font());

		auto bottom_label_rotation = gs.get_bottom_label_angle();
		if (bottom_label_rotation != 0) {
			axis_rect->axis(QCPAxis::atBottom)->setTickLabelRotation(bottom_label_rotation);
		}

		axis_rect->axis(QCPAxis::atBottom)->setUpperEnding(QCPLineEnding::esSpikeArrow);
		axis_rect->axis(QCPAxis::atLeft)->setUpperEnding(QCPLineEnding::esSpikeArrow);

		axis_rect->axis(QCPAxis::atBottom)->setLabel(use_setting_title ? gs.get_bottom_title(bottom_title) : bottom_title);
		axis_rect->axis(QCPAxis::atBottom)->setLabelFont(gs.get_bottom_title_font());

		axis_rect->axis(QCPAxis::atLeft)->setLabel(use_setting_title ? gs.get_left_title(left_title) : left_title);
		axis_rect->axis(QCPAxis::atLeft)->setLabelFont(gs.get_left_title_font());

	};

	void set_simple_axis_no_title(
		QCPAxisRect* axis_rect,
		const GraphSettings& gs
	) {
		QPen pen(Qt::black);
		pen.setWidth(3);
		axis_rect->axis(QCPAxis::atLeft)->setBasePen(pen);
		axis_rect->axis(QCPAxis::atLeft)->setTickPen(pen);
		axis_rect->axis(QCPAxis::atLeft)->setSubTicks(false);

		axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 6);
		axis_rect->axis(QCPAxis::atLeft)->setTickLabelFont(gs.get_left_label_font());

		axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
		axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
		axis_rect->axis(QCPAxis::atBottom)->setSubTicks(false);

		axis_rect->axis(QCPAxis::atBottom)->setTickLength(0, 6);
		axis_rect->axis(QCPAxis::atBottom)->setTickLabelFont(gs.get_bottom_label_font());

		axis_rect->axis(QCPAxis::atBottom)->setUpperEnding(QCPLineEnding::esSpikeArrow);
		axis_rect->axis(QCPAxis::atLeft)->setUpperEnding(QCPLineEnding::esSpikeArrow);
	};

	void add_gradient_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		double minimum_value,
		double maximum_value,
		const QString& title,
		const GraphSettings& gs,
		QColor low_color,
		QColor middle_color,
		QColor high_color
	) {
		QCPLayoutGrid* sub_legend_layout = new QCPLayoutGrid;
		QCPLayoutGrid* color_scale_layout = new QCPLayoutGrid;
		QCPAxisRect* color_scale_left_rect = new QCPAxisRect(draw_area, false);
		QCPAxisRect* color_scale_right_rect = new QCPAxisRect(draw_area, false);
		sub_legend_layout->setRowSpacing(30);
		auto [row, col] = custom_plot::patch::find_next_empty_position(legend_layout);
		legend_layout->addElement(row, col, sub_legend_layout);

		QCPColorGradient gradient;
		gradient.setColorStopAt(0.0, low_color);
		gradient.setColorStopAt(0.5, middle_color);
		gradient.setColorStopAt(1.0, high_color);
		QCPColorScale* color_scale = new QCPColorScale(draw_area);
		color_scale->setType(QCPAxis::atRight);

		if (minimum_value == maximum_value) {
			maximum_value += 0.01;
		}
		color_scale->setDataRange(QCPRange(minimum_value, maximum_value));
		color_scale->setGradient(gradient);
		color_scale->axis()->setLabelFont(gs.get_legend_tick_label_font());

		color_scale->mAxisRect->axis(QCPAxis::atBottom)->setBasePen(Qt::NoPen);
		color_scale->mAxisRect->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
		color_scale->mAxisRect->axis(QCPAxis::atTop)->setBasePen(Qt::NoPen);
		color_scale->mAxisRect->axis(QCPAxis::atRight)->setBasePen(Qt::NoPen);

		if (!gs.is_legend_tick_shown()) {
			color_scale->axis()->setTickPen(Qt::NoPen);
			color_scale->axis()->setSubTickPen(Qt::NoPen);
		}

		color_scale_layout->addElement(0, 0, color_scale_left_rect);
		color_scale_layout->addElement(0, 1, color_scale);
		color_scale_layout->addElement(0, 2, color_scale_right_rect);
		color_scale->setBarWidth(28);
		color_scale_layout->setMinimumSize(200, 200);
		color_scale_layout->setMaximumSize(200, 200);
		sub_legend_layout->addElement(0, 0, color_scale_layout);
		custom_plot::patch::add_title(draw_area, sub_legend_layout, gs.get_legend_title(title), gs.get_legend_title_font());
	};

	void add_gradient_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		double minimum_value,
		double maximum_value,
		const QString& title,
		const GraphSettings& gs
	) {
		QCPLayoutGrid* sub_legend_layout = new QCPLayoutGrid;
		QCPLayoutGrid* color_scale_layout = new QCPLayoutGrid;
		QCPAxisRect* color_scale_left_rect = new QCPAxisRect(draw_area, false);
		QCPAxisRect* color_scale_right_rect = new QCPAxisRect(draw_area, false);
		sub_legend_layout->setRowSpacing(30);
		auto [row, col] = custom_plot::patch::find_next_empty_position(legend_layout);
		legend_layout->addElement(row, col, sub_legend_layout);

		QCPColorGradient gradient;
		gradient.setColorStopAt(0.0, gs.get_gradient_low_color());
		gradient.setColorStopAt(0.5, gs.get_gradient_middle_color());
		gradient.setColorStopAt(1.0, gs.get_gradient_high_color());
		QCPColorScale* color_scale = new QCPColorScale(draw_area);
		color_scale->setType(QCPAxis::atRight);

		color_scale->setDataRange(QCPRange(minimum_value, maximum_value));
		color_scale->setGradient(gradient);
		color_scale->axis()->setLabelFont(gs.get_legend_tick_label_font());

		color_scale->mAxisRect->axis(QCPAxis::atBottom)->setBasePen(Qt::NoPen);
		color_scale->mAxisRect->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
		color_scale->mAxisRect->axis(QCPAxis::atTop)->setBasePen(Qt::NoPen);
		color_scale->mAxisRect->axis(QCPAxis::atRight)->setBasePen(Qt::NoPen);

		if (!gs.is_legend_tick_shown()) {
			color_scale->axis()->setTickPen(Qt::NoPen);
			color_scale->axis()->setSubTickPen(Qt::NoPen);
		}

		color_scale_layout->addElement(0, 0, color_scale_left_rect);
		color_scale_layout->addElement(0, 1, color_scale);
		color_scale_layout->addElement(0, 2, color_scale_right_rect);
		color_scale->setBarWidth(28);
		color_scale_layout->setMinimumSize(200, 200);
		color_scale_layout->setMaximumSize(200, 200);
		sub_legend_layout->addElement(0, 0, color_scale_layout);
		custom_plot::patch::add_title(draw_area, sub_legend_layout, gs.get_legend_title(title), gs.get_legend_title_font());
	};

	void proportion_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		const QString& legend_title,
		const GraphSettings& gs
	) {
		QCPLayoutGrid* sub_legend_layout = new QCPLayoutGrid;
		auto [row, col] = custom_plot::patch::find_next_empty_position(legend_layout);
		legend_layout->addElement(row, col, sub_legend_layout);
		QCPAxisRect* legend_rect = new QCPAxisRect(draw_area, true);
		sub_legend_layout->addElement(0, 0, legend_rect);
		sub_legend_layout->setRowSpacing(10);
		custom_plot::patch::set_single_round_legend(draw_area, legend_rect, "0%", Qt::black, gs.get_legend_label_font(), 10, 0, 1);
		custom_plot::patch::set_single_round_legend(draw_area, legend_rect, "25%", Qt::black, gs.get_legend_label_font(), 10, 1, 6);
		custom_plot::patch::set_single_round_legend(draw_area, legend_rect, "50%", Qt::black, gs.get_legend_label_font(), 10, 2, 11);
		custom_plot::patch::set_single_round_legend(draw_area, legend_rect, "75%", Qt::black, gs.get_legend_label_font(), 10, 3, 16);
		custom_plot::patch::set_single_round_legend(draw_area, legend_rect, "100%", Qt::black, gs.get_legend_label_font(), 10, 4, 21);
		custom_plot::patch:: remove_left_bottom_axis(legend_rect);
		int legend_column_width = gs.get_legend_column_width();
		if (legend_column_width == 0) {
			legend_column_width = 40;
		}
		custom_plot::patch::set_fixed_size(legend_rect, legend_column_width, gs.get_legend_row_width() * 5);
		custom_plot::patch::set_range(legend_rect, QCPRange(-10, legend_column_width - 10), QCPRange(-0.5, 4.5));
		custom_plot::patch::add_title(draw_area, sub_legend_layout, legend_title, gs.get_legend_title_font());
	};

	void add_square_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		const QStringList& levels,
		const QVector<QColor>& colors,
		QFont font,
		const QString& legend_title,
		QFont legend_title_font,
		int legend_column_width,
		int legend_row_width
	) {
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
		td.setDefaultFont(font);

		for (int i = 0; i < ncol; ++i) {
			int column_width = 0;
			int row = rows[i];
			QCPAxisRect* legend_rect = new QCPAxisRect(draw_area, true);
			item_layout_->addElement(0, i, legend_rect);
			int j;
			for (j = 0; j < row; ++j, ++index) {
				custom_plot::patch::set_single_square_legend(draw_area, legend_rect, levels[index], colors[index], font, 15, (row - j - 1) * 25 + 10, 10);
				td.setHtml(levels[index]);
				auto [width, _] = td.size().toSize();
				column_width = column_width > width ? column_width : width;
			}
			if (legend_column_width == 0) {
				legend_column_width = column_width + 25;
			}
			custom_plot::patch::remove_left_bottom_axis(legend_rect);
			custom_plot::patch::set_fixed_size(legend_rect, legend_column_width, legend_row_width * j);
			custom_plot::patch::set_range(legend_rect, QCPRange(-10, legend_column_width - 10), QCPRange(0, j * 25));
		}
		custom_plot::patch::add_title(draw_area, sub_legend_layout, legend_title, legend_title_font);
	};

	void add_round_legend(
		QCustomPlot* draw_area,
		QCPLayoutGrid* legend_layout,
		const QStringList& levels,
		const QList<QColor>& colors,
		const QString& legend_title,
		const GraphSettings& gs,
		int legend_index
	) {
		custom_plot::patch::add_round_legend(
			draw_area,
			legend_layout,
			levels,
			colors,
			gs.get_legend_title(legend_title, legend_index),
			gs.get_legend_column_width(),
			gs.get_legend_row_width(),
			gs.get_legend_title_font(),
			gs.get_legend_label_font());
	};

}
