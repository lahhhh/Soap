#pragma once

#include "Identifier.h"

#include "CoverageTrack.h"

struct COVERAGE_PLOT_ELEMENTS {

	QString region_name;
	QString plot_title;
	Location region;
	QStringList group_factors;
	Eigen::MatrixXd normalized_matrix;
	QMap<QString, QPair<char, QList<std::tuple<int, int, QString>>>> gene_structure;
	QList< std::pair<int, int>> peak_locations;
	QList< std::tuple<int, int, double>> peak_links;
	QList< std::pair<int, int>> ccan_locations;
	bool draw_legend{ false };
	bool draw_bac{ false };
	bool draw_ac{ true };

	void clear() {
		region_name.clear();
		plot_title.clear();
		group_factors.clear();
		normalized_matrix.resize(0, 0);
		gene_structure.clear();
		peak_locations.clear();
		peak_links.clear();
		ccan_locations.clear();
	}

};

struct STREAM_PLOT_ELEMENTS {
	Eigen::MatrixXd embedding;
	Eigen::ArrayXd x;
	Eigen::ArrayXd y;
	Eigen::MatrixXd u;
	Eigen::MatrixXd v;
	Eigen::MatrixX<bool> mask;
	QStringList embedding_names;
	QStringList graph_settings;
};

struct VELO_GRID_PLOT_ELEMENTS {
	Eigen::MatrixXd embedding;
	Eigen::MatrixXd arrows_start;
	Eigen::MatrixXd direction;
	QStringList embedding_names;
	QStringList graph_settings;
};

struct HEATMAP_PLOT_ELEMENTS {

	Eigen::MatrixXd mat;
	QString legend_title;
	QStringList column_names;
	QStringList row_names;

	QList<QPair<QString, QStringList>> column_annotations;
	QList<QPair<QString, QStringList>> row_annotations;

	bool show_annotation_legend{ true };

	bool annotate_column_at_top{ true };
	bool annotate_row_at_left{ true };

	bool show_column_names{ false };
	bool show_row_names{ false };
	bool show_column_names_at_bottom{ true };
	bool show_row_names_at_left{ true };
	bool show_column_annotation_name{ false };
	bool show_row_annotation_name{ false };
	bool show_column_annotation_name_at_left{ true };
	bool show_row_annotation_name_at_bottom{ true };

	QMap<QString, QString> info;
};