#pragma once

#include "Identifier.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

#include "CustomPlot.h"

class CoveragePlotWorker :
	public QObject
{
	Q_OBJECT
public:

	CoveragePlotWorker(
		const SingleCellMultiome* single_cell_multiome,
		const Fragments* fragments,
		const QStringList& factors,
		const QStringList& levels,
		const QString& region,
		bool draw_gene,
		bool draw_link,
		double link_cutoff,
		bool draw_legend);

	CoveragePlotWorker(
		const SingleCellMultiome* single_cell_multiome,
		const Fragments* fragments,
		const QStringList& factors,
		const QStringList& levels,
		const QString& region,
		const QString& ccan_name,
		const QVector<std::pair<int, int>>& ccan_locs,
		bool draw_gene,
		bool draw_link,
		double link_cutoff,
		bool draw_legend);

	CoveragePlotWorker(
		const SingleCellAtac* single_cell_atac,
		const Fragments* fragments,
		const QStringList& factors,
		const QStringList& levels,
		const QString& region,
		bool draw_gene,
		bool draw_link,
		double link_cutoff,
		bool draw_legend);

	CoveragePlotWorker(
		const SingleCellAtac* single_cell_atac,
		const Fragments* fragments,
		const QStringList& factors,
		const QStringList& levels,
		const QString& region,
		const QString& ccan_name,
		const QVector<std::pair<int, int>>& ccan_locs,
		bool draw_gene,
		bool draw_link,
		double link_cutoff,
		bool draw_legend);

	enum class WorkMode {
		MultiomeFragmentsObject,
		MultiomeFragmentsObjectCcan,
		AtacFragmentsObject,
		AtacFragmentsObjectCcan
	};

	WorkMode mode_{ WorkMode::MultiomeFragmentsObject };

	const SingleCellMultiome* single_cell_multiome_{ nullptr };
	const SingleCellAtac* single_cell_atac_{ nullptr };

	const Fragments* fragments_{ nullptr };

	QStringList factors_;
	QStringList levels_;
	QString region_;
	QString ccan_name_;

	QVector<std::pair<int, int>> ccan_locs_;

	bool draw_gene_{ false };
	bool draw_link_{ false };
	bool draw_legend_{ false };

	double link_cutoff_{ 0.0 };

	GenomicRange genome_;

	Location location_;

	bool get_location_by_gene_name_ = false;

	std::vector<int> cell_index_;

	std::vector<double> cell_size_;

	QMap<QString, int> group_distribution_;

	QMap<QString, QPair<char, QList<std::tuple<int, int, QString> > > > gene_structure_;

	QList< std::pair<int, int> > peak_locations_;
	QVector<int> peak_index_;

	QList<std::tuple<int, int, double>> peak_links_;

	Eigen::MatrixXd insertion_matrix_;

	Eigen::MatrixXd normalized_matrix_;

	QVector<double> draw_bottom_axis_;

	bool load_genome();

	bool find_gene_in_region();
	bool find_peak_in_region();
	bool find_link_in_region();

	void expand_plot_region();

	void build_index();

	bool calculate_fragments_size();

	bool calculate_matrix();

	void smooth_matrix();

	bool get_location();

	bool get_location_by_name();

	bool prepare_draw_matrix();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_plot_ready(COVERAGE_PLOT_ELEMENTS);

};

