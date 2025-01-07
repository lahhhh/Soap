#pragma once

#include "Identifier.h"

#include "SingleCellMultiome.h"

struct ATAC_LANDSCAPE_PLOT_ELEMENTS {

	Eigen::MatrixXd mat;

	QList<std::tuple<QString, int, int>> cell_loc;

	QList<std::tuple<QString, int, int>> chr_loc;

	bool scale{ false };
};

class AtacLandscapePlotWorker
	: public QObject
{
	Q_OBJECT
public:

	AtacLandscapePlotWorker(
		const Fragments* fragments,
		const QString& normalize_method,
		const QStringList& levels,
		const QStringList& factors,
		const QStringList& chromosomes,
		bool scale,
		bool log_transform
	):
		fragments_(fragments),
		normalize_method_(normalize_method),
		levels_(levels),
		factors_(factors),
		chromosomes_(chromosomes),
		scale_(scale),
		log_transform_(log_transform)
	{};

	const Fragments* fragments_{ nullptr };

	QString normalize_method_;

	QStringList levels_;

	QStringList factors_;

	QStringList chromosomes_;

	QList<std::tuple<QString, int, int>> cell_loc_;
	QVector<int> cell_index_;

	QList<std::tuple<QString, int, int>> chr_loc_;
	QMap<QString, int> chr_index_;
	QMap<QString, int> chr_sizes_;

	int nrow_{ 0 };
	int ncol_{ 0 };
	Eigen::MatrixXd plot_mat_;

	int resolution_{ 100000 };

	bool scale_{ false };
	bool log_transform_{ false };

	bool create_index1();

	bool create_index2();

	bool calculate_matrix();


public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_plot_ready(ATAC_LANDSCAPE_PLOT_ELEMENTS);
};

