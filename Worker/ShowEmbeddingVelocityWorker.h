#pragma once

#include "Identifier.h"

#include "VelocityEstimate.h"
#include "Embedding.h"

#include "CustomPlot.h"

class ShowEmbeddingVelocityWorker : 
	public QObject
{
	Q_OBJECT
public:

	ShowEmbeddingVelocityWorker(
		const VelocityEstimate* estimate,
		const Embedding& embedding,
		const QStringList& graph_settings,
		int graph_mode
	):
		estimate_(estimate),
		embedding_(embedding.data_.mat_(Eigen::all, std::vector<int>{0,1})),
		embedding_names_(embedding.data_.colnames_),
		graph_settings_(graph_settings),
		graph_mode_(graph_mode)
	{}

	const VelocityEstimate* estimate_ = nullptr;

	Eigen::MatrixXd embedding_;

	QStringList embedding_names_;

	QStringList graph_settings_;

	// 0 : grid; 1 : stream
	int graph_mode_{ 0 };

	int n_cell_neighbor_{ 100 };

	void show_velocity();

	Eigen::MatrixXd col_delta_cor_log10(
		const Eigen::MatrixXd& e,
		const Eigen::MatrixXd& d,
		const double pseudo_count = 1.0
	);

	Eigen::MatrixXd emb_arrows(
		const Eigen::MatrixXd& emb,
		const Eigen::SparseMatrix<double>& tp,
		double arrow_scale = 1.0
	);

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_grid_graph_ready(VELO_GRID_PLOT_ELEMENTS);

	void x_stream_graph_ready(STREAM_PLOT_ELEMENTS);

	
};

