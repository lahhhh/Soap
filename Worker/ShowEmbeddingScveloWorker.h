#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"
#include "CustomPlot.h"

class ShowEmbeddingScveloWorker : public QObject
{
	Q_OBJECT
public:
	
	ShowEmbeddingScveloWorker(
		const ScveloEstimate* estimate,
		const Embedding& embedding,
		const QStringList& graph_settings,
		int graph_mode
	) :
		estimate_(estimate),
		embedding_(embedding.data_.mat_(Eigen::all, std::vector<int>{ 0,1 })),
		embedding_names_(embedding.data_.colnames_),
		graph_settings_(graph_settings),
		graph_mode_(graph_mode)
	{}

	const ScveloEstimate* estimate_ = nullptr;

	Eigen::MatrixXd embedding_;
	Eigen::MatrixXd velocity_embedding_;

	QStringList embedding_names_;

	QStringList graph_settings_;

	// 0 : grid; 1 : stream
	int graph_mode_{ 0 };

	void velocity_embedding();

	void compute_velocity_on_grid();

	Eigen::SparseMatrix<double> transition_matrix();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_grid_graph_ready(VELO_GRID_PLOT_ELEMENTS);

	void x_stream_graph_ready(STREAM_PLOT_ELEMENTS);
		
};

