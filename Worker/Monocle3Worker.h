#pragma once

#include "Identifier.h"

#include "Embedding.h"
#include "Monocle3.h"

#include <igraph.h>

/*
* modified from Monocle3 software package
* Qiu, Xiaojie et al.¡°Reversed graph embedding resolves complex single - cell trajectories.¡± 
* Nature methods vol. 14, 10 (2017) : 979 - 982. doi : 10.1038 / nmeth.4402
*/

class Monocle3Worker
	: public QObject
{
	Q_OBJECT

public:

	Monocle3Worker(
		const Embedding& embedding,
		const Eigen::ArrayX<bool>& cell_included
	)
		:
		original_embedding_(embedding),
		cell_choosed_(cell_included)
	{};

	Embedding original_embedding_;

	// [ndim, n_cell]
	Eigen::MatrixXd embedding_;

	QStringList clusters_;

	Eigen::ArrayX<bool> cell_choosed_;

	bool close_loop_{ true };

	int max_iter_{ 10 };

	double l1_gamma_{ 0.5 };
	double l1_sigma_{ 0.01 };
	double epsilon_{ 1e-5 };

	igraph_t pr_graph_;
	igraph_t cell_graph_;

	std::unique_ptr<Monocle3> res_{ nullptr };

	void learn_graph();

	std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, QVector<double> > calculate_principal_graph(
		const Eigen::MatrixXd& data,
		const Eigen::MatrixXd& centroids,
		int max_iter = 10,
		double epsilon = 1e-5,
		double l1_gamma = 0.5,
		double l1_sigma = 0.05);

	void project_to_mst(
		bool orthogonal_proj_tip,
		const Eigen::MatrixXd& pr_node_embedding);

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_monocle3_ready(Monocle3*);
};

