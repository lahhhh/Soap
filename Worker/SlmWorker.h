#pragma once

#include "Identifier.h"

#include <QString>

#include "SLM/VOSClusteringTechnique.h"

class SlmWorker:
	public QObject
{
	Q_OBJECT

public:
	SlmWorker(
		const Eigen::MatrixXd& mat,
		const QString& method,
		const QString& nn_method,
		int modularity_function,
		int n_random_start,
		int n_iterations,
		int random_seed,
		int n_neighbors,
		int n_trees,
		double resolution
	);

	Eigen::MatrixXd mat_;

	QString method_;
	QString nn_method_;

	int modularity_function_ = 0;
	int n_random_start_ = 0;
	int n_iterations_ = 0;
	int random_seed_ = 0;
	int n_neighbors_ = 0;
	int n_trees_ = 0;
	int n_nodes_ = 0;

	double resolution_ = 0.0;

	Eigen::MatrixXi knn_;

	Eigen::SparseMatrix<double> snn_;

	Network network_;

	Clustering clustering_;

	std::default_random_engine dre_;
	
public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_cluster_ready(std::vector<int>);
	
};

