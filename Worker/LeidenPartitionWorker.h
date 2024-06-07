#pragma once

#include "Identifier.h"

#include <QString>

#include <igraph.h>
#include <leiden/Optimiser.h>

static Graph* create_graph(const Eigen::SparseMatrix<double>& snn);

QVector<int> leiden_cluster(
	const Eigen::MatrixXd& mat,
	const QString& method,
	const QString& nn_method,
	int n_neighbors,
	int n_trees,
	double resolution);

class LeidenPartitionWorker:
	public QObject
{
	Q_OBJECT
public:

	LeidenPartitionWorker(
		const Eigen::MatrixXd& mat,
		const QString& method,
		const QString& nn_method,
		int n_neighbors,
		int n_trees,
		double resolution
	);

	Eigen::MatrixXd mat_;

	QString method_;
	QString nn_method_;

	int n_neighbors_ = 0;
	int n_trees_ = 0;

	double resolution_ = 0.0;

	Eigen::MatrixXi knn_;

	Eigen::SparseMatrix<double> snn_;

	// modified from seurat
	void create_shared_nearest_neighbors_matrix();

	void find_partition();


public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_leiden_ready(QVector<int>);
	
};

