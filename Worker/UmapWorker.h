#pragma once

#include "Identifier.h"

class UmapWorker : public QObject
{
	Q_OBJECT
public:
	UmapWorker(
		const Eigen::MatrixXd& mat,
		int n_neighbors = 30,
		const QString& metric = "Angular",
		double learning_rate = 1.0,
		const QString& init = "spectral",
		double minimum_distance = 0.1,
		double spread = 1.0,
		double set_op_mix_ratio = 1.0,
		double repulsion_strength = 1.0,
		int negative_sample_rate = 5,
		int random_state = 1997,
		int n_trees = 50
	);

	Eigen::MatrixXd mat_;

	QString metric_;
	QString init_;

	double learning_rate_{ 1.0 };
	double minimum_distance_{ 0.1 };
	double spread_{ 1.0 };
	double set_op_mix_ratio_{ 1.0 };
	double repulsion_strength_{ 1.0 };

	int n_epochs_{ 0 };
	int n_neighbors_{ 30 };
	int negative_sample_rate_{ 5 };
	int random_state_{ 1997 };
	int n_trees_{ 50 };

	Eigen::MatrixXd res_;

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_umap_ready(Eigen::MatrixXd);
};

