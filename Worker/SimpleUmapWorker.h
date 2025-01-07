#pragma once

#include "Identifier.h"

class SimpleUmapWorker :
	public QObject
{
	Q_OBJECT
public:
	SimpleUmapWorker(
		const Eigen::MatrixXd& mat) :
		mat_(mat)
	{};

	Eigen::MatrixXd mat_;

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
