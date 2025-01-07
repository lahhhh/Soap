#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

class VelocytoDownstreamWorker
	: public QObject
{
	Q_OBJECT
public:

	VelocytoDownstreamWorker(const VelocytoBase* velocyto_base): velocyto_base_(velocyto_base){}

	const VelocytoBase* velocyto_base_ = nullptr;

	int n_cell_ = 10;
	int n_gene_ = 1;
	int scale_factor_ = 1000;
	
	double pseudo_count_ = 1.0;

	Eigen::MatrixXd ko_;

	Eigen::MatrixXd conv_nmat_norm_1;
	Eigen::MatrixXd conv_emat_norm_1;

	Eigen::MatrixXd conv_nmat_norm_2;
	Eigen::MatrixXd conv_emat_norm_2;

	QVector<int> ko_index_;
	QVector<int> gamma_index_;

	Eigen::ArrayXd emat_cs_;

	std::unique_ptr<VelocityEstimate> res_{nullptr};

	bool gene_relative_velocity_estimates();

	void get_conv_mat();

	// use [0, 5] | [95, 100]
	bool fit_linear();

	void calculate_velocity_shift();

	void calculate_extrapolated_cell_state();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_estimate_ready(VelocityEstimate*);
	
};

