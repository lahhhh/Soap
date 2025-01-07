#pragma once

#include "Identifier.h"

#include "SingleCellRna.h"

/*
* modified from scvelo package
* Bergen, Volker et al. ¡°Generalizing RNA velocity to transient cell states through dynamical modeling.¡±
* Nature biotechnology vol. 38,12 (2020): 1408-1414. doi:10.1038/s41587-020-0591-3
*/

class ScveloWorker
	: public QObject
{
	Q_OBJECT
public:

	ScveloWorker(
		const VelocytoBase* velocyto_base,
		const Eigen::MatrixXd& mat
	) : 
		velocyto_base_(velocyto_base),
		mat_(mat)
	{}

	Eigen::MatrixXd mat_;
	QString metric_ = "Euclidean";

	const VelocytoBase* velocyto_base_ = nullptr;

	int n_neighbors_ = 30;

	Eigen::MatrixXd knn_distance_;
	Eigen::MatrixXi knn_indices_;

	std::vector<std::vector<int>> connectivities_;

	Eigen::SparseMatrix<double> spliced_normalized_;
	Eigen::SparseMatrix<double> unspliced_normalized_;

	Eigen::MatrixXd ms_;
	Eigen::MatrixXd mu_;
	
	Eigen::MatrixXd residual_; // velocity
	//Eigen::MatrixXd residual2_; // variance velocity

	Eigen::SparseMatrix<double> strengths_;

	std::unique_ptr<ScveloEstimate> res_{ nullptr };

	bool normalize_counts();

	bool moments();

	std::pair<Eigen::ArrayXd, Eigen::ArrayXd> smooth_knn_distance(double local_connectivity = 1.0, int n_iter = 64);

	void compute_membership_strengths(const Eigen::ArrayXd& sigmas, const Eigen::ArrayXd& rho);

	void fuzzy_simplicial_set();

	void get_connectivities();

	bool get_velocity();

	bool compute_deterministic();

	void compute_stochastic();

	void velocity_graph();

	static Eigen::ArrayXd cosine_correlation(const Eigen::MatrixXd& mat, const Eigen::ArrayXd& arr);

public:

	bool work();

public slots:

	void run();	

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_scvelo_ready(ScveloEstimate*);	
};

