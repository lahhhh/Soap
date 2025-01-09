#include "UmapWorker.h"
#include "annoylib.h"
#include "kissrandom.h"
#include "mman.h"
#include <QDebug>
#include "Umap.h"

UmapWorker::UmapWorker(
	const Eigen::MatrixXd& mat,
	int n_neighbors,
	const QString& metric,
	double learning_rate,
	const QString& init,
	double minimum_distance,
	double spread,
	double set_op_mix_ratio,
	double repulsion_strength,
	int negative_sample_rate,
	int random_state,
	int n_trees
) :
	mat_(mat),
	n_neighbors_(n_neighbors),
	metric_(metric),
	learning_rate_(learning_rate),
	init_(init),
	minimum_distance_(minimum_distance),
	spread_(spread),
	set_op_mix_ratio_(set_op_mix_ratio),
	repulsion_strength_(repulsion_strength),
	negative_sample_rate_(negative_sample_rate),
	random_state_(random_state),
	n_trees_(n_trees)
{}

bool UmapWorker::work() {

	UMAP umap(
		&this->mat_,
		this->n_neighbors_,
		this->metric_,
		this->learning_rate_,
		this->init_,
		this->minimum_distance_,
		this->spread_,
		this->set_op_mix_ratio_,
		this->repulsion_strength_,
		this->negative_sample_rate_,
		this->random_state_,
		this->n_trees_
	);

	umap.fit();

	this->res_ = umap.embedding_;

	return true;
};

void UmapWorker::run() {

	G_TASK_LOG("UMAP start...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_umap_ready(this->res_);

	G_TASK_LOG("UMAP finished.");

	G_TASK_END;
}
