#include "SimpleUmapWorker.h"

#include "TruncatedSvd.h"
#include "Umap.h"

bool SimpleUmapWorker::work() {

	auto [U, S, V] = tsvd(&this->mat_, 50, 1997, 0.00001, 1000);

	Eigen::MatrixXd emb = U * S.asDiagonal();

	Eigen::MatrixXd umap_input = emb.block(0, 0, emb.rows(), 20);

	UMAP umap(
		&umap_input
	);

	umap.fit();

	this->res_ = umap.embedding_;

	return true;
};

void SimpleUmapWorker::run() {

	G_TASK_LOG("Start UMAP...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_umap_ready(this->res_);

	G_TASK_LOG("UMAP finished.");

	G_TASK_END;
};