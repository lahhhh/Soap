#include "SimpleUmapWorker.h"

#include "TruncatedSvd.h"
#include "Umap.h"

void SimpleUmapWorker::run() {

	auto [U, S, V] = tsvd(&this->mat_, 50, 1997, 0.00001, 1000);

	Eigen::MatrixXd emb = U * S.asDiagonal();

	Eigen::MatrixXd umap_input = emb.block(0, 0, emb.rows(), 20);

	UMAP umap(
		&umap_input
	);

	umap.fit();

	emit x_umap_ready(umap.embedding_);

	G_TASK_END;
};