#include "LeidenPartitionWorker.h"

#include <leiden/CPMVertexPartition.h>
#include <leiden/RBConfigurationVertexPartition.h>
#include <leiden/RBERVertexPartition.h>
#include <leiden/SignificanceVertexPartition.h>
#include <leiden/SurpriseVertexPartition.h>
#include <leiden/ModularityVertexPartition.h>

#include "Custom.h"

QVector<int> leiden_cluster(
	const Eigen::MatrixXd& mat,
	const QString& method,
	const QString& nn_method,
	int n_neighbors,
	int n_trees,
	double resolution) {

	Eigen::MatrixXi knn;

	if (nn_method == "Euclidean") {
		knn = custom::get_knn_mt<Euclidean>(mat, n_neighbors, n_trees);
	}
	else if (nn_method == "Angular") {
		knn = custom::get_knn_mt<Angular>(mat, n_neighbors, n_trees);
	}
	else /* if (nn_method == "Manhattan") */ {
		knn = custom::get_knn_mt<Manhattan>(mat, n_neighbors, n_trees);
	}
	Eigen::SparseMatrix<double> snn = custom::create_shared_nearest_neighbors_matrix(knn);

	Optimiser otm;

	std::unique_ptr<MutableVertexPartition> partition;

	auto graph = create_graph(snn);
	if (method == "CPM") {
		partition.reset(otm.find_partition<CPMVertexPartition>(graph, resolution));
	}
	else if (method == "RBConfiguration") {
		partition.reset(otm.find_partition<RBConfigurationVertexPartition>(graph, resolution));
	}
	else if (method == "RBER") {
		partition.reset(otm.find_partition<RBERVertexPartition>(graph, resolution));
	}
	else if (method == "Significance") {
		partition.reset(otm.find_partition<SignificanceVertexPartition>(graph));
	}
	else if (method == "Surprise") {
		partition.reset(otm.find_partition<SurpriseVertexPartition>(graph));
	}
	else /* if (method == "Modularity") */ {
		partition.reset(otm.find_partition<ModularityVertexPartition>(graph));
	}
	QVector<int> clusters = custom::sapply(partition->membership(), [](std::size_t val) { return static_cast<int>(val); });
	return clusters;
};

Graph* create_graph(const Eigen::SparseMatrix<double>& snn) {
	igraph_sparsemat_t trans_mat, graph_mat;
	const int nrow = snn.rows(), ncol = snn.cols(); // nrow == ncol
	const Eigen::Index size = snn.nonZeros();

	igraph_sparsemat_init(&trans_mat, nrow, ncol, size);

	for (Eigen::Index i = 0; i < ncol; ++i) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
			igraph_sparsemat_entry(&trans_mat, it.row(), i, it.value());
		}
	}

	igraph_sparsemat_compress(&trans_mat, &graph_mat);
	igraph_sparsemat_destroy(&trans_mat);

	igraph_t* graph = new igraph_t();
	igraph_vector_t weights;
	igraph_vector_init(&weights, 0);
	igraph_sparse_weighted_adjacency(graph, &graph_mat, IGRAPH_ADJ_UNDIRECTED, &weights, IGRAPH_NO_LOOPS);
	igraph_sparsemat_destroy(&graph_mat);

	const int n_edge = igraph_ecount(graph);
	std::vector<double> edge_weights(n_edge);
	for (int i = 0; i < n_edge; ++i) {
		edge_weights[i] = VECTOR(weights)[i];
	}
	igraph_vector_destroy(&weights);

	return Graph::GraphFromEdgeWeights(graph, edge_weights);
};

void LeidenPartitionWorker::create_shared_nearest_neighbors_matrix() {
	if (this->nn_method_ == "Euclidean") {
		this->knn_ = custom::get_knn_mt<Euclidean>(this->mat_, this->n_neighbors_, this->n_trees_);
	}
	else if (this->nn_method_ == "Angular") {
		this->knn_ = custom::get_knn_mt<Angular>(this->mat_, this->n_neighbors_, this->n_trees_);
	}
	else /* if (this->nn_method_ == "Manhattan") */ {
		this->knn_ = custom::get_knn_mt<Manhattan>(this->mat_, this->n_neighbors_, this->n_trees_);
	}

	this->snn_ = custom::create_shared_nearest_neighbors_matrix(this->knn_);
};

LeidenPartitionWorker::LeidenPartitionWorker(
	const Eigen::MatrixXd& mat,
	const QString& method,
	const QString& nn_method,
	int n_neighbors,
	int n_trees,
	double resolution
) :
	mat_(mat),
	method_(method),
	nn_method_(nn_method),
	n_neighbors_(n_neighbors),
	n_trees_(n_trees),
	resolution_(resolution)
{};

void LeidenPartitionWorker::run() {

	this->create_shared_nearest_neighbors_matrix();

	this->find_partition();

	G_TASK_END;
};

void LeidenPartitionWorker::find_partition() {

	Optimiser otm;

	std::unique_ptr<MutableVertexPartition> partition;

	auto graph = create_graph(this->snn_);
	
	if (this->method_ == "CPM") {
		partition.reset(otm.find_partition<CPMVertexPartition>(graph, this->resolution_));
	}
	else if (this->method_ == "RBConfiguration") {
		partition.reset(otm.find_partition<RBConfigurationVertexPartition>(graph, this->resolution_));
	}
	else if (this->method_ == "RBER") {
		partition.reset(otm.find_partition<RBERVertexPartition>(graph, this->resolution_));
	}
	else if (this->method_ == "Significance") {
		partition.reset(otm.find_partition<SignificanceVertexPartition>(graph));
		G_TASK_NOTICE("Parameter : resolution is not used in Significance Vertex Partition");
	}
	else if (this->method_ == "Surprise") {
		partition.reset(otm.find_partition<SurpriseVertexPartition>(graph));
		G_TASK_NOTICE("Parameter : resolution is not used in Surprise Vertex Partition");
	}
	else /* if (this->method_ == "Modularity") */ {
		partition.reset(otm.find_partition<ModularityVertexPartition>(graph));
		G_TASK_NOTICE("Parameter : resolution is not used in Modularity Vertex Partition");
	}

	igraph_t* g = graph->get_igraph();
	igraph_destroy(g);
	delete g;

	QVector<int> classification = custom::sapply(partition->membership(), [](std::size_t val) { return static_cast<int>(val); });

	G_TASK_LOG(QString::number(custom::unique_element_number(classification)) + " clusters were found in Leiden Partition.");

	emit x_leiden_ready(classification);

};