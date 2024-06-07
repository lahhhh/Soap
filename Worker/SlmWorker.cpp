#include "SlmWorker.h"

#include "Custom.h"

#include "Louvain.h"

SlmWorker::SlmWorker(
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
) :
    mat_(mat),
    method_(method),
    nn_method_(nn_method),
	modularity_function_(modularity_function),
	n_random_start_(n_random_start),
	n_iterations_(n_iterations),
	random_seed_(random_seed),
    n_neighbors_(n_neighbors),
    n_trees_(n_trees),
    n_nodes_(mat.rows()),
    resolution_(resolution)
{
    this->dre_.seed(random_seed);
}

void SlmWorker::run() {

    auto cluster = louvain_cluster(
        this->mat_,
        this->method_,
        this->nn_method_,
        this->modularity_function_,
        this->n_random_start_,
        this->n_iterations_,
        this->random_seed_,
        this->n_neighbors_,
        this->n_trees_,
        this->resolution_
    );

    emit x_cluster_ready(cluster);

    G_TASK_LOG(QString::number(_Cs unique_element_number(cluster)) + " clusters were found in partition.");

    G_TASK_END;
}