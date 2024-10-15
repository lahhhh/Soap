#include "Monocle3Worker.h"

#include <igraph.h>
#include "Custom.h"
#include "pval.h"

// from rigraph
int igraph_i_neighbors(const igraph_t* graph, igraph_vector_t* neis, igraph_integer_t pnode,
	igraph_neimode_t mode, igraph_loops_t loops, igraph_multiple_t multiple) {
#define DEDUPLICATE_IF_NEEDED(vertex, n)                                                 \
    if (should_filter_duplicates) {                                                        \
        if (                                                                               \
            (loops == IGRAPH_NO_LOOPS && vertex == pnode) ||                               \
            (loops == IGRAPH_LOOPS_ONCE && vertex == pnode && last_added == pnode)         \
        ) {                                                                                \
            length -= n;                                                                   \
            if (loops == IGRAPH_LOOPS_ONCE) {                                              \
                last_added = -1;                                                           \
            }                                                                              \
            continue;                                                                      \
        } else if (multiple == IGRAPH_NO_MULTIPLE && vertex == last_added) {               \
            length -= n;                                                                   \
            continue;                                                                      \
        } else {                                                                           \
            last_added = vertex;                                                           \
        }                                                                                  \
    }

	long int length = 0, idx = 0;
	long int i, j;

	long int node = pnode;
	igraph_real_t last_added = -1;
	igraph_bool_t should_filter_duplicates;

	//if (node < 0 || node > igraph_vcount(graph) - 1) {
	//	IGRAPH_ERROR("Given vertex is not in the graph.", IGRAPH_EINVVID);
	//}
	//if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
	//	mode != IGRAPH_ALL) {
	//	IGRAPH_ERROR("Mode should be either IGRAPH_OUT, IGRAPH_IN or IGRAPH_ALL.", IGRAPH_EINVMODE);
	//}

	if (!igraph_is_directed(graph)) {
		mode = IGRAPH_ALL;
	}

	//if (mode != IGRAPH_ALL && loops == IGRAPH_LOOPS_TWICE) {
	//	IGRAPH_ERROR("For a directed graph (with directions not ignored), "
	//		"IGRAPH_LOOPS_TWICE does not make sense.\n", IGRAPH_EINVAL);
	//}
	/* Calculate needed space first & allocate it */
	/* Note that 'mode' is treated as a bit field here; it's okay because
	 * IGRAPH_ALL = IGRAPH_IN | IGRAPH_OUT, bit-wise */
	if (mode & IGRAPH_OUT) {
		length += (VECTOR(graph->os)[node + 1] - VECTOR(graph->os)[node]);
	}
	if (mode & IGRAPH_IN) {
		length += (VECTOR(graph->is)[node + 1] - VECTOR(graph->is)[node]);
	}

	igraph_vector_resize(neis, length);

	/* The loops below produce an ordering what is consistent with the
	 * ordering returned by igraph_neighbors(), and this should be preserved.
	 * We are dealing with two sorted lists; one for the successors and one
	 * for the predecessors. If we have requested only one of them, we have
	 * an easy job. If we have requested both, we need to merge the two lists
	 * to ensure that the output is sorted by the vertex IDs of the "other"
	 * endpoint of the affected edges. We don't need to merge if the graph
	 * is undirected, because in that case the data structure guarantees that
	 * the "out-edges" contain only (u, v) pairs where u <= v and the
	 * "in-edges" contains the rest, so the result is sorted even without
	 * merging. */
	if (!igraph_is_directed(graph) || mode != IGRAPH_ALL) {
		/* graph is undirected or we did not ask for both directions in a
		 * directed graph; this is the easy case */

		should_filter_duplicates = !(multiple == IGRAPH_MULTIPLE &&
			((!igraph_is_directed(graph) && loops == IGRAPH_LOOPS_TWICE) ||
				(igraph_is_directed(graph) && loops != IGRAPH_NO_LOOPS)));

		if (mode & IGRAPH_OUT) {
			j = (long int)VECTOR(graph->os)[node + 1];
			for (i = (long int)VECTOR(graph->os)[node]; i < j; i++) {
				igraph_real_t to = VECTOR(graph->to)[(long int)VECTOR(graph->oi)[i]];
				DEDUPLICATE_IF_NEEDED(to, 1);
				VECTOR(*neis)[idx++] = to;
			}
		}

		if (mode & IGRAPH_IN) {
			j = (long int)VECTOR(graph->is)[node + 1];
			for (i = (long int)VECTOR(graph->is)[node]; i < j; i++) {
				igraph_real_t from = VECTOR(graph->from)[(long int)VECTOR(graph->ii)[i]];
				DEDUPLICATE_IF_NEEDED(from, 1);
				VECTOR(*neis)[idx++] = from;
			}
		}
	}
	else {
		/* Both in- and out- neighbors in a directed graph,
		   we need to merge the two 'vectors' so the result is
		   correctly ordered. */
		long int j1 = (long int)VECTOR(graph->os)[node + 1];
		long int j2 = (long int)VECTOR(graph->is)[node + 1];
		long int i1 = (long int)VECTOR(graph->os)[node];
		long int i2 = (long int)VECTOR(graph->is)[node];
		long int eid1, eid2;
		long int n1, n2;

		should_filter_duplicates = !(multiple == IGRAPH_MULTIPLE &&
			loops == IGRAPH_LOOPS_TWICE);

		while (i1 < j1 && i2 < j2) {
			eid1 = (long int)VECTOR(graph->oi)[i1];
			eid2 = (long int)VECTOR(graph->ii)[i2];
			n1 = (long int)VECTOR(graph->to)[eid1];
			n2 = (long int)VECTOR(graph->from)[eid2];
			if (n1 < n2) {
				i1++;
				DEDUPLICATE_IF_NEEDED(n1, 1);
				VECTOR(*neis)[idx++] = n1;
			}
			else if (n1 > n2) {
				i2++;
				DEDUPLICATE_IF_NEEDED(n2, 1);
				VECTOR(*neis)[idx++] = n2;
			}
			else {
				i1++;
				i2++;
				DEDUPLICATE_IF_NEEDED(n1, 2);
				VECTOR(*neis)[idx++] = n1;
				if (should_filter_duplicates && ((loops == IGRAPH_LOOPS_ONCE && n1 == pnode && last_added == pnode) ||
					(multiple == IGRAPH_NO_MULTIPLE))) {
					length--;
					if (loops == IGRAPH_LOOPS_ONCE) {
						last_added = -1;
					}
					continue;
				}
				VECTOR(*neis)[idx++] = n2;
			}
		}

		while (i1 < j1) {
			eid1 = (long int)VECTOR(graph->oi)[i1++];
			igraph_real_t to = (long int)VECTOR(graph->to)[eid1];
			DEDUPLICATE_IF_NEEDED(to, 1);
			VECTOR(*neis)[idx++] = to;
		}

		while (i2 < j2) {
			eid2 = (long int)VECTOR(graph->ii)[i2++];
			igraph_real_t from = (long int)VECTOR(graph->from)[eid2];
			DEDUPLICATE_IF_NEEDED(from, 1);
			VECTOR(*neis)[idx++] = from;
		}

	}
	//IGRAPH_CHECK(igraph_vector_resize(neis, length));

	return IGRAPH_SUCCESS;
#undef DEDUPLICATE_IF_NEEDED
}

std::pair<Eigen::MatrixXd, double> soft_assignment(const Eigen::MatrixXd& X, const Eigen::MatrixXd& C, const double sigma) {

	int D = X.rows();
	int N = X.cols();
	int K = C.cols();

	Eigen::ArrayXd sq = X.colwise().squaredNorm();

	Eigen::MatrixXd norm_X_sq(N, K);
	norm_X_sq.colwise() = sq.matrix();

	sq = C.colwise().squaredNorm();

	Eigen::MatrixXd norm_C_sq(N, K);
	norm_C_sq.rowwise() = sq.matrix().transpose();

	Eigen::MatrixXd dist_XC = norm_X_sq + norm_C_sq - 2 * X.transpose() * C;

	Eigen::ArrayXd min_dist = dist_XC.rowwise().minCoeff();

	dist_XC.colwise() -= min_dist.matrix();
	Eigen::MatrixXd phi_XC = exp(-dist_XC.array() / sigma);

	Eigen::ArrayXd phi_rs = phi_XC.rowwise().sum();

	Eigen::MatrixXd P = phi_XC.array().colwise() / phi_rs;

	double obj = -sigma * (log(phi_rs) - min_dist / sigma).sum();

	return std::make_pair(P, obj);
}

Eigen::MatrixXd generate_centers(const Eigen::MatrixXd& X, const Eigen::MatrixXd& W, const Eigen::MatrixXd& P, const double gamma) {

	int D = X.rows();

	Eigen::MatrixXd wdiag = W.colwise().sum().asDiagonal();
	Eigen::MatrixXd pdiag = P.colwise().sum().asDiagonal();

	Eigen::MatrixXd Q = 2.0 * (wdiag - W) + gamma * pdiag;
	Eigen::MatrixXd B = gamma * X * P;

	return B * Q.inverse();
}

Eigen::MatrixXd mst(const Eigen::MatrixXd& phi) {

	int n_vertice = phi.rows();

	igraph_matrix_t mat, adjacency;
	igraph_t g, mst;
	igraph_vector_t weights;
	igraph_integer_t i, j;

	igraph_vector_init(&weights, 0);

	igraph_matrix_init(&mat, n_vertice, n_vertice);
	igraph_matrix_init(&adjacency, n_vertice, n_vertice);
	for (j = 0; j < n_vertice; ++j){
		for (i = 0; i < n_vertice; ++i) {
			MATRIX(mat, i, j) = phi(i, j);
		}
	}

	igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_LOWER, &weights, IGRAPH_NO_LOOPS);

	igraph_minimum_spanning_tree_prim(&g, &mst, &weights);
	igraph_get_adjacency(&mst, &adjacency, IGRAPH_GET_ADJACENCY_LOWER, &weights, IGRAPH_NO_LOOPS);

	Eigen::MatrixXd mst_res(n_vertice, n_vertice);
	for (j = 0; j < n_vertice; ++j) {
		for (i = 0; i < n_vertice; ++i) {
			mst_res(i, j) = MATRIX(adjacency, i, j);
		}
	}

	igraph_matrix_destroy(&mat);
	igraph_matrix_destroy(&adjacency);
	igraph_vector_destroy(&weights);
	igraph_destroy(&g);
	igraph_destroy(&mst);

	return mst_res;
}

Eigen::MatrixXd jaccard_coeff(const Eigen::MatrixXi& idx, bool weight) {
	int nrow = idx.rows(), ncol = idx.cols(), r = 0;
	Eigen::MatrixXd weights(nrow * ncol, 3);

	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			int k = idx(i, j);

			weights(r, 0) = i; // i + 1 in R
			weights(r, 1) = k; // k + 1 in R
			weights(r, 2) = 1;

			if (weight) {

				Eigen::ArrayXi nodei = idx.row(i);
				Eigen::ArrayXi nodej = idx.row(k);

				int u = custom::intersect_length(nodei, nodej);  // count intersection number
				int v = 2 * ncol - u;  // count union number

				if (u > 0) {

					weights(r, 2) = (double)u / (double)v;  // normalize the values
				}
			}

			r++;

		}
	}

	weights.col(2) /= weights.col(2).maxCoeff();

	return weights;
}

// data : length_var * n_var
igraph_t cluster_cells_make_graph(
	const Eigen::MatrixXd& data,
	int k
) {
	if (k > data.rows() - 2) {
		k = data.rows() - 2;
	}

	auto [nn_index, nn_distance] = custom::get_knn_mt<Euclidean, true>(data, k);

	nn_index = nn_index(Eigen::all, Eigen::seq(1, nn_index.cols() - 1)).eval();

	nn_distance = nn_distance(Eigen::all, Eigen::seq(1, nn_distance.cols() - 1)).eval();

	Eigen::MatrixXd links = jaccard_coeff(nn_index, false);

	int n_vertice = nn_index.rows();

	igraph_t g;
	igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
	igraph_add_vertices(&g, n_vertice, NULL);

	int n_edge = links.rows();

	igraph_vector_int_t edges;
	igraph_vector_int_init(&edges, 2 * n_edge);


	for (int i = 0; i < n_edge; ++i) {
		VECTOR(edges)[i * 2] = links(i, 0);
		VECTOR(edges)[i * 2 + 1] = links(i, 1);
	}

	igraph_add_edges(&g, &edges, NULL);

	igraph_vector_int_destroy(&edges);

	return g;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_partitions(igraph_t g, igraph_vector_int_t membership, double q_val_threshold) {

	int n_vertice = igraph_vcount(&g);

	igraph_matrix_t adjacency;
	igraph_matrix_init(&adjacency, n_vertice, n_vertice);

	igraph_get_adjacency(&g, &adjacency, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_NO_LOOPS);
	
	Eigen::MatrixXi adj(n_vertice, n_vertice);
	for (int i = 0; i < n_vertice; ++i) {
		for (int j = 0; j < n_vertice; ++j) {
			adj(i, j) = MATRIX(adjacency, i, j);
		}
	}
	igraph_matrix_destroy(&adjacency);

	int n_louvain_cluster{ 0 };
	for (int i = 0; i < n_vertice; ++i) {
		if (VECTOR(membership)[i] > n_louvain_cluster) {
			n_louvain_cluster = VECTOR(membership)[i];
		}
	}
	++n_louvain_cluster;

	Eigen::MatrixXi membership_matrix = Eigen::MatrixXi::Zero(n_vertice, n_louvain_cluster);
	QVector<int> m;
	for (int i = 0; i < n_vertice; ++i) {
		membership_matrix(i, VECTOR(membership)[i]) = 1;
		m << VECTOR(membership)[i];
	}

	Eigen::MatrixXd num_link = (membership_matrix.transpose().eval() * adj * membership_matrix).cast<double>();
	num_link.diagonal().array() = 0.0;
	adj.resize(0, 0);

	Eigen::ArrayXd edges_per_module = num_link.rowwise().sum();

	double total_edges = edges_per_module.sum();
	Eigen::MatrixXd theta = (edges_per_module / total_edges).matrix() *
		(edges_per_module / total_edges).matrix().transpose();

	Eigen::MatrixXd var_null_num_links = theta.array() * (1.0 - theta.array()) / total_edges;
	Eigen::MatrixXd num_link_ij = num_link.array() / total_edges - theta.array();

	Eigen::MatrixXd cluster_mat(n_louvain_cluster, n_louvain_cluster);
	for (int i = 0; i < n_louvain_cluster; ++i) {
		for (int j = 0; j < n_louvain_cluster; ++j) {

			if (var_null_num_links(i, j) == 0.0) {
				cluster_mat(i, j) = 1.0;
			}
			else {
				cluster_mat(i, j) = p_normal(num_link_ij(i, j), 0.0, sqrt(var_null_num_links(i, j)), false);
			}
		}
	}

	Eigen::ArrayXd p_val(n_louvain_cluster * n_louvain_cluster);
	for (int i = 0; i < n_louvain_cluster; ++i) {
		for (int j = 0; j < n_louvain_cluster; ++j) {
			p_val[i * n_louvain_cluster + j] = cluster_mat(i, j);
		}
	}
	p_val = custom::adjust_p_value(p_val, "fdr");
	for (int i = 0; i < n_louvain_cluster; ++i) {
		for (int j = 0; j < n_louvain_cluster; ++j) {
			cluster_mat(i, j) = p_val[i * n_louvain_cluster + j];
		}
	}

	num_link = num_link_ij / total_edges;
	Eigen::MatrixXd sig_links = num_link;

	for (int i = 0; i < n_louvain_cluster; ++i) {
		for (int j = 0; j < n_louvain_cluster; ++j) {

			if (cluster_mat(i, j) > q_val_threshold) {
				sig_links(i, j) = 0.0;
			}
		}
	}

	sig_links.diagonal().array() = 0.0;

	return std::make_pair(cluster_mat, sig_links);
}

std::pair<igraph_t, igraph_vector_int_t> louvain_clustering(const Eigen::MatrixXd& data) {

	igraph_t g = cluster_cells_make_graph(data.transpose(), 25);

	igraph_real_t resolution = 1.0;
	igraph_vector_int_t membership;
	igraph_matrix_int_t memberships;
	igraph_vector_t modularity;

	igraph_vector_init(&modularity, 0);
	igraph_vector_int_init(&membership, 0);
	igraph_matrix_int_init(&memberships, 0, 0);

	igraph_community_multilevel(&g, 0, resolution, &membership, &memberships, &modularity);

	igraph_matrix_int_destroy(&memberships);
	igraph_vector_destroy(&modularity);

	return std::make_pair(g, membership);
}

igraph_t connect_tips(
	const Eigen::MatrixXd& stree_ori, 
	const Eigen::ArrayXi& kmeans_cluster, 
	const Eigen::MatrixXd& pr_node_embedding,
	const Eigen::MatrixXd& cell_embedding,
	int minimal_branch_len
){

	bool prune_tree{ true };

	int n_vertice = stree_ori.rows();

	igraph_t mst_g_old;
	igraph_matrix_t adjacency_ori;
	igraph_integer_t i, j;

	igraph_matrix_init(&adjacency_ori, n_vertice, n_vertice);
	for (j = 0; j < n_vertice; ++j) {
		for (i = 0; i < n_vertice; ++i) {
			MATRIX(adjacency_ori, i, j) = stree_ori(i, j);
		}
	}

	igraph_adjacency(&mst_g_old, &adjacency_ori, IGRAPH_ADJ_UNDIRECTED, IGRAPH_NO_LOOPS);

	igraph_matrix_destroy(&adjacency_ori);

	igraph_vector_t vertice_id;
	igraph_vector_init(&vertice_id, n_vertice);
	for (int j = 0; j < n_vertice; ++j) {
		VECTOR(vertice_id)[j] = j;
	}

	SETVANV(&mst_g_old, "vid", &vertice_id);

	igraph_vector_destroy(&vertice_id);

	igraph_vector_int_t degree;
	igraph_vector_int_init(&degree, 0);
	igraph_degree(&mst_g_old, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

	std::vector<int> degree2(n_vertice);

	for (int i = 0; i < n_vertice; ++i) {
		degree2[i] = VECTOR(degree)[i];
	}

	igraph_vector_int_destroy(&degree);

	auto tip_pc_points = custom::match(degree2, 1);
	
	auto [louvain_g, louvain_membership] = louvain_clustering(cell_embedding);

	int n_cell = kmeans_cluster.size();
	for (int i = 0; i < n_cell; ++i) {
		VECTOR(louvain_membership)[i] = kmeans_cluster[i];
	}

	constexpr double q_val_threshold = 0.05;

	auto [cluster_mat, sig_links] = compute_partitions(louvain_g, louvain_membership, q_val_threshold);

	igraph_vector_int_destroy(&louvain_membership);

	std::vector<std::pair<int, int>> valid_connections;

	int n_louvain_cluster = cluster_mat.rows();

	for (int i = 0; i < n_louvain_cluster; ++i) {
		for (int j = 0; j < n_louvain_cluster; ++j) {

			if (cluster_mat(i, j) < q_val_threshold) {
				if (tip_pc_points.contains(i) && tip_pc_points.contains(j) && i != j) {
					valid_connections.emplace_back(i, j);
				}
			}
		}
	}

	for (int i = 0; i < n_louvain_cluster; ++i) {
		for (int j = 0; j < n_louvain_cluster; ++j) {

			if (cluster_mat(i, j) < q_val_threshold) {
				cluster_mat(i, j) = 1.0;
			}
			else {
				cluster_mat(i, j) = 0.0;
			}
		}
	}

	if (cluster_mat.sum() < 1.0 || valid_connections.empty()) {
		return mst_g_old;
	}

	igraph_real_t diameter_dis;
	igraph_diameter_dijkstra(&mst_g_old, 0, &diameter_dis, 0, 0, 0, 0, IGRAPH_DIRECTED, true);

	auto distance = custom::euclidean_distance_mt(pr_node_embedding);

	igraph_matrix_t adj_mat;
	igraph_matrix_init(&adj_mat, n_vertice, n_vertice);

	for (int i = 0; i < n_vertice; ++i) {
		for (int j = 0; j < n_vertice; ++j) {
			MATRIX(adj_mat, i, j) = distance(i, j);
		}
	}

	igraph_t tmp_g, mst, mst_g;
	igraph_copy(&mst_g ,&mst_g_old);
	igraph_vector_t weights;
	igraph_vector_init(&weights, 0);

	igraph_weighted_adjacency(&tmp_g, &adj_mat, IGRAPH_ADJ_UNDIRECTED, &weights, IGRAPH_NO_LOOPS);
	igraph_minimum_spanning_tree_prim(&tmp_g, &mst, &weights);

	igraph_matrix_destroy(&adj_mat);

	int n_edge = igraph_vector_size(&weights);
	double max_node_dist = VECTOR(weights)[0];
	for (int i = 1; i < n_edge; ++i) {

		double dist = VECTOR(weights)[i];
		if (dist > max_node_dist) {
			max_node_dist = dist;
		}
	}

	constexpr double euclidean_distance_ratio = 1.0;
	constexpr double geodesic_distance_ratio = 0.33;

	std::vector<std::pair<int, int>> added_edges;
	std::vector<int> vertice_to_keep;


	int n_valid_connection = valid_connections.size();
	for (int i = 0; i < n_valid_connection; ++i) {
		auto [v1, v2] = valid_connections[i];
		igraph_matrix_t dist;
		igraph_matrix_init(&dist, 0, 0);

		igraph_distances(&mst_g, &dist, igraph_vss_1(v1), igraph_vss_1(v2), IGRAPH_ALL);
		double distance = MATRIX(dist, 0, 0);
		double distance2 = (pr_node_embedding.col(v1).array() - pr_node_embedding.col(v2).array()).matrix().norm();

		if (distance >= geodesic_distance_ratio * diameter_dis 
			&& distance2 < euclidean_distance_ratio * max_node_dist
			&& distance > minimal_branch_len)
		{
			igraph_add_edge(&mst_g, v1, v2);
			added_edges.emplace_back(v1, v2);

			igraph_vector_int_t path;
			igraph_vector_int_init(&path, 0);

			igraph_get_shortest_path(&mst_g_old, &path, 0, v1, v2, IGRAPH_ALL);
			int len_path = igraph_vector_int_size(&path);

			for (int j = 0; j < len_path; ++j) {
				vertice_to_keep.push_back(VECTOR(path)[j]);
			}
		}
	}

	vertice_to_keep = custom::unique(vertice_to_keep);

	igraph_vector_int_t neighborhood_size;
	igraph_vector_int_init(&neighborhood_size, 0);

	igraph_neighborhood_size(&mst_g_old, &neighborhood_size, igraph_vss_all(), 1, IGRAPH_ALL, 0);

	int root_cell{ 0 };

	for (int i = 0; i < n_vertice; ++i) {
		if (VECTOR(neighborhood_size)[i] == 2) {
			root_cell = i;
			break;
		}
	}

	igraph_vector_int_t parents, order;
	igraph_vector_int_init(&parents, 0);
	igraph_vector_int_init(&order, 0);
	igraph_dfs(&mst_g_old, root_cell, IGRAPH_ALL, false, &order, 0, &parents, 0, 0, 0, 0);

	std::vector<int> vertice_to_delete;

	if (igraph_vector_int_size(&order) > 0) {
		for (int i = 0; i < n_vertice; ++i) {

			int parent = VECTOR(parents)[i];
			if (parent != -1) {
				if (degree2[parent] > 2) {

					igraph_vector_t parent_neighbors;
					igraph_vector_init(&parent_neighbors, 0);
					igraph_i_neighbors(&mst_g_old, &parent_neighbors, parent, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);

					int parent1 = VECTOR(parent_neighbors)[1];
					int parent2 = VECTOR(parent_neighbors)[2];

					igraph_t tmp;
					igraph_copy(&tmp, &mst_g_old);

					int n_neighbor = igraph_vector_size(&parent_neighbors);

					igraph_vector_int_t neighbor_edge_id;
					igraph_vector_int_init(&neighbor_edge_id, n_neighbor);
					for (int j = 0; j < n_neighbor; ++j) {
						igraph_get_eid(&tmp, &VECTOR(neighbor_edge_id)[j], parent, VECTOR(parent_neighbors)[j], IGRAPH_UNDIRECTED, 0);
					}
					igraph_delete_edges(&tmp, igraph_ess_vector(&neighbor_edge_id));

					igraph_graph_list_t components;
					igraph_graph_list_init(&components, 0);
					igraph_decompose(&tmp, &components, IGRAPH_WEAK, -1, 0);

					igraph_t* comp1;
					igraph_t* comp2;

					bool found1{ false }, found2{ false };
					int n_components = igraph_graph_list_size(&components);

					igraph_vector_t sub_vertice_id;
					igraph_vector_init(&sub_vertice_id, n_vertice);

					for (int j = 0; j < n_components; ++j) {
						VANV(&VECTOR(components)[j], "vid", &sub_vertice_id);
						int n_comp_vertice = igraph_vector_size(&sub_vertice_id);

						for (int k = 0; k < n_comp_vertice; ++k) {
							if (!found1) {
								if (VECTOR(sub_vertice_id)[k] == parent1) {
									comp1 = &VECTOR(components)[j];
									found1 = true;
								}
							}

							if (!found2) {
								if (VECTOR(sub_vertice_id)[k] == parent2) {
									comp2 = &VECTOR(components)[j];
									found2 = true;
								}
							}

							if (found1 && found2) {
								break;
							}
						}

						if (found1 && found2) {
							break;
						}
					}

					igraph_real_t diameter1;
					igraph_diameter_dijkstra(comp1, 0, &diameter1, 0, 0, 0, 0, IGRAPH_DIRECTED, true);
					igraph_real_t diameter2;
					igraph_diameter_dijkstra(comp2, 0, &diameter2, 0, 0, 0, 0, IGRAPH_DIRECTED, true);

					++diameter1;
					++diameter2;

					if (diameter1 < minimal_branch_len) {

						VANV(comp1, "vid", &sub_vertice_id);
						int n_comp_vertice = igraph_vector_size(&sub_vertice_id);

						for (int k = 0; k < n_comp_vertice; ++k) {
							vertice_to_delete.push_back(VECTOR(sub_vertice_id)[k]);
						}
					}

					if (diameter2 < minimal_branch_len) {

						VANV(comp2, "vid", &sub_vertice_id);
						int n_comp_vertice = igraph_vector_size(&sub_vertice_id);

						for (int k = 0; k < n_comp_vertice; ++k) {
							vertice_to_delete.push_back(VECTOR(sub_vertice_id)[k]);
						}
					}

					for (int j = 0; j < n_components; ++j) {
						igraph_destroy(&VECTOR(components)[j]);
					}

					igraph_vector_destroy(&sub_vertice_id);
					igraph_vector_destroy(&parent_neighbors);
					igraph_graph_list_destroy(&components);
					igraph_destroy(&tmp);
				}
			}
		}
	}

	vertice_to_delete = custom::unique(vertice_to_delete);
	vertice_to_delete = custom::set_difference(vertice_to_delete, vertice_to_keep);
	
	if (!vertice_to_delete.empty()) {
		int n_vertice_delete = vertice_to_delete.size();
		igraph_vector_int_t deleted_vertices;
		igraph_vector_int_init(&deleted_vertices, n_vertice_delete);

		for (int i = 0; i < n_vertice_delete; ++i) {
			VECTOR(deleted_vertices)[i] = vertice_to_delete[i];
		}

		igraph_delete_vertices(&mst_g, igraph_vss_vector(&deleted_vertices));

		igraph_vector_int_destroy(&deleted_vertices);
	}

	igraph_destroy(&mst_g_old);

	return mst_g;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, QVector<double> > Monocle3Worker::calculate_principal_graph(
	const Eigen::MatrixXd& data,
	const Eigen::MatrixXd& centroids,
	int max_iter,
	double epsilon,
	double l1_gamma,
	double l1_sigma
) {

	Eigen::MatrixXd c = centroids;
	int k = c.cols();

	QVector<double> objs;

	Eigen::MatrixXd W;
	Eigen::MatrixXd P;
	double obj_P;

	for (int iter = 0; iter < max_iter; ++iter) {

		Eigen::ArrayXd sq = c.colwise().squaredNorm();

		Eigen::MatrixXd norm_sq(k, k);
		norm_sq.rowwise() = sq.matrix().transpose();

		Eigen::MatrixXd phi = norm_sq + norm_sq.transpose() - 2 * c.transpose() * c;

		auto stree = mst(phi);

		stree += stree.transpose().eval();

		W = (stree.array() != 0.0).cast<double>();

		double obj_W = stree.sum();

		std::tie(P, obj_P) = soft_assignment(data, c, l1_sigma);

		double obj = obj_W + l1_gamma * obj_P;

		objs << obj;

		if (iter > 0) {
			double diff = abs(objs[iter] - objs[iter - 1]) / abs(objs[iter - 1]);

			if (diff < epsilon) {
				G_TASK_LOG("epsilon = " + QString::number(diff) + " converge.");
				break;
			}

			if (iter == max_iter - 1) {
				G_TASK_NOTICE("epsilon = " + QString::number(diff) + " reach max iteration.");
				break;
			}
		}

		c = generate_centers(data, W, P, l1_gamma);
	}

	return std::make_tuple(c, W, P, objs);
}

void Monocle3Worker::run() {

	int n_cell = this->cell_choosed_.count();

	if (n_cell < 10) {
		G_TASK_WARN("Two few cells for Monocle3.");
		G_TASK_END;
	}


	igraph_set_attribute_table(&igraph_cattribute_table);

	this->learn_graph();

	G_TASK_END;
};

void Monocle3Worker::learn_graph() {

	int n_cell = this->cell_choosed_.count();

	this->embedding_ = custom::row_sliced(this->original_embedding_.data_.mat_, this->cell_choosed_).transpose();

	int n_center = std::round( 150 * std::log10(n_cell));
	if (n_center >= n_cell) {
		n_center = n_cell - 1;
	}

	Eigen::MatrixXd centers = this->embedding_(Eigen::all, custom::integer_linspaced(n_center, 0, n_cell - 1))
		 + Eigen::MatrixXd::Random(this->embedding_.rows(), n_center) * 1e-8;

	auto [kmeans_center, kmeans_cluster] = custom::kmeans_lloyd(this->embedding_, centers, true, 100);

	int n_kmeans_cluster = kmeans_center.cols();

	auto nearest_center = custom::find_nearest(kmeans_center, this->embedding_);

	constexpr int k = 25;

	auto [nn_index, nn_distance] = custom::get_knn_mt<Euclidean, true>(this->embedding_, std::min(k, n_cell - 1));

	nn_index = nn_index(Eigen::all, Eigen::seq(1, nn_index.cols() - 1)).eval();

	nn_distance = nn_distance(Eigen::all, Eigen::seq(1, nn_distance.cols() - 1)).eval();

	Eigen::ArrayXd rho = -nn_distance.rowwise().mean();

	rho = exp(rho).eval();

	QVector<int> high_density_loc;

	for (int i = 0; i < n_kmeans_cluster; ++i) {

		auto index = custom::match(kmeans_cluster, i);

		auto min_index = custom::argmax(custom::reordered(rho, index));

		high_density_loc << index[min_index];
	}

	auto medioids = this->embedding_(Eigen::all, high_density_loc);

	auto [pr_node_embedding, graph_mat, probability, objective_vals] = this->calculate_principal_graph(
		this->embedding_,
		medioids,
		this->max_iter_,
		this->epsilon_,
		this->l1_gamma_,
		this->l1_sigma_
		);
		

	this->pr_graph_ = connect_tips(
		graph_mat,
		kmeans_cluster,
		pr_node_embedding,
		this->embedding_,
		10
	);

	int n_vertice_remain = igraph_vcount(&this->pr_graph_);
	QVector<int> remain_id(n_vertice_remain);
	igraph_vector_t remain_id_i;
	igraph_vector_init(&remain_id_i, n_vertice_remain);
	VANV(&this->pr_graph_, "vid", &remain_id_i);

	for (int i = 0; i < n_vertice_remain; ++i) {
		remain_id[i] = VECTOR(remain_id_i)[i];
	}
	igraph_vector_destroy(&remain_id_i);

	pr_node_embedding = pr_node_embedding(Eigen::all, remain_id).eval();
	//probability = probability(Eigen::all, remain_id).eval();
	//medioids = medioids(Eigen::all, remain_id).eval();

	this->project_to_mst(false, pr_node_embedding);

};

static Eigen::ArrayXd proj_point_to_line_segment(const Eigen::ArrayXd& p, const Eigen::ArrayXd& from, const Eigen::ArrayXd& to) {
	
	Eigen::ArrayXd ab = to - from;
	double ab_squared = ab.matrix().squaredNorm();

	if (ab_squared == 0) {
		return from;
	}
	else {
		Eigen::ArrayXd ap = p - from;

		double t = (ap * ab).sum() / ab_squared;

		if (t < 0.0) {
			return from;
		}

		if (t > 1.0) {
			return to;
		}

		return from + t * ab;
	}

}

static Eigen::ArrayXd proj_point_on_line(const Eigen::ArrayXd& p, const Eigen::ArrayXd& from, const Eigen::ArrayXd& to) {
	Eigen::ArrayXd ap = p - from;
	Eigen::ArrayXd ab = to - from;

	return from + (ap * ab).sum() / (ab * ab).sum() * ab;
}

void Monocle3Worker::project_to_mst(
	bool orthogonal_proj_tip, 
	const Eigen::MatrixXd& pr_node_embedding) 
{

	auto closest_vertex = custom::find_nearest(this->embedding_, pr_node_embedding);
	int n_cell = closest_vertex.size();

	igraph_vector_int_t degree;
	igraph_vector_int_init(&degree, 0);
	igraph_degree(&this->pr_graph_, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

	Eigen::MatrixXd P(2, n_cell);
	Eigen::MatrixXi nearest_edges(n_cell, 2);

	for (int i = 0; i < n_cell; ++i) {

		igraph_vector_t vertex_neighbors;
		igraph_vector_init(&vertex_neighbors, 0);
		igraph_i_neighbors(&this->pr_graph_, &vertex_neighbors, closest_vertex[i], IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
				
		int n_neighbor = igraph_vector_size(&vertex_neighbors);

		Eigen::MatrixXd projection(n_neighbor, 2);
		Eigen::ArrayXd distances(n_neighbor);
		for (int j = 0; j < n_neighbor; ++j) {
			if (VECTOR(degree)[closest_vertex[i]] == 1) {
				if (orthogonal_proj_tip) {
					projection.row(j) = proj_point_on_line(
						this->embedding_.col(i),
						pr_node_embedding.col(closest_vertex[i]),
						pr_node_embedding.col(VECTOR(vertex_neighbors)[j]));
				}
				else {
					projection.row(j) = proj_point_to_line_segment(
						this->embedding_.col(i),
						pr_node_embedding.col(closest_vertex[i]),
						pr_node_embedding.col(VECTOR(vertex_neighbors)[j]));
				}
			}
			else {
				projection.row(j) = proj_point_to_line_segment(
					this->embedding_.col(i),
					pr_node_embedding.col(closest_vertex[i]),
					pr_node_embedding.col(VECTOR(vertex_neighbors)[j]));
			}

			distances[j] = (projection.row(j) - this->embedding_.col(i)).squaredNorm();
		}
		auto which_min = custom::argmin(distances);
		P.col(i) = projection.row(which_min);
		nearest_edges(i, 0) = closest_vertex[i];
		nearest_edges(i, 1) = VECTOR(vertex_neighbors)[which_min];
		igraph_vector_destroy(&vertex_neighbors);
	}

	Eigen::ArrayXi source = nearest_edges.col(0);
	Eigen::ArrayXi target = nearest_edges.col(1);

	for (int i = 0; i < n_cell; ++i) {
		if (source[i] > target[i]) {
			std::swap(source[i], target[i]);
		}
	}

	Eigen::ArrayXd distance_to_source(n_cell);
	for (int i = 0; i < n_cell; ++i) {
		distance_to_source[i] = (P.col(i) - pr_node_embedding.col(source[i])).squaredNorm();
	}

	QVector<std::pair<int, int>> edges;
	for (int i = 0; i < n_cell; ++i) {
		edges << std::make_pair(source[i], target[i]);
	}

	Eigen::ArrayXi new_cell_order(n_cell);
	int loc_now = 0;

	auto pr_edge = custom::unique(edges);
	for (auto&& edge : pr_edge) {

		auto edge_index = custom::match(edges, edge);
		Eigen::ArrayXd distance = distance_to_source(edge_index);
		auto dis_order = custom::order(distance);

		int n_sub_cell = edge_index.size();
		new_cell_order.segment(loc_now, n_sub_cell) = custom::cast<Eigen::ArrayX>(custom::reordered(edge_index, dis_order));
		
		loc_now += n_sub_cell;
	}

	edges = custom::reordered(edges, new_cell_order);

	source = source(new_cell_order).eval();
	target = target(new_cell_order).eval();

	int n_pr_vertice = igraph_vcount(&this->pr_graph_);
	Eigen::ArrayXi cell_graph_id = new_cell_order + n_pr_vertice;


	QVector<int> cell_edges;
	QVector<double> cell_weights;

	std::map<std::pair<int, int>, int> pr_edge_weight_adj;

	for (auto&& edge : pr_edge) {
		auto edge_index = custom::match(edges, edge);
		Eigen::ArrayXi cell_id = cell_graph_id(edge_index);
		Eigen::ArrayXi cell_order = new_cell_order(edge_index);
		int n_cell_edge = cell_id.size();
		int s = source[edge_index[0]];
		int t = target[edge_index[0]];

		cell_edges << s << cell_id[0];
		cell_weights << (pr_node_embedding.col(s) - P.col(cell_order[0])).norm();

		for (int i = 0; i < n_cell_edge - 1; ++i) {
			cell_edges << cell_id[i] << cell_id[i + 1];
			cell_weights << (P.col(cell_order[i]) - P.col(cell_order[i + 1])).norm();
		}

		cell_edges << cell_id[n_cell_edge - 1] << t;
		cell_weights << (P.col(cell_order[n_cell_edge - 1]) - pr_node_embedding.col(t)).norm();

		pr_edge_weight_adj[std::make_pair(s, t)] = n_cell_edge + 1;
	}

	igraph_vector_int_t pr_graph_edges;
	igraph_vector_int_init(&pr_graph_edges, 0);
	igraph_get_edgelist(&this->pr_graph_, &pr_graph_edges, 0);

	int n_cell_edge = cell_weights.size();
	int n_pr_graph_edge = igraph_vector_int_size(&pr_graph_edges) / 2;

	QVector<double> pr_weights;
	for (int i = 0; i < n_pr_graph_edge; ++i) {
		pr_weights << (pr_node_embedding.col(VECTOR(pr_graph_edges)[2 * i]) - pr_node_embedding.col(VECTOR(pr_graph_edges)[2 * i + 1])).norm();
	}
	int n_all_edge = n_pr_graph_edge + n_cell_edge;

	igraph_vector_int_t all_edges;
	igraph_vector_int_init(&all_edges, n_all_edge * 2);

	igraph_vector_t all_weights;
	igraph_vector_init(&all_weights, n_all_edge);

	for (int i = 0; i < n_pr_graph_edge; ++i) {
		VECTOR(all_edges)[2 * i] = VECTOR(pr_graph_edges)[2 * i];
		VECTOR(all_edges)[2 * i + 1] = VECTOR(pr_graph_edges)[2 * i + 1];
		VECTOR(all_weights)[i] = pr_weights[i];
	}

	for (int i = 0; i < n_cell_edge; ++i) {
		VECTOR(all_edges)[2 * i + n_pr_graph_edge * 2] = cell_edges[2 * i];
		VECTOR(all_edges)[2 * i + n_pr_graph_edge * 2 + 1] = cell_edges[2 * i + 1];
		VECTOR(all_weights)[i + n_pr_graph_edge] = cell_weights[i];
	}
	
	//double min_weight{ 0.0 };
	//for (int i = 0; i < n_all_edge; ++i) {
	//	double weight = VECTOR(all_weights)[i];
	//	if (weight > 0.0) {
	//		if (min_weight == 0.0 || min_weight > weight) {
	//			min_weight = weight;
	//		}
	//	}
	//}

	//for (int i = 0; i < n_pr_graph_edge; ++i) {
	//	VECTOR(all_weights)[i] += min_weight * pr_edge_weight_adj[std::make_pair(VECTOR(all_edges)[2 * i], VECTOR(all_edges)[2 * i + 1])];
	//}

	//for (int i = n_pr_graph_edge; i < n_all_edge; ++i) {
	//	VECTOR(all_weights)[i] += min_weight;
	//}

	igraph_empty(&this->cell_graph_, n_pr_vertice + n_cell, IGRAPH_UNDIRECTED);

	igraph_add_edges(&this->cell_graph_, &all_edges, NULL);

	Monocle3* m = new Monocle3();
	igraph_vector_init_copy(&m->cell_graph_weights_, &all_weights);
	m->original_embedding_ = this->original_embedding_;
	m->cell_embedding_ = this->embedding_;
	m->pr_embedding_ = pr_node_embedding;
	m->cell_included_ = this->cell_choosed_;
	igraph_copy(&m->cell_graph_, &this->cell_graph_);
	igraph_copy(&m->pr_graph_, &this->pr_graph_);

	emit x_monocle3_ready(m);
};