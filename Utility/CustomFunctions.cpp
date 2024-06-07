#include "CustomFunctions.h"

#include <sstream>
#include <QDebug>
#include <iostream>
#include <zlib.h>
#include <QRegularExpression>

#include "CustomTemplates.h"

#include "pval.h"

namespace custom {

	Eigen::ArrayXd row_var_mt(const Eigen::MatrixXd& mat) {
		const int ncol = mat.cols(), nrow = mat.rows();
		Eigen::ArrayXd means = mat.rowwise().mean();
		Eigen::ArrayXd vars(nrow);
	#pragma omp parallel for
		for (int i = 0; i < nrow; ++i) {
			vars[i] = (mat.row(i).array() - means[i]).square().sum() / (ncol - 1);
		}
		return vars;
	}

	Eigen::ArrayXd col_var_mt(const Eigen::MatrixXd& mat) {
		const int ncol = mat.cols(), nrow = mat.rows();
		Eigen::ArrayXd means = mat.colwise().mean();
		Eigen::ArrayXd vars(ncol);
	#pragma omp parallel for
		for (int i = 0; i < ncol; ++i) {
			vars[i] = (mat.col(i).array() - means[i]).square().sum() / (nrow - 1);
		}
		return vars;
	}

	double sd(const Eigen::MatrixXd& mat) {
	
		const int nrow = mat.rows();
		const int ncol = mat.cols();

		double mean = mat.mean();
		double val{ 0.0 };

		for (int j = 0; j < ncol; ++j) {
			for (int i = 0; i < nrow; ++i) {
				val += (mat(i, j) - mean) * (mat(i, j) - mean);
			}
		}

		return std::sqrt(val / (nrow * ncol - 1));
	};

	StringSplitView split_view(const std::string& s, char delimiter) {
	
		return StringSplitView(s, delimiter);
	};

	std::vector<std::string> split(std::string_view str, char delimiter) {
		std::vector<std::string> result;
		size_t start = 0;
		size_t end = 0;

		while ((end = str.find(delimiter, start)) != std::string_view::npos) {
			result.emplace_back(str.substr(start, end - start));
			start = end + 1;
		}
		result.emplace_back(str.substr(start));

		return result;
	}

	Eigen::ArrayXi cluster_louvain_igraph(const Eigen::SparseMatrix<double>& adj, double cutoff) {

		Eigen::SparseMatrix<double> filtered = _Cs eliminate_less_than(adj, cutoff);

		int n_con = filtered.nonZeros();

		int n_vertice = filtered.rows();

		if (n_con == 0) {

			return Eigen::ArrayXi::LinSpaced(n_vertice, 0, n_vertice - 1);
		}

		int n_edge = (n_con + n_vertice) / 2;

		igraph_t g;
		igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
		igraph_add_vertices(&g, n_vertice, NULL);

		igraph_vector_t weights;
		igraph_vector_init(&weights, 0);

		igraph_sparsemat_t tmp, graph_mat;
		igraph_sparsemat_init(&tmp, n_vertice, n_vertice, n_edge);

		for (int i = 0; i < n_vertice; ++i) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(filtered, i); it; ++it) {

				int row = it.row();
				int col = it.col();
				double val = it.value();

				if (row >= col) {
					break;
				}

				igraph_sparsemat_entry(&tmp, row, col, val);
				igraph_sparsemat_entry(&tmp, col, row, val);
			}
		}
		igraph_sparsemat_compress(&tmp, &graph_mat);
		igraph_sparsemat_destroy(&tmp);

		igraph_sparse_weighted_adjacency(&g, &graph_mat, IGRAPH_ADJ_UNDIRECTED, &weights, IGRAPH_LOOPS_ONCE);

		igraph_real_t resolution = 1.0;
		igraph_vector_int_t membership;
		igraph_matrix_int_t memberships;
		igraph_vector_t modularity;

		igraph_vector_init(&modularity, 0);
		igraph_vector_int_init(&membership, 0);
		igraph_matrix_int_init(&memberships, 0, 0);

		igraph_community_multilevel(&g, &weights, resolution, &membership, &memberships, &modularity);

		igraph_matrix_int_destroy(&memberships);
		igraph_vector_destroy(&modularity);
		igraph_destroy(&g);
		igraph_sparsemat_destroy(&graph_mat);
		igraph_vector_destroy(&weights);

		Eigen::ArrayXi cluster(n_vertice);

		for (int i = 0; i < n_vertice; ++i) {
			cluster[i] = VECTOR(membership)[i];
		}

		igraph_vector_int_destroy(&membership);

		return cluster;
	}

	void sort_by_first(Eigen::ArrayXd& first, Eigen::MatrixXd& second, bool decrease) {

		Eigen::ArrayXi ind = _Cs order(first, decrease);

		first = first(ind).eval();

		second = second(ind, Eigen::all).eval();
	};

	Eigen::SparseMatrix<double> set_from_triplets_mean(const std::vector<Eigen::Triplet<double>>& triplets, int nrow, int ncol) {

		std::vector<Eigen::Triplet<int>> count_triplets;
		count_triplets.reserve(triplets.size());
		for (auto&& t : triplets) {
			count_triplets.emplace_back(t.row(), t.col(), 1);
		}

		Eigen::SparseMatrix<int> count_matrix(nrow, ncol);
		Eigen::SparseMatrix<double> val_matrix(nrow, ncol);

		val_matrix.setFromTriplets(triplets.cbegin(), triplets.cend());
		count_matrix.setFromTriplets(count_triplets.cbegin(), count_triplets.cend());

		for (int i = 0; i < ncol; ++i) {
			Eigen::SparseMatrix<double>::InnerIterator val_it(val_matrix, i);
			Eigen::SparseMatrix<int>::InnerIterator count_it(count_matrix, i);
			for (; val_it; ++val_it, ++count_it) {
				int count = count_it.value();

				if (count != 1) {
					val_it.valueRef() /= count;
				}
			}
		}

		return val_matrix;
	}

	Eigen::ArrayXi find_nearest(const Eigen::MatrixXd& query, const Eigen::MatrixXd& data, bool by_col) {

		Eigen::ArrayXi nearest;

		if (by_col) {

			int n_query = query.cols(), n_data = data.cols();

			nearest = Eigen::ArrayXi::Zero(n_query);

			for (int i = 0; i < n_query; ++i) {

				Eigen::ArrayXd distance = (data.colwise() - query.col(i)).colwise().squaredNorm();

				nearest[i] = _Cs argmin(distance);

			}
		}
		else {
			int n_query = query.rows(), n_data = data.rows();

			nearest = Eigen::ArrayXi::Zero(n_query);

			for (int i = 0; i < n_query; ++i) {

				Eigen::ArrayXd distance = (data.rowwise() - query.row(i)).rowwise().squaredNorm();

				nearest[i] = _Cs argmin(distance);

			}
		}

		return nearest;

	}

	Eigen::MatrixXd cov2cor(const Eigen::MatrixXd& mat) {
		Eigen::ArrayXd Is = mat.diagonal();

		Is = sqrt(1.0 / Is);

		Eigen::MatrixXd cor = (mat.array().colwise() * Is).rowwise() * Is.transpose();

		cor.diagonal().array() = 1.0;

		return cor;
	}

	Eigen::MatrixXd cov(const Eigen::MatrixXd& mat) {
	
		int n_ob = mat.cols();
		int n_var = mat.rows();

		Eigen::MatrixXd ret(n_ob, n_ob);

		Eigen::MatrixXd tmp = mat.rowwise() - mat.colwise().mean();

		for (int i = 0; i < n_ob; ++i) {
			for (int j = 0; j <= i; ++j) {
				double val = (tmp.col(i).array() * tmp.col(j).array()).sum() / (n_var - 1);

				ret(i, j) = val;
				ret(j, i) = val;
			}
		}

		return ret;
	};

	Eigen::ArrayXd mean_compress(const Eigen::ArrayXd& arr, int time) {
		int size = arr.size();
		int compress_size{ 0 };

		Eigen::ArrayXd res;

		if (size % time == 0) {
			compress_size = size / time;

			res = Eigen::ArrayXd::Zero(compress_size);

			for (int i = 0; i < size; ++i) {
				res[i / time] += arr[i];
			}

			return res / (double)time;

		}
		else {
			compress_size = size / time + 1;

			res = Eigen::ArrayXd::Zero(compress_size);

			for (int i = 0; i < compress_size - 1; ++i) {
				for (int j = 0; j < time; ++j) {
					res[i] += arr[i * time + j];
				}

				res[i] = res[i] / (double)time;
			}

			res[compress_size - 1] = arr.segment((compress_size - 1) * time, size - (compress_size - 1) * time).mean();

			return res;
		}
	};

	void igraph_vector_copy(igraph_vector_t& to, const igraph_vector_t& from) {
		if (to.stor_begin != NULL) {
			igraph_vector_destroy(&to);
		}

		igraph_vector_init_copy(&to, &from);
	};

	void igraph_vector_int_copy(igraph_vector_int_t& to, const igraph_vector_int_t& from) {
		if (to.stor_begin != NULL) {
			igraph_vector_int_destroy(&to);
		}

		igraph_vector_int_init_copy(&to, &from);
	};

	double sign(double val) {
		return val > 0.0 ? 1.0 : (val < 0.0 ? -1.0 : 0.0);
	};

	QString standardize_windows_file_name(const QString& file_name) {
		QString standardized(file_name);

		for (auto&& ch : standardized) {
			if (ch == '|'
				|| ch == '>'
				|| ch == '<'
				|| ch == '*'
				|| ch == ':'
				|| ch == '\\'
				|| ch == '/'
				|| ch == '?') {
				ch = '_';
			}
		}

		return standardized;
	};

	QString complementary_strand(const QString& sequence) {
		QString ret(sequence);

		for (auto& ch : ret) {
			if (ch == 'A') {
				ch = 'T';
			}
			else if (ch == 'T') {
				ch = 'A';
			}
			else if (ch == 'G') {
				ch = 'C';
			}
			else if (ch == 'C') {
				ch = 'G';
			}
			else if (ch == 'a') {
				ch = 't';
			}
			else if (ch == 't') {
				ch = 'a';
			}
			else if (ch == 'g') {
				ch = 'c';
			}
			else if (ch == 'c') {
				ch = 'g';
			}
		}
		return ret;
	};

	QString standardize_chromosome_name(const QString& chr) {
		if (chr.size() > 3 && (chr.startsWith("chr") || chr.startsWith("Chr") || chr.startsWith("CHR"))) {
			return chr.sliced(3);
		}
		else {
			return chr;
		}
	};

	std::string standardize_chromosome_name(const std::string& chr) {
		if (chr.size() > 3 && (chr.starts_with("chr") || chr.starts_with("Chr") || chr.starts_with("CHR"))) {
			return chr.substr(3);
		}
		else {
			return chr;
		}
	};

	QStringList paste(const QString& p1, const QStringList& p2, const QString& concat) {

		bool has_concat = !concat.isEmpty();

		const std::size_t size = p2.size();
		QStringList ret;
		ret.reserve(size);

		auto iter = p2.cbegin();
		auto end = p2.cend();

		if (has_concat) {
			for (; iter != end; ++iter) {
				ret.emplace_back(p1 + concat + *iter);
			}
		}
		else {
			for (; iter != end; ++iter) {
				ret.emplace_back(p1 + *iter);
			}
		}

		return ret;
	};

	QStringList paste(const QStringList& p1, const QString& p2, const QString& concat) {

		bool has_concat = !concat.isEmpty();

		const std::size_t size = p1.size();
		QStringList ret;
		ret.reserve(size);

		auto iter = p1.cbegin();
		auto end = p1.cend();

		if (has_concat) {
			for (; iter != end; ++iter) {
				ret.emplace_back(*iter + concat + p2);
			}
		}
		else {
			for (; iter != end; ++iter) {
				ret.emplace_back(*iter + p2);
			}
		}

		return ret;
	};

	QStringList paste(const QStringList& p1, const QStringList& p2, const QString& concat) {

		bool has_concat = !concat.isEmpty();

		const std::size_t size = p1.size();
		QStringList ret;
		ret.reserve(size);

		auto iter = p1.cbegin();
		auto iter2 = p2.cbegin();
		auto end = p1.cend();

		if (has_concat) {
			for (; iter != end; ++iter, ++iter2) {
				ret.emplace_back(*iter + concat + *iter2);
			}
		}
		else {
			for (; iter != end; ++iter, ++iter2) {
				ret.emplace_back(*iter + *iter2);
			}
		}

		return ret;
	};

	QVector<int> seq_n(int start, int size, int incr) {

		QVector<int> ret(size);

		for (int i = 0; i < size; ++i) {

			ret[i] = start;
			
			start += incr;
		}

		return ret;
	};

	void rowwise_divide_in_place(Eigen::SparseMatrix<double>& dividend, const Eigen::ArrayXd& divisor) {
		for (int k = 0; k < dividend.outerSize(); ++k) {

			double div = divisor[k];

			if (div == 0) {
				continue;
			}

			for (Eigen::SparseMatrix<double>::InnerIterator it(dividend, k); it; ++it) {
				it.valueRef() /= div;
			}
		}
	};

	Eigen::MatrixXi balanced_knn_mt(
		const Eigen::MatrixXd& data,
		int n_neighbor,
		int maxl,
		bool var_by_column,
		const QString& method
	) {
		Eigen::MatrixXd distance;
		if (method == "cor") {
			distance = _Cs correlation_pearson_distance_mt(data, var_by_column);
			distance = 1.0 - distance.array();
		}
		else {
			distance = _Cs euclidean_distance_mt(data, var_by_column);
		}
		const int n_var = distance.cols();
		Eigen::ArrayXi l = Eigen::ArrayXi::Zero(n_var); // reciprocal neighbor count
		Eigen::MatrixXi dsi = Eigen::MatrixXi::Zero(n_var, n_var); // will store column sort indices of the d matrix

	#pragma omp parallel for
		for(int i = 0; i < n_var; ++i) {
			auto order = _Cs order(distance.col(i));
			dsi.col(i) = order;

		#pragma omp critical
			{
				for (int j = 0; j < n_neighbor; ++j) {
					++l[order[j]];
				}
			}
		};

		Eigen::ArrayXi lsi = _Cs order(l, true);// greedy order of column considerations
		l.setZero(); // reset so that l can be used to keep track of reciprocal counts in the kNN being constucted
		// greedy knn (construct row index vector)

		Eigen::MatrixXi knn(n_var, n_neighbor);

		// run regular knn, record number of reciprocal neighbors
		for (int i = 0; i < n_var; ++i) {
			int el = lsi[i]; // element for which the neighbors are being found
			Eigen::ArrayXi si = dsi.col(el);
			int p = 0;
			int j;
			for (j = 0; j < n_var && p < n_neighbor; ++j) {
				// consider si[j] for p-th neighbor position
				int m = si[j];
				if (el == m) { continue; } // dont record or count self-relationships
				if (l[m] >= maxl) { continue; } // element has already maxed out its neighbors
				knn(el, p) = m;
				++l[m];
				++p; // record neighbor
			}
			if (j == n_var && p < n_neighbor) {
				// fill in last element(s) with self-identity if there were not enough spare elements 
				while (p < n_neighbor) {
					knn(el, p) = el;
					++p;
				}
			}
		}

		return knn;
	}

	Eigen::ArrayXd cor(const Eigen::MatrixXd& mat, const Eigen::ArrayXd& arr) {
		const int ncol = mat.cols();
		Eigen::ArrayXd correlation(ncol);

		for (int i = 0; i < ncol; ++i) {
			correlation[i] = _Cs correlation_pearson(mat.col(i), arr);
		}

		return correlation;
	};

	Eigen::ArrayXd cor_mt(const Eigen::MatrixXd& mat, const Eigen::ArrayXd& arr) {
		const int ncol = mat.cols();
		Eigen::ArrayXd correlation(ncol);

	#pragma omp parallel for
		for (int i = 0; i < ncol; ++i) {
			correlation[i] = _Cs correlation_pearson(mat.col(i), arr);
		}

		return correlation;
	};

	Eigen::MatrixXd distance(const QVector<double>& vec) {
		
		int n_node = vec.size();
		Eigen::MatrixXd dist(n_node, n_node);
		for (int i = 0; i < n_node; ++i) {
			for (int j = 0; j < n_node; ++j) {
				dist(i, j) = std::abs(vec[i] - vec[j]);
			}
		}

		return dist;
	};

	Eigen::MatrixXd cor(const Eigen::MatrixXd& data, bool var_by_column) {
		if (var_by_column) {
			const int n_var = data.cols(), var_length = data.rows();
			Eigen::MatrixXd correlation = Eigen::MatrixXd::Zero(n_var, n_var);

			std::size_t n_element = std::size_t(n_var) * var_length;
			for (int i = 0; i < n_var; ++i) {
				correlation(i, i) = 0.5;
				for (int j = 0; j < i; ++j) {
					correlation(i, j) = _Cs correlation_pearson(data.col(i), data.col(j));
				}
			}

			return correlation + correlation.transpose();
		}
		else {
			const int n_var = data.rows(), var_length = data.cols();
			Eigen::MatrixXd correlation = Eigen::MatrixXd::Zero(n_var, n_var);

			std::size_t n_element = std::size_t(n_var) * var_length;
			for (int i = 0; i < n_var; ++i) {
				correlation(i, i) = 0.5;
				for (int j = 0; j < i; ++j) {
					correlation(i, j) = _Cs correlation_pearson(data.row(i), data.row(j));
				}
			}

			return correlation + correlation.transpose();
		}
	};
	
	Eigen::MatrixXd correlation_pearson_distance_mt(const Eigen::MatrixXd& data, bool var_by_column) {
		if (var_by_column) {
			const int n_var = data.cols(), var_length = data.rows();
			Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(n_var, n_var);

			std::size_t n_element = std::size_t(n_var) * var_length;
			if (n_element < 1e6) {
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = 1 - _Cs correlation_pearson(data.col(i), data.col(j));
					}
				}
			}
			else {

			#pragma omp parallel for
				for(int i = 0; i < n_var; ++i){
					for (int j = 0; j < i; ++j) {
						distance(i, j) = 1 - _Cs correlation_pearson(data.col(i), data.col(j));
					}
				}
			}

			return distance + distance.transpose();
		}
		else {
			const int n_var = data.rows(), var_length = data.cols();
			Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(n_var, n_var);

			std::size_t n_element = std::size_t(n_var) * var_length;
			if (n_element < 1e6) {
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = 1 - _Cs correlation_pearson(data.row(i), data.row(j));
					}
				}
			}
			else {

			#pragma omp parallel for
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = 1 - _Cs correlation_pearson(data.row(i), data.row(j));
					}
				};
			}

			return distance + distance.transpose();
		}
	}

	Eigen::MatrixXd euclidean_distance_mt(const Eigen::MatrixXd& data, bool var_by_column) {
		if (var_by_column) {
			const int n_var = data.cols(), var_length = data.rows();
			Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(n_var, n_var);

			std::size_t n_element = std::size_t(n_var) * var_length;
			if (n_element < 1e6) {
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = std::sqrt((data.col(i).array() - data.col(j).array()).square().sum());
					}
				}
			}
			else {

			#pragma omp parallel for
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = std::sqrt((data.col(i).array() - data.col(j).array()).square().sum());
					}
				};
			}

			return distance + distance.transpose();
		}
		else {
			const int n_var = data.rows(), var_length = data.cols();
			Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(n_var, n_var);

			std::size_t n_element = std::size_t(n_var) * var_length;
			if (n_element < 1e6) {
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = std::sqrt((data.row(i).array() - data.row(j).array()).square().sum());
					}
				}
			}
			else {
			#pragma omp parallel for
				for (int i = 0; i < n_var; ++i) {
					for (int j = 0; j < i; ++j) {
						distance(i, j) = std::sqrt((data.row(i).array() - data.row(j).array()).square().sum());
					}
				};
			}

			return distance + distance.transpose();
		}
	}

	Eigen::SparseMatrix<int> create_matrix_from_knn_index(const Eigen::MatrixXi& knn) {
		std::vector<Eigen::Triplet<int> > triplets;
		const int nrow = knn.rows(), ncol = knn.cols(), n_nodes = knn.rows(), n_neighbors = knn.cols();

		triplets.reserve(nrow * ncol);

		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {
				triplets.emplace_back(i, knn(i, j), 1);
			}
		}

		Eigen::SparseMatrix<int> nn(n_nodes, n_nodes);
		nn.reserve(nrow * ncol);
		nn.setFromTriplets(triplets.cbegin(), triplets.cend());

		triplets.clear();

		return nn;
	};

	Eigen::SparseMatrix<double> create_shared_nearest_neighbors_matrix(const Eigen::MatrixXi& knn) {
		std::vector<Eigen::Triplet<double> > triplets;
		const int nrow = knn.rows(), ncol = knn.cols(), n_nodes = knn.rows(), n_neighbors = knn.cols();

		triplets.reserve(nrow * ncol);

		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {
				triplets.emplace_back(i, knn(i, j), 1);
			}
		}

		Eigen::SparseMatrix<double> nn(n_nodes, n_nodes);
		nn.reserve(nrow * ncol);
		nn.setFromTriplets(triplets.cbegin(), triplets.cend());

		triplets.clear();

		constexpr double prune_snn = 1.0 / 15;
		Eigen::SparseMatrix<double> snn = nn * (nn.transpose());

		for (int i = 0; i < n_nodes; ++i) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
				it.valueRef() = it.value() / (n_neighbors + (n_neighbors - it.value()));
			}
		}

		snn = _Cs eliminate_less_than(snn, prune_snn);

		return snn;
	};

	std::vector<int> generate_random_permutation(std::size_t n_element, std::default_random_engine& re) {
		std::vector<int> permutation(n_element);
		for (std::size_t i = 0; i < n_element; ++i) {
			permutation[i] = i;
		}
		std::shuffle(permutation.begin(), permutation.end(), re);
		return permutation;
	};

	QStringList multi_split(const QString& src, const QList<QChar>& splits) {
		QStringList res;
		const qsizetype size = src.size();
		if (size == 0) {
			return res;
		}

		qsizetype last_index = 0;
		for (qsizetype i = 0; i < size; ++i) {
			if (splits.contains(src[i])) {
				res << src.sliced(last_index, i - last_index);
				last_index = i + 1;
			}
		}
		res << src.sliced(last_index, size - last_index);
		return res;
	};

	QVector<double> empirical_cumulative_distribution(const QVector<double>& vec) {

		const qsizetype length = vec.size();

		Eigen::ArrayXi index = Eigen::ArrayXi::LinSpaced(length, 0, length - 1);

		QVector<double> ecd(length);

		std::sort(index.begin(), index.end(), [&vec](int i, int j)-> bool {
			return vec[i] < vec[j];
		});

		qsizetype count = 1;

		double current_value = vec[index[0]];
		for (qsizetype i = 1; i < length; ++i) {
			double this_value = vec[index[i]];

			if (current_value != this_value) {
				if (count > 1) {
					double rank_mean = i - (count + 1) / 2.0;
					for (qsizetype j = 0; j < count; ++j) {
						ecd[index[i - j - 1]] = i;
					}
				}
				else {
					ecd[index[i - 1]] = i;
				}
				count = 1;
			}
			else {
				++count;
			}
			current_value = this_value;
		}
		if (count > 1) {
			for (int j = 0; j < count; ++j) {
				ecd[index[length - j - 1]] = length;
			}
		}
		else {
			ecd[index[length - 1]] = length;
		}
		return _Cs divide(ecd, length);
	};

	bool is_integer(const QString& str) {
		bool is_int{ false };
		str.toInt(&is_int);
		return is_int;
	}

	bool is_numeric(const QString& str) {
		bool is_num{ false };
		str.toDouble(&is_num);
		return is_num;
	}

	bool is_integer(const char* str) {

		bool is_int{ false };

		int index = 0;
		for (; *str != '\0'; ++str, ++index)
		{
			switch (*str)
			{
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				is_int = true;
				break;
			case '-':
			case '+':
				if (index != 0)
				{
					return false;
				}
				break;
			default:
				return false;
			}
		}

		return is_int;
	}

	int detect_type(const QStringList& vec) {
		
		bool test_integer{ true };
		bool test_numeric{ true };

		for (auto&& i : vec) {

			if (!test_integer && !test_numeric) {
				break;
			}

			if (test_integer) {
				if (!_Cs is_integer(i)) {
					test_integer = false;
				}
			}

			if (test_numeric) {
				if (!_Cs is_numeric(i)) {
					test_numeric = false;
				}
			}
		}

		if (test_integer) {
			return 2;
		}

		if (test_numeric) {
			return 1;
		}

		return 0;
	};

	bool is_integer(const QStringList& vec) {

		for (auto&& i : vec) {
			if (!_Cs is_integer(i)) {
				return false;
			}
		}

		return true;
	}

	bool is_numeric(const QStringList& vec) {
		for (auto&& i : vec) {
			if (!_Cs is_numeric(i)) {
				return false;
			}
		}

		return true;
	}

	bool is_numeric(const char* str) {
		bool is_e = false,
			is_point = false,
			number_before_e = false,
			number_behind_e = false;

		int index = 0;
		for (; *str != '\0'; ++str, ++index)
		{
			switch (*str)
			{
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				if (is_e)
				{
					number_behind_e = true;
				}
				else
				{
					number_before_e = true;
				}
				break;
			case '+':
			case '-':
				if (index != 0)
				{
					return false;
				}
				break;
			case 'e':
			case 'E':
				if (is_e || !number_before_e)
				{
					return false;
				}
				else
				{
					is_point = true;
					index = -1;
					is_e = true;
				}
				break;
			case '.':
				if (is_point)
				{
					return false;
				}
				else
				{
					is_point = true;
				}
				break;
			default:
				return false;
			}
		}

		if (!number_before_e)
		{
			return false;
		}
		else if (is_e && !number_behind_e)
		{
			return false;
		}

		return true;
	}

	int line_length(const char* str) {
		int len = 0;
		while (!std::isspace(*str)) {
			++len;
			++str;
		}
		return len;
	};

	QVector<int> integer_divide(const QVector<int>& dividend, const int divisor) {
		const qsizetype size = dividend.size();
		QVector<int> ret(size, 0);
		for (qsizetype i = 0; i < size; ++i) {
			ret[i] = dividend[i] / divisor;
		}
		return ret;
	}

	int atoi_specialized(const char* src) {
		int res = 0;
		char** str = (char**)&src;
		while (std::isdigit(**str))
		{

			res = (res * 10) + **str - '0';
			++(*str);
		}
		return res;
	};

	int atoi_specialized(const char** src) {
		int res = 0;
		while (!std::isdigit(**src)) {
			++(*src);
		}
		while (std::isdigit(**src))
		{

			res = (res * 10) + **src - '0';
			++(*src);
		}
		return res;
	};

	Eigen::ArrayXXd row_compressed(const Eigen::ArrayXXd& mat, int rate) {
		int ncol = mat.cols(), nrow = ceil(mat.rows() / (double)rate);

		Eigen::ArrayXXd ret(nrow, ncol);

		for (int i = 0; i < nrow - 1; ++i) {
			auto index = _Cs seq_n(i * rate, rate);
			ret.row(i) = mat(index, Eigen::all).colwise().mean();
		}

		auto index = _Cs integer_linspaced((mat.rows() - 1 - (nrow - 1) * rate) + 1, (nrow - 1) * rate, mat.rows() - 1);

		ret.row(nrow - 1) = mat(index, Eigen::all).colwise().mean();

		return ret;
	};

	Eigen::ArrayXXd column_compressed(const Eigen::ArrayXXd& mat, int rate) {
		int nrow = mat.rows(), ncol = ceil(mat.cols() / (double)rate);

		Eigen::ArrayXXd ret(nrow, ncol);

		for (int i = 0; i < ncol - 1; ++i) {
			auto index = _Cs seq_n(i * rate, rate);
			ret.col(i) = mat(Eigen::all, index).rowwise().mean();
		}

		auto index = _Cs integer_linspaced((mat.cols() - 1 - (ncol - 1) * rate) + 1, (ncol - 1) * rate, mat.cols() - 1);

		ret.col(ncol - 1) = mat(Eigen::all, index).rowwise().mean();

		return ret;
	};

	Eigen::ArrayXXd column_compressed(const Eigen::ArrayXXd& mat, const QVector<int>& segment_width) {
		int nrow = mat.rows(), ncol = segment_width.size();

		Eigen::ArrayXXd ret(nrow, ncol);
		int current_col_start = 0;

		for (int i = 0; i < ncol; ++i) {
			int width = segment_width[i];
			ret.col(i) = mat(Eigen::all, Eigen::seqN(current_col_start, width)).rowwise().mean();
			current_col_start += width;
		}

		return ret;
	}

	Eigen::ArrayXXd row_compressed(const Eigen::ArrayXXd& mat, const QVector<int>& segment_width) {
		int ncol = mat.cols(), nrow = segment_width.size();

		Eigen::ArrayXXd ret(nrow, ncol);
		int current_row_start = 0;

		for (int i = 0; i < nrow; ++i) {
			int width = segment_width[i];
			ret.row(i) = mat(Eigen::seqN(current_row_start, width), Eigen::all).colwise().mean();
			current_row_start += width;
		}

		return ret;
	}


	QColor random_color() {
		static std::default_random_engine e(1997);
		static std::uniform_int_distribution<unsigned> u(0, 255);
		return QColor(u(e), u(e), u(e));
	};

	QString string_next_line(const QFontMetrics& font, const QString& text, int label_size, bool is_English)
	{
		const qsizetype size = text.size();

		if (size < 2) {
			return text;
		}
		
		int text_size = font.boundingRect(text).width();
		
		if (text_size > label_size)
		{
			qsizetype position = 0;
			int offset = 0;
			bool segment = false;
			if (is_English) {

				int last_loc{ 0 };

				for (qsizetype i = 0; i < size; ++i) {

					offset += font.boundingRect(text[i]).width();

					if (i > 0) {
						if (text[i] == ' ' || text[i - 1] == ' ') {							
							if (offset > label_size)
							{
								if (last_loc > 0) {
									position = last_loc;
								}
								else {
									position = i;
								}
								segment = true;
								break;
							}

							last_loc = i;
						}
					}
					
				}

				if (!segment) {
					return text;
				}
			}
			else {
				for (qsizetype i = 0; i < size; ++i)
				{
					offset += font.boundingRect(text.at(i)).width();
					if (offset > label_size)
					{
						position = i;
						break;
					}
				}
			}

			if (position == 0) {
				position = 1;
			}

			if (position == size - 1) {
				return text;
			}

			QString string_left_data = text.sliced(0, position);
			QString string_middle_data = text.sliced(position);
			return string_left_data + "\n" + string_next_line(font, string_middle_data, label_size, is_English);
		}
		return text;
	}

	Eigen::SparseMatrix<double> eliminate_less_than(const Eigen::SparseMatrix<double>& mat, double threshold) {
		std::vector<Eigen::Triplet<double> > triplets;
		triplets.reserve(mat.nonZeros());
		for (Eigen::Index i = 0; i < mat.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
				if (it.value() >= threshold) {
					triplets.emplace_back(it.row(), i, it.value());
				}
			}
		}
		Eigen::SparseMatrix<double> f(mat.rows(), mat.cols());
		f.setFromTriplets(triplets.begin(), triplets.end());
		return f;
	}

	[[nodiscard]] Eigen::SparseMatrix<double> eliminate_zero(const Eigen::SparseMatrix<double>& mat) {
		std::vector<Eigen::Triplet<double> > triplets;
		triplets.reserve(mat.nonZeros());
		for (Eigen::Index i = 0; i < mat.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
				if (it.value() != 0.0) {
					triplets.emplace_back(it.row(), i, it.value());
				}
			}
		}
		Eigen::SparseMatrix<double> f(mat.rows(), mat.cols());
		f.setFromTriplets(triplets.begin(), triplets.end());
		return f;
	}

	void digest(
		const std::string& line, 
		const char delimiter, 
		std::size_t start, 
		std::size_t end, 
		std::vector<std::pair<std::size_t, std::size_t>>& line_loc) 
	{

		if (delimiter == '\"') { // not allowed
			return;
		}

		std::size_t loc{ start };
		bool in_quote = (line[start] == '\"');
		for (std::size_t i = start; i < end; ++i) {
			if (line[i] == delimiter) {
				if (in_quote) {
					if (line[i - 1] == '\"') {
						if (i < loc + 2) { // Ensure only a quote is in one segment
							return;
						}
						line_loc.emplace_back(loc + 1, i - 1);
						loc = i + 1;
						if (loc < end && line[loc] == '\"') {
							in_quote = true;
							++i;
						}
						else {
							in_quote = false;
						}
					}
				}
				else {
					line_loc.emplace_back(loc, i);
					loc = i + 1;
					if (loc < end && line[loc] == '\"') {
						in_quote = true;
						++i;
					}
				}
			}
		}

		if (loc == end) {
			line_loc.emplace_back(end, end); // last segment is empty
		}

		if (loc < end) {
			if (in_quote) {
				if (line[end - 1] == '\"') {
					if (end < loc + 2) { // only " in one segment
						return;
					}
					line_loc.emplace_back(loc + 1, end - 1);
				}
				else {
					line_loc.emplace_back(loc + 1, end); //still append last incomplete segment
				}
			}
			else {
				line_loc.emplace_back(loc, end);
			}
		}
	};

	void digest(
		const QString& line,
		const QChar delimiter, 
		qsizetype start, 
		qsizetype end, 
		std::vector<std::pair<qsizetype, qsizetype>>& line_loc
	) {

		if (delimiter == '\"') { // not allowed
			return;
		}

		qsizetype loc{ start };
		bool in_quote = (line[start] == '\"');
		for (qsizetype i = start; i < end; ++i) {
			if (line[i] == delimiter) {
				if (in_quote) {
					if (line[i - 1] == '\"') {
						if (i < loc + 2) { // Ensure only a quote is in one segment
							return;
						}
						line_loc.emplace_back(loc + 1, i - 1);
						loc = i + 1;
						if (loc < end && line[loc] == '\"') {
							in_quote = true;
							++i;
						}
						else {
							in_quote = false;
						}
					}
				}
				else {
					line_loc.emplace_back(loc, i);
					loc = i + 1;
					if (loc < end && line[loc] == '\"') {
						in_quote = true;
						++i;
					}
				}
			}
		}

		if (loc == end) {
			line_loc.emplace_back(end, end); // last segment is empty
		}

		if (loc < end) {
			if (in_quote) {
				if (line[end - 1] == '\"') {
					if (end < loc + 2) { // only " in one segment
						return;
					}
					line_loc.emplace_back(loc + 1, end - 1);
				}
				else {
					line_loc.emplace_back(loc + 1, end); //still append last incomplete segment
				}
			}
			else {
				line_loc.emplace_back(loc, end);
			}
		}
	};

	QStringList digest(const QString& line, const QChar delimiter) {
		QStringList raw = line.split(delimiter);
		QStringList ret;
		QString processing_quote_field;
		bool in_quote = false;

		for (auto&& field : raw) {
			if (in_quote) {
				if (field.endsWith('\"')) {
					in_quote = false;
					processing_quote_field += delimiter + field;
					ret << processing_quote_field.sliced(1, processing_quote_field.size() - 2);
					processing_quote_field.clear();
				}
				else {
					processing_quote_field += delimiter + field;
				}
			}
			else {
				if (field.startsWith('\"')) {
					if (field.endsWith('\"')) {
						ret << field.sliced(1, field.size() - 2);
					}
					else {
						in_quote = true;
						processing_quote_field = field;
					}
				}
				else {
					ret << field;
				}
			}
		}

		if (!processing_quote_field.isEmpty()) {
			ret << processing_quote_field;
		}
		return ret;
	};

	QStringList digest(const QString& line, const QString& delimiter) {
		QStringList raw = line.split(delimiter);
		QStringList ret;
		QString processing_quote_field;
		bool in_quote = false;

		for (auto&& field : raw) {
			if (in_quote) {
				if (field.endsWith('\"')) {
					in_quote = false;
					processing_quote_field += delimiter + field;
					ret << processing_quote_field.sliced(1, processing_quote_field.size() - 2);
					processing_quote_field.clear();
				}
				else {
					processing_quote_field += delimiter + field;
				}
			}
			else {
				if (field.startsWith('\"')) {
					if (field.endsWith('\"')) {
						ret << field.sliced(1, field.size() - 2);
					}
					else {
						in_quote = true;
						processing_quote_field = field;
					}
				}
				else {
					ret << field;
				}
			}
		}

		if (!processing_quote_field.isEmpty()) {
			ret << processing_quote_field;
		}
		return ret;
	};

	QString merge_to_string(const QStringList& vec, const QString& delimiter) {
		QString merged;

		if (vec.isEmpty()) {
			return merged;
		}

		merged = vec[0];
		const qsizetype size = vec.size();
		for (qsizetype i = 1; i < size; ++i) {
			merged.append(delimiter).append(vec[i]);
		}

		return merged;
	}

	QString merge_to_string(const QStringList& vec, const QChar delimiter, bool quote) {

		QString merged;

		if (vec.isEmpty()) {
			return merged;
		}

		if (quote) {
			merged = '\"' + vec[0] + '\"';
			const qsizetype size = vec.size();
			for (qsizetype i = 1; i < size; ++i) {
				merged.append(delimiter).append('\"' + vec[i] + '\"');
			}
		}
		else {
			merged = vec[0];
			const qsizetype size = vec.size();
			for (qsizetype i = 1; i < size; ++i) {
				merged.append(delimiter).append(vec[i]);
			}
		}
		return merged;
	};

	QStringList to_upper(const QStringList& list) {
		const qsizetype size = list.size();
		QStringList ret(size);
		for (qsizetype i = 0; i < size; ++i) {
			ret[i] = list[i].toUpper();
		}
		return ret;
	};

	double signal_to_noise(const Eigen::ArrayXd& X, const Eigen::ArrayXd& Y) {
		double 
			mean_x = X.mean(), 
			mean_y = Y.mean(), 
			std_x = sqrt((X - mean_x).square().sum() / X.size()),
			std_y = sqrt((Y - mean_y).square().sum() / Y.size());

		if (std_x + std_y == 0)
		{
			if (mean_x == mean_y) {
				return 0;
			}
			else if (mean_x > mean_y) {
				return 1;
			}
			else {
				return -1;
			}
		}
		return (mean_x - mean_y) / (std_x + std_y);
	};

	double max(const Eigen::SparseMatrix<double>& mat) {
		const Eigen::Index size = mat.nonZeros();
		if (size == 0) {
			return 0;
		}

		const double* data = mat.valuePtr();
		double max_value = data[0];
		for (Eigen::Index i = 1; i < size; ++i) {
			if (max_value < data[i]) {
				max_value = data[i];
			}
		}
		return max_value;
	};

	QStringList make_unique(const QStringList& list) {
		Eigen::ArrayXi ind = order(list);
		QStringList trans = _Cs reordered(list, ind);
		bool equal = false;
		qsizetype equal_length = 0;
		const qsizetype length = list.size();

		for (qsizetype i = 0; i < length - 1; ++i) {
			if (trans[i] != trans[i + 1]) {
				if (equal) {
					for (qsizetype j = 0; j < equal_length; ++j) {
						bool duplicate = true;
						qsizetype counter = j + 1;
						while (duplicate) {
							QString tmp = trans[i + j + 1 - equal_length] + "-" + QString::number(counter++);
							if (!trans.contains(tmp)) {
								trans[i + j + 1 - equal_length] = tmp;
								duplicate = false;
							}
						}
					}
					equal = false;
					equal_length = 0;
				}
			}
			else {
				++equal_length;
				equal = true;
			}
		}
		if (equal) {
			for (qsizetype j = 0; j < equal_length; ++j) {
				bool duplicate = true;
				qsizetype counter = j + 1;
				while (duplicate) {
					QString tmp = trans[length + j - equal_length] + "-" + QString::number(counter++);
					if (!trans.contains(tmp)) {
						trans[length + j - equal_length] = tmp;
						duplicate = false;
					}
				}
			}
		}
		return _Cs reordered(trans, order(ind));
	};

	QVector<double> linspaced(int length, double begin, double end) {
		Eigen::ArrayXd a = Eigen::ArrayXd::LinSpaced(length, begin, end);
		return _Cs cast<QVector>(a);
	};

	QVector<int> integer_linspaced(int length, int begin, int end) {
		Eigen::ArrayXi a = Eigen::ArrayXi::LinSpaced(length, begin, end);
		return _Cs cast<QVector>(a);
	};


	Eigen::ArrayXd set_upper_bound(const Eigen::ArrayXd& vec, double bound) {
		const Eigen::Index size = vec.size();
		Eigen::ArrayXd ret(size);
		for (qsizetype i = 0; i < size; ++i) {
			double val = vec(i);
			ret[i] = val < bound ? val : bound;
		}
		return ret;
	};

	Eigen::ArrayXd adjust_p_value(const Eigen::ArrayXd& p_value, const QString& method) {
		if (method == "Bonferroni") {
			return _Cs set_upper_bound(p_value * p_value.size(), 1);
		}
		else {
			// fdr
			auto p_order = _Cs order(p_value);
			auto original_order = _Cs order(p_order);
			const Eigen::Index size = p_value.size();
			Eigen::ArrayXd tmp(size);
			for (Eigen::Index i = 0; i < size; ++i) {
				tmp[i] = (double)size / (i + 1);
			}
			return _Cs set_upper_bound((tmp * p_value(p_order))(original_order), 1);
		}
	};

	QVector<double> adjust_p_value(const QVector<double>& p_value, const QString& method) {
		return _Cs cast<QVector>(_Cs adjust_p_value(_Cs cast<Eigen::ArrayX>(p_value), method));
	};

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_lloyd(
		const Eigen::MatrixXd& mat, 
		Eigen::MatrixXd centers, 
		bool remove_empty_cluster,
		int max_iteration
		) 
	{

		int nrow = mat.rows(), ncol = mat.cols(), index_new = 0, n_center = centers.cols();
		bool updated{ false };
		Eigen::ArrayXi cluster = Eigen::ArrayXi::Zero(ncol), cluster_count = Eigen::ArrayXi::Zero(n_center);
		double best, dd;

		for (int iter = 0; iter < max_iteration; ++iter) {

			updated = false;
			
			for (int i = 0; i < ncol; ++i) {

				best = INT_MAX;
				
				for (int j = 0; j < n_center; ++j) {

					dd = (mat.col(i) - centers.col(j)).squaredNorm();
					
					if (dd < best) {
						best = dd;
						index_new = j;
					}
				}

				if (cluster[i] != index_new) {

					updated = true;
					
					cluster[i] = index_new;
				}
			}

			if (!updated) {
			
				break;
			}

			centers.setZero();
			cluster_count.setZero();
			
			for (int i = 0; i < ncol; ++i) {
				
				int index = cluster[i];
				
				++cluster_count[index];

				centers.col(index) += mat.col(i);
			}

			for (int i = 0; i < n_center; ++i) {
			
				if (cluster_count[i] == 0) {
					continue;
				}

				centers.col(i) /= cluster_count[i];
			}

		}

		if (remove_empty_cluster) {
			auto filter = _Cs not_equal(cluster_count, 0);

			if (!_Cs all(filter, true)) {

				auto valid_cluster = _Cs which(filter);

				centers = centers(Eigen::all, valid_cluster).eval();

				int n_valid_cluster = valid_cluster.size();

				QMap<int, int> cluster_map;

				for (int i = 0; i < n_valid_cluster; ++i) {
					cluster_map[valid_cluster[i]] = i;
				}

				std::ranges::for_each(cluster, [&cluster_map](auto& clu) { clu = cluster_map[clu]; });

			}
		}

		return std::make_pair(centers, cluster);
	};

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_lloyd(
		const Eigen::MatrixXd& mat, 
		int k, 
		bool remove_empty_cluster,
		int max_iteration, 
		int random_state) 
	{

		Eigen::MatrixXd centers = mat(Eigen::all, sample_integer(0, mat.cols() - 1, k, random_state));

		return _Cs kmeans_lloyd(mat, centers, remove_empty_cluster, max_iteration);
	};

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_hartigan_wong_once(
		const Eigen::MatrixXd& mat, 
		int k, 
		int max_iteration, 
		int random_state) 
	{
		int ncol = mat.cols(), nrow = mat.rows();
		Eigen::ArrayXi cluster = Eigen::ArrayXi::Zero(ncol);
		Eigen::ArrayXi cluster_count = Eigen::ArrayXi::Zero(k);
		Eigen::ArrayXi cluster2 = Eigen::ArrayXi::Zero(ncol);
		Eigen::MatrixXd centers;

		bool flag = true;
		int loop = 0;

		while (flag) {
			flag = false;
			centers = mat(Eigen::all, sample_integer(0, ncol - 1, k, random_state + loop));

			for (int i = 0; i < ncol; ++i) {
				auto& best = cluster[i];
				best = 0;
				double best_distance = (mat.col(i) - centers.col(best)).squaredNorm();

				auto& second = cluster2[i];
				second = 1;
				double second_distance = (mat.col(i) - centers.col(second)).squaredNorm();

				if (best_distance > second_distance) {
					std::swap(best, second);
					std::swap(best_distance, second_distance);
				}

				for (int j = 2; j < k; ++j) {
					double candidate_dist = (mat.col(i) - centers.col(j)).squaredNorm();
					if (candidate_dist < second_distance) {
						second_distance = candidate_dist;
						second = j;
						if (candidate_dist < best_distance) {
							std::swap(best_distance, second_distance);
							std::swap(best, second);
						}
					}
				}
			}
			centers.setZero();
			cluster_count.setZero();
			for (int i = 0; i < ncol; ++i) {
				int index = cluster[i];
				cluster_count[index] ++;
				centers.col(index) += mat.col(i);
			}

			for (int i = 0; i < k; ++i) {
				if (cluster_count[i] == 0) {
					flag = true;
					loop++;
					break;
				}
				centers.col(i) /= cluster_count[i];
			}
		}

		Eigen::ArrayXd an1(k), an2(k), d(ncol);

		double big = 1e30;

		for (int i = 0; i < k; ++i) {
			double count = cluster_count[i];
			an2[i] = count / (count + 1);
			an1[i] = (count > 1 ? count / (count - 1) : big);
		}

		int imaxqtr = ncol * 50;

		Eigen::ArrayXi ncp = Eigen::ArrayXi::Zero(k), live = ncp;
		Eigen::ArrayX<bool> itran = Eigen::ArrayX<bool>::Constant(k, true);

		int index = 0;
		for (int iter = 0; iter < max_iteration; ++iter) {
			for (int i = 0; i < k; ++i) {
				if (itran[i]) {
					live[i] = ncol;
				}
			}

			for (int i = 0; i < ncol; ++i) {
				++index;
				int l1 = cluster[i];
				if (cluster_count[l1] != 1) {
					if (ncp[l1] != 1) {
						d[i] = (mat.col(i) - centers.col(l1)).squaredNorm() * an1[l1];
					}
					int l2 = cluster2[i];
					int ll = l2;
					double r2 = (mat.col(i) - centers.col(l2)).squaredNorm() * an2[l2];

					for (int j = 0; j < k; ++j) {
						if ((i >= live[l1] && i >= live[j]) || j == l1 || j == ll) {
							continue;
						}
						double rr = r2 / an2[j];
						double dc = (mat.col(i) - centers.col(j)).squaredNorm();
						if (dc < rr) {
							r2 = dc * an2[j];
							l2 = j;
						}
					}

					if (r2 > d[i]) {
						cluster2[i] = l2;
					}
					else {
						index = 0;
						live[l1] = live[l2] = ncol + i;
						ncp[l1] = ncp[l2] = i + 2;

						double n1 = cluster_count[l1], n2 = cluster_count[l2], alw = n1 - 1, alt = n2 + 1;
						centers.col(l1) = (centers.col(l1) * n1 - mat.col(i)) / alw;
						centers.col(l2) = (centers.col(l2) * n2 + mat.col(i)) / alt;
						++cluster_count[l2];
						--cluster_count[l1];
						cluster[i] = l2;
						cluster2[i] = l1;
						an2[l1] = alw / n1;
						an1[l1] = (alw > 1 ? alw / (alw - 1) : big);
						an1[l2] = alt / n2;
						an2[l2] = alt / (alt + 1);
					}
				}
				if (index == ncol)break;
			}
			if (index == ncol)break;

			for (int i = 0; i < k; ++i) {
				itran[i] = false;
				live[i] -= ncol;
			}

			int icoun = 0, istep = 0;
			while (true) {
				for (int i = 0; i < ncol; ++i) {
					++icoun;
					int l1 = cluster[i];
					if (cluster_count[l1] != 1) {

						if (ncp[l1] >= istep + 2) {
							d[i] = (mat.col(i) - centers.col(l1)).squaredNorm() * an1[l1];
						}

						int l2 = cluster2[i];
						if ((ncp[l1] > istep + 2) || (ncp[l2] > istep + 2)) {
							if ((mat.col(i) - centers.col(l2)).squaredNorm() < d[i] / an2[l2]) {
								icoun = 0;
								index = 0;
								itran[l1] = itran[l2] = true;
								ncp[l1] = ncp[l2] = istep + ncol + 2;

								double n1 = cluster_count[l1], n2 = cluster_count[l2], alw = n1 - 1, alt = n2 + 1;
								centers.col(l1) = (centers.col(l1) * n1 - mat.col(i)) / alw;
								centers.col(l2) = (centers.col(l2) * n2 + mat.col(i)) / alt;
								++cluster_count[l2];
								--cluster_count[l1];
								cluster[i] = l2;
								cluster2[i] = l1;
								an2[l1] = alw / n1;
								an1[l1] = (alw > 1 ? alw / (alw - 1) : big);
								an1[l2] = alt / n2;
								an2[l2] = alt / (alt + 1);
							}
						}
					}
					if (icoun == ncol)break;
					++istep;
					if (istep >= imaxqtr) {
						imaxqtr = -1;
						break;
					}
				}
				if (icoun == ncol)break;
				if (istep >= imaxqtr) {
					imaxqtr = -1;
					break;
				}
			}
			if (imaxqtr < 0) {
				break;
			}
			if (k == 2) {
				break;
			}
			ncp = Eigen::ArrayXi::Ones(k);
		}
		centers.setZero();

		for (int i = 0; i < ncol; ++i) {
			centers.col(cluster[i]) += mat.col(i);
		}

		for (int i = 0; i < k; ++i) {
			if (cluster_count[i] == 0) {
				break;
			}
			centers.col(i) = centers.col(i) / cluster_count[i];
		}

		return std::make_pair(centers, cluster);
	};

	double calculate_total_withinss(const Eigen::MatrixXd& mat, const Eigen::MatrixXd& centers, const Eigen::ArrayXi& cluster) {
		int m = mat.cols();
		Eigen::ArrayXi clu = _Cs unique(cluster);
		double total_within_ss = 0;
		for (int i = 0; i < m; ++i) {
			total_within_ss += (mat.col(i) - centers.col(cluster[i])).squaredNorm();
		}
		return total_within_ss;
	};

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_hartigan_wong_mt(
		const Eigen::MatrixXd& mat, 
		int k, 
		int max_iteration, 
		int random_state, 
		int n_start) 
	{
		if (n_start < 2) {
			return kmeans_hartigan_wong_once(mat, k, max_iteration, random_state);
		}
		else {

			auto [final_centers, final_clusters] = kmeans_hartigan_wong_once(mat, k, max_iteration, random_state);
			double total_within_ss = calculate_total_withinss(mat, final_centers, final_clusters);


		#pragma omp parallel for
			for(int i = 1; i < n_start; ++i){
				auto [centers, cluster] = _Cs kmeans_hartigan_wong_once(mat, k, max_iteration, random_state + i);
				double within_ss = calculate_total_withinss(mat, centers, cluster);

			#pragma omp critical
				{
					if (within_ss < total_within_ss) {
						final_centers = centers;
						final_clusters = cluster;
						total_within_ss = within_ss;
					}
				}
			};
			return std::make_pair(final_centers, final_clusters);
		}
	};

	Eigen::MatrixXd norm2(const Eigen::MatrixXd& mat, bool by_column) {
		Eigen::MatrixXd ret = mat;
		if (by_column) {
			Eigen::ArrayXd col_val = mat.array().colwise().norm();
			int ncol = ret.cols();
			for (int i = 0; i < ncol; ++i) {
				double val = col_val[i];
				if (val != 0) {
					ret.col(i) /= val;
				}
			}
		}
		else {
			Eigen::ArrayXd row_val = mat.array().rowwise().norm();
			int nrow = ret.rows();
			for (int i = 0; i < nrow; ++i) {
				double val = row_val[i];
				if (val != 0) {
					ret.row(i) /= val;
				}
			}
		}
		return ret;
	}

	void normalize_in_place(Eigen::MatrixXd& mat, double scale_factor, bool by_column) {

		if (by_column) {
			Eigen::ArrayXd col_sum = mat.colwise().sum();
			int ncol = mat.cols();

			if (scale_factor <= 0.0) {
				scale_factor = std::abs(_Cs median(col_sum));
			}

			for (int i = 0; i < ncol; ++i) {
				double sum = col_sum[i];
				if (sum != 0) {
					mat.col(i) *= (scale_factor / sum);
				}
			}
		}
		else {
			Eigen::ArrayXd row_sum = mat.rowwise().sum();
			int nrow = mat.rows();

			if (scale_factor <= 0.0) {
				scale_factor = std::abs(_Cs median(row_sum));
			}

			for (int i = 0; i < nrow; ++i) {
				double sum = row_sum[i];
				if (sum != 0) {
					mat.row(i) *= (scale_factor / sum);
				}
			}
		}
	};

	Eigen::MatrixXd normalize(const Eigen::MatrixXd& mat, double scale_factor, bool by_column) {

		Eigen::MatrixXd ret = mat;

		_Cs normalize_in_place(ret, scale_factor, by_column);

		return ret;
	};

	QVector<int> sample_integer(int low, int high, int n, int seed, bool unique) {

		if (n > (high - low + 1) / 10 && unique && n < high - low + 1) {
			QVector<int> index = integer_linspaced(high - low + 1, low, high);
			std::default_random_engine e(seed);
			std::shuffle(index.begin(), index.end(), e);
			return index.sliced(0, n);
		}

		QVector<int> ret(n);
		std::default_random_engine e(seed);
		std::uniform_int_distribution<int> u(low, high);

		for (int i = 0; i < n; ++i) {
			ret[i] = u(e);
		}
		if (!unique || n > high - low + 1) {
			return ret;
		}
		ret = _Cs stable_unique(ret);
		int current_size = ret.size();
		while (current_size < n) {
			int tail_size = n - current_size;
			if (tail_size > 100) {
				QVector<int> tail(tail_size);
				for (int i = 0; i < tail_size; ++i) {
					tail[i] = u(e);
				}
				ret << tail;
				ret = _Cs stable_unique(ret);
				current_size = ret.size();
			}
			else {
				while (tail_size > 0) {
					int tmp = u(e);
					if (!ret.contains(tmp)) {
						ret << tmp;
						--tail_size;
						++current_size;
					}
				}
			}
		}
		return ret;
	};

	std::pair<bool, double> threshold_minimum(const Eigen::ArrayXd& arr) {
		auto [counts, _, bin_centers] = _Cs histogram(arr);

		Eigen::ArrayXd smooth_histogram = counts.cast<double>();

		QList<int> maximum_index;

		auto find_local_maxima_index = [&maximum_index](const Eigen::ArrayXd& hist) {
			maximum_index.clear();
			int direction = 1, length = hist.size();
			for (int i = 0; i < length - 1; ++i) {
				if (direction > 0) {
					if (hist[i + 1] < hist[i]) {
						direction = -1;
						maximum_index.append(i);
					}
				}
				else {
					if (hist[i + 1] > hist[i]) {
						direction = 1;
					}
				}
			}
		};

		int maximum_iteration = 10000;
		int i;
		for (i = 0; i < maximum_iteration; ++i) {
			smooth_histogram = uniform_filter_1d(smooth_histogram, 3);
			find_local_maxima_index(smooth_histogram);
			if (maximum_index.size() < 3)break;
		}

		if (maximum_index.size() != 2 || i == maximum_iteration - 1) {
			return std::make_pair(false, 0);
		}
		int start_index = maximum_index[0];
		int end_index = maximum_index[1];
		int threshold_index;
		smooth_histogram.segment(start_index, end_index - start_index + 1).minCoeff(&threshold_index);
		return std::make_pair(true, bin_centers[start_index + threshold_index]);
	};

	Eigen::ArrayXi uniform_filter_1d(const Eigen::ArrayXi& arr, int size) {
		const Eigen::Index original_size = arr.size();
		const Eigen::Index expanded_size = original_size + size;
		Eigen::ArrayXi expanded(expanded_size);

		const Eigen::Index prefix_length = size / 2;
		const Eigen::Index tail_length = size - prefix_length - 1;

		for (Eigen::Index i = 0; i < prefix_length; ++i) {
			expanded[prefix_length - 1 - i] = arr[i];
		}
		for (Eigen::Index i = 0; i < tail_length; ++i) {
			expanded[original_size + i + prefix_length] = arr[original_size - i - 1];
		}
		expanded.segment(prefix_length, original_size) = arr;

		Eigen::ArrayXi ret(original_size);
		for (Eigen::Index i = 0; i < original_size; ++i) {
			ret[i] = floor(expanded.segment(i, size).mean());
		}
		return ret;
	};

	Eigen::ArrayXd uniform_filter_1d(const Eigen::ArrayXd& arr, int size) {
		const Eigen::Index original_size = arr.size();
		const Eigen::Index expanded_size = original_size + size;
		Eigen::ArrayXd expanded(expanded_size);

		const Eigen::Index prefix_length = size / 2;
		const Eigen::Index tail_length = size - prefix_length - 1;

		for (Eigen::Index i = 0; i < prefix_length; ++i) {
			expanded[prefix_length - 1 - i] = arr[i];
		}
		for (Eigen::Index i = 0; i < tail_length; ++i) {
			expanded[original_size + i + prefix_length] = arr[original_size - i - 1];
		}

		expanded.segment(prefix_length, original_size) = arr;
		Eigen::ArrayXd ret(original_size);
		for (Eigen::Index i = 0; i < original_size; ++i) {
			ret[i] = expanded.segment(i, size).mean();
		}
		return ret;
	};

	QStringList split_lines(const QString& s) {
		return s.split(QRegularExpression("[\r\n]"));
	};

	double linear_percentile(const Eigen::ArrayXd& arr, double p) {

		Eigen::ArrayXd tmp = _Cs sorted(arr);

		if (p <= 0) {
			return tmp[0];
		}

		if (p >= 100) {
			return tmp[tmp.size() - 1];
		}

		const Eigen::Index length = arr.size();

		if (length == 1) {
			return tmp[0];
		}

		if (std::fmod(length - 1, 100 / p) == 0) {
			return tmp[(int)((length - 1) * p / 100)];
		}
		else {
			int low = std::floor((length - 1) * p / 100);
			double res = (length - 1) * p / 100 - low;
			return tmp[low] + res * (tmp[low + 1] - tmp[low]);
		}
	};

	double linear_percentile_no_sort(const Eigen::ArrayXd& arr, double p) {

		if (p <= 0) {
			return arr[0];
		}

		if (p >= 100) {
			return arr[arr.size() - 1];
		}

		const Eigen::Index length = arr.size();

		if (length == 1) {
			return arr[0];
		}

		if (std::fmod(length - 1, 100 / p) == 0) {
			return arr[(int)((length - 1) * p / 100)];
		}
		else {
			int low = std::floor((length - 1) * p / 100);
			double res = (length - 1) * p / 100 - low;
			return arr[low] + res * (arr[low + 1] - arr[low]);
		}
	};

	double iqr(const Eigen::ArrayXd& arr) {
		
		auto arr_sorted = _Cs sorted(arr);

		return (_Cs linear_percentile_no_sort(arr_sorted, 75) - _Cs linear_percentile_no_sort(arr_sorted, 25));
	};

	Eigen::ArrayXd quantile(const Eigen::ArrayXd& arr, const Eigen::ArrayXd& p) {
		const Eigen::Index size = p.size();
		Eigen::ArrayXd ret(size);

		auto arr_sorted = _Cs sorted(arr);

		for (Eigen::Index i = 0; i < size; ++i) {
			ret[i] = _Cs linear_percentile_no_sort(arr_sorted, p[i]);
		}
		return ret;
	};

	double trimean(const Eigen::ArrayXd& arr) {
		return _Cs quantile(arr, Eigen::Array4d{ 25, 50, 50, 75 }).mean();
	};

	Eigen::ArrayXd tricubic_weighting(const Eigen::ArrayXd& x, double loc, double maximum_distance) {

		if (maximum_distance == 0) {
			return Eigen::ArrayXd::Constant(x.size(), 1);
		}

		const int size = x.size();

		Eigen::ArrayXd ret(size);

		for (int i = 0; i < size; ++i) {
			double a = std::abs(x[i] - loc) / maximum_distance;
			a = a * a * a;
			ret[i] = (1 - a) * (1 - a) * (1 - a);
		}

		return ret;
	}

	Eigen::MatrixXd pseudo_inverse(const Eigen::MatrixXd& mat)
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
		double tolerance = 1.e-8;
		int row = mat.rows();
		int col = mat.cols();
		int k = std::min(row, col);

		Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(col, row);
		Eigen::MatrixXd singular_values_inv = svd.singularValues();
		Eigen::MatrixXd singular_values_inv_mat = Eigen::MatrixXd::Zero(col, row);
		for (int i = 0; i < k; ++i) {
			if (singular_values_inv(i) > tolerance) {
				singular_values_inv(i) = 1.0 / singular_values_inv(i);
			}
			else {
				singular_values_inv(i) = 0;
			}
		}
		for (int i = 0; i < k; ++i)
		{
			singular_values_inv_mat(i, i) = singular_values_inv(i);
		}
		ret = (svd.matrixV()) * (singular_values_inv_mat) * (svd.matrixU().transpose());

		return ret;
	}

	Eigen::ArrayXd least_square(const Eigen::ArrayXd& y, const Eigen::ArrayXd& x, const int degree) {

		const int size = x.size();

		Eigen::ArrayXXd tmp = Eigen::ArrayXXd::Ones(size, degree + 1);
		Eigen::MatrixXd left(degree + 1, degree + 1);
		Eigen::VectorXd right(degree + 1);
		
		if (degree < 100) {
			for (int i = 1; i < degree + 1; ++i) {
				for (int j = 0; j < size; ++j) {
					double mul = x[j];
					double val = mul;

					for (int k = 1; k < i; ++k) {
						val *= mul;
					}

					tmp(j, i) = val;
				}
			}

		}
		else {
			for (int i = 1; i < degree + 1; ++i) {
				tmp.col(i) = pow(x, i);
			}
		}
		
		for (int i = 0; i < degree + 1; ++i) {
			for (int j = 0; j < degree + 1; ++j) {
				left(i, j) = (tmp.col(i) * tmp.col(j)).sum();
			}
		}

		for (int i = 0; i < degree + 1; ++i) {
			right(i) = (tmp.col(i) * y).sum();
		}
		return left.inverse() * right;
	};

	Eigen::ArrayXd weighted_least_square(const Eigen::ArrayXd& y, const Eigen::ArrayXd& x, const Eigen::ArrayXd& weight, const int degree) {

		const int size = x.size();

		Eigen::MatrixXd left(degree + 1, degree + 1);
		Eigen::ArrayXXd tmp = Eigen::ArrayXXd::Ones(size, degree + 1);
		Eigen::VectorXd right(degree + 1);

		if (degree < 100) {
			for (int i = 1; i < degree + 1; ++i) {
				for (int j = 0; j < size; ++j) {
					double mul = x[j];
					double val = mul;

					for (int k = 1; k < i; ++k) {
						val *= mul;
					}

					tmp(j, i) = val;
				}
			}

		}
		else {
			for (int i = 1; i < degree + 1; ++i) {
				tmp.col(i) = pow(x, i);
			}
		}

		for (int i = 0; i < degree + 1; ++i) {
			for (int j = 0; j < degree + 1; ++j) {
				left(i, j) = (tmp.col(i) * weight * tmp.col(j)).sum();
			}
		}

		for (int i = 0; i < degree + 1; ++i) {
			right(i) = (tmp.col(i) * weight * y).sum();
		}

		return left.inverse() * right;
	};

	Eigen::ArrayXd loess_mt(const Eigen::ArrayXd& _y, const Eigen::ArrayXd& _x, const int degree, const double span, const int delta) {

		auto index = _Cs order(_x);
		Eigen::ArrayXd y = _y(index);
		Eigen::ArrayXd x = _x(index);

		const int length = x.rows();
		const int width = floor(span * length);
		QVector<int> exact_point;
		for (int i = 0; i < length; i += delta) {
			exact_point << i;
		}
		if (exact_point.last() != length - 1) {
			exact_point << length - 1;
		}
		int exact_count = exact_point.size();

		QVector<double> val(exact_count);

		auto fun = [&x, &y, length, width, degree](int i)->double {
			int left = i;
			int right = i;
			int n = 1;
			double maximum_distance = 0;
			while (true) {
				if (left > 0) {
					if (x[left - 1] == x[left]) {
						left -= 1;
						n += 1;
						continue;
					}
				}
				if (right < length - 1) {
					if (x[right + 1] == x[right]) {
						right += 1;
						n += 1;
						continue;
					}
				}
				if (n >= width)break;
				if (left == 0) {
					right += 1;
					n += 1;
					continue;
				}
				if (right == length - 1) {
					left -= 1;
					n += 1;
					continue;
				}
				if (x[i] - x[left - 1] > x[right + 1] - x[i]) {
					right += 1;
					n += 1;
					continue;
				}
				else {
					left -= 1;
					n += 1;
					continue;
				}
				break;
			}
			maximum_distance = x[right] - x[i] > x[i] - x[left] ? x[right] - x[i] : x[i] - x[left];
			int slice_length = right - left + 1;
			Eigen::ArrayXd local_y = y.segment(left, slice_length);
			Eigen::ArrayXd local_x = x.segment(left, slice_length);
			Eigen::ArrayXd local_weight = tricubic_weighting(local_x, x[i], maximum_distance);
			Eigen::ArrayXd coef = weighted_least_square(local_y, local_x, local_weight, degree);
			if (degree == 1) {
				return coef[1] * x[i] + coef[0];
			}
			else if (degree == 2) {
				return coef[2] * x[i] * x[i] + coef[1] * x[i] + coef[0];
			}
			else {
				return 0;
			}
		};
	#pragma omp parallel for
		for (int i = 0; i < exact_count; ++i) {
			val[i] = fun(exact_point[i]);
		}

		Eigen::ArrayXd fitted(length);
		const qsizetype exact_size = exact_point.size();
		for (qsizetype i = 0; i < exact_size; ++i) {
			fitted[exact_point[i]] = val[i];
		}
		for (qsizetype i = 0; i < exact_size - 1; ++i) {
			qsizetype first = exact_point[i], second = exact_point[i + 1];
			for (qsizetype j = first + 1; j < second; ++j) {
				if (x[j] == x[first]) {
					fitted[j] = fitted[first];
					continue;
				}
				fitted[j] = fitted[first] + (fitted[second] - fitted[first]) * ((x[j] - x[first]) / (x[second] - x[first]));
			}
		}

		Eigen::ArrayXd ret(length);
		for (int i = 0; i < length; ++i) {
			ret[index[i]] = fitted[i];
		}
		return ret;
	};

	QChar detect_delimiter(const QString& line, const QString& line2) {

		int comma_count = line.count(','), tab_count = line.count('\t'), space_count = line.count(' ');

		int comma_count2 = line2.count(','), tab_count2 = line2.count('\t'), space_count2 = line2.count(' ');

		int c = std::min(comma_count, comma_count2);

		int t = std::min(tab_count, tab_count2);

		int s = std::min(space_count, space_count2);

		return c > t ? (c > s ? ',' : ' ') : (t > s ? '\t' : ' ');
	};

	QString capitalize_first(const QString& str) {
		if (str.isEmpty())return QString();
		QChar f = str[0];
		if (f.isLetter()) {
			if (f.isLower()) {
				return f.toUpper() + str.sliced(1, str.size() - 1);
			}
		}
		return str;
	}

	Eigen::ArrayXd cumsum(const Eigen::ArrayXd& vec) {

		const std::size_t size = std::ranges::size(vec);
		Eigen::ArrayXd ret(size);

		auto iter = std::cbegin(vec);
		auto end = std::cend(vec);

		auto to_iter = std::begin(ret);

		double sum = 0;
		for (; iter != end; ++iter, ++to_iter) {

			sum += *iter;
			*to_iter = sum;
		}

		return ret;
	}

	Eigen::MatrixXd row_slice(const Eigen::MatrixXd& mat, const Eigen::ArrayX<bool>& slice) {
		const int ret_row = slice.count();
		const int ncol = mat.cols(), nrow = mat.rows();
		Eigen::MatrixXd res(ret_row, ncol);
		int index = 0;
		for (int i = 0; i < nrow; ++i) {
			if (slice[i]) {
				res.row(index) = mat.row(i);
				++index;
			}
		}
		return res;
	};

	void scale_in_place(Eigen::MatrixXd& mat, bool by_column) {

		if (by_column) {

			Eigen::ArrayXd col_mean = mat.colwise().mean();
			Eigen::ArrayXd col_sdev = mat.colwise().squaredNorm();
			col_sdev = sqrt(col_sdev / (mat.rows() - 1));

			for (auto&& d : col_sdev) {
				if (d == 0.0) {
					d = 1.0;
				}
			}

			mat.array().rowwise() -= col_mean.transpose();
			mat.array().rowwise() /= col_sdev.transpose();
		}
		else {

			Eigen::ArrayXd row_mean = mat.rowwise().mean();
			Eigen::ArrayXd row_sdev = mat.rowwise().squaredNorm();
			row_sdev = sqrt(row_sdev / (mat.cols() - 1));

			for (auto&& d : row_sdev) {
				if (d == 0.0) {
					d = 1.0;
				}
			}

			mat.array().colwise() -= row_mean;
			mat.array().colwise() /= row_sdev;
		}
	};

	Eigen::MatrixXd row_scale_mt(const Eigen::SparseMatrix<double>& mat) {
		const int nrow = mat.rows();
		const int ncol = mat.cols();
		Eigen::SparseMatrix<double, Eigen::RowMajor> tmp(mat);
		Eigen::MatrixXd scaled(nrow, ncol);
		QVector<int> params = integer_linspaced(nrow, 0, nrow - 1);

	#pragma omp parallel for
		for(int i = 0; i < nrow; ++i){
			double row_mean = 0;
			double row_sdev = 0;
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, i); it; ++it)
			{
				row_mean += it.value();
			}
			row_mean /= ncol;
			int n_not_zero = 0;
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(tmp, i); it; ++it)
			{
				n_not_zero += 1;
				row_sdev += std::pow((it.value() - row_mean), 2);
			}
			row_sdev += std::pow(row_mean, 2) * (ncol - n_not_zero);
			row_sdev = sqrt(row_sdev / (ncol - 1));
			if (row_sdev == 0) {
				continue;
			}
			Eigen::VectorXd row = Eigen::VectorXd(tmp.row(i));
			scaled.row(i) = (row.array() - row_mean) / row_sdev;
		}

		return scaled;
	};

	Eigen::SparseMatrix<double> normalize(const Eigen::SparseMatrix<double>& mat, double scale_factor) {
		Eigen::SparseMatrix<double> res = mat;
		Eigen::ArrayXd col_sum = _Cs col_sum(mat);
		for (int k = 0; k < res.outerSize(); ++k) {
			double column_sum = col_sum[k];
			if (column_sum != 0) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(res, k); it; ++it) {
					it.valueRef() = it.value() / column_sum * scale_factor;
				}
			}
		}
		return res;
	}

	Eigen::SparseMatrix<double> row_normalize(const Eigen::SparseMatrix<double>& mat, double scale_factor) {
		Eigen::SparseMatrix<double> res = mat;
		Eigen::ArrayXd row_sum = _Cs row_sum(mat);

		for (int k = 0; k < res.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(res, k); it; ++it) {
				double sum = row_sum[it.row()];
				if (sum != 0) {
					it.valueRef() = it.value() / sum * scale_factor;
				}
			}
		}
		return res;
	};

	Eigen::SparseMatrix<double> row_normalize2(const Eigen::SparseMatrix<double>& mat, double scale_factor) {
		Eigen::SparseMatrix<double> res = mat;
		Eigen::ArrayXd row_sum = _Cs row_sum_abs(mat);

		for (int k = 0; k < res.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(res, k); it; ++it) {
				double sum = row_sum[it.row()];
				if (sum != 0) {
					it.valueRef() = it.value() / sum * scale_factor;
				}
			}
		}
		return res;
	};

	void quiver_scale(const Eigen::MatrixXd& x, Eigen::MatrixXd& v, int factor) {

		auto [x_min, x_max] = std::ranges::minmax(x.col(0));
		auto [y_min, y_max] = std::ranges::minmax(x.col(1));

		double xv_max = v.col(0).array().abs().maxCoeff(), yv_max = v.col(1).array().abs().maxCoeff();
		double x_scale = (x_max - x_min) / xv_max / factor, y_scale = (y_max - y_min) / yv_max / factor;

		v.array() *= x_scale > y_scale ? y_scale : x_scale;
	};

	Eigen::SparseMatrix<double> normalize(const Eigen::SparseMatrix<int>& mat, double scale_factor) {
		Eigen::SparseMatrix<double> res = mat.cast<double>();
		Eigen::ArrayXd col_sum = _Cs col_sum(mat).cast<double>();
		for (int k = 0; k < res.outerSize(); ++k) {
			if (col_sum[k] != 0) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(res, k); it; ++it) {
					it.valueRef() = it.value() / col_sum[k] * scale_factor;
				}
			}
		}
		return res;
	}

	void normalize_in_place(Eigen::SparseMatrix<double>& mat, double scale_factor) {
		Eigen::ArrayXd col_sum = _Cs col_sum(mat);

		if (scale_factor < 0.0) {
			scale_factor = std::abs(_Cs median(col_sum));
		}

		for (int k = 0; k < mat.outerSize(); ++k) {
			if (col_sum[k] != 0) {

				double m = scale_factor / col_sum[k];

				for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
					it.valueRef() = it.value() * m;
				}
			}
		}
	}

}
