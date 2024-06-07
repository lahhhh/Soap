#pragma once

#ifndef _Cs
#define _Cs ::custom::
#endif // !_Cs

#include "Identifier.h"

#include <string>
#include <random>
#include <type_traits>
#include <vector>
#include <ranges>
#include <unordered_set>

#include <QFontMetrics>
#include <QColor>
#include <QStringList>
#include <igraph.h>

#include "StringSplitView.h"

namespace custom {

	Eigen::ArrayXd col_var_mt(const Eigen::MatrixXd& mat);
	Eigen::ArrayXd row_var_mt(const Eigen::MatrixXd& mat);

	double sd(const Eigen::MatrixXd& mat);

	StringSplitView split_view(const std::string& s, char delimiter);

	std::vector<std::string> split(std::string_view str, char delimiter);

	Eigen::ArrayXi cluster_louvain_igraph(const Eigen::SparseMatrix<double>& adj, double cutoff);

	void sort_by_first(Eigen::ArrayXd& first, Eigen::MatrixXd& second, bool decrease = false);

	Eigen::SparseMatrix<double> set_from_triplets_mean(const std::vector<Eigen::Triplet<double>>& triplets, int nrow, int ncol);

	Eigen::ArrayXi find_nearest(const Eigen::MatrixXd& query, const Eigen::MatrixXd& data, bool by_col = true);

	Eigen::MatrixXd cov(const Eigen::MatrixXd& mat);
	Eigen::MatrixXd cov2cor(const Eigen::MatrixXd& mat);

	void igraph_vector_copy(igraph_vector_t& to, const igraph_vector_t& from);

	void igraph_vector_int_copy(igraph_vector_int_t& to, const igraph_vector_int_t& from);

	Eigen::ArrayXd mean_compress(const Eigen::ArrayXd& arr, int time);

	double sign(double val);

	QString standardize_windows_file_name(const QString& file_name);

	QString complementary_strand(const QString& sequence);

	QString standardize_chromosome_name(const QString& chr);
	std::string standardize_chromosome_name(const std::string& chr);

	QStringList paste(const QString& p1, const QStringList& p2, const QString& concat = "");
	QStringList paste(const QStringList& p1, const QString& p2, const QString& concat = "");
	QStringList paste(const QStringList& p1, const QStringList& p2, const QString& concat = "");

	QVector<int> seq_n(int start, int size, int incr = 1);

	void rowwise_divide_in_place(Eigen::SparseMatrix<double>& dividend, const Eigen::ArrayXd& divisor);

	// modified from Velocyto
	Eigen::MatrixXi balanced_knn_mt(
		const Eigen::MatrixXd& distance,
		int n_neighbor,
		int maxl,
		bool var_by_column = true,
		const QString& method = "cor"
	);

	Eigen::ArrayXd cor(const Eigen::MatrixXd& mat, const Eigen::ArrayXd& arr);
	Eigen::ArrayXd cor_mt(const Eigen::MatrixXd& mat, const Eigen::ArrayXd& arr);
	Eigen::MatrixXd cor(const Eigen::MatrixXd& data, bool var_by_column = true);

	Eigen::MatrixXd euclidean_distance_mt(const Eigen::MatrixXd& data, bool var_by_column = true);
	Eigen::MatrixXd correlation_pearson_distance_mt(const Eigen::MatrixXd& data, bool var_by_column = true);
	Eigen::MatrixXd distance(const QVector<double>& vec);

	// modified from seurat
	Eigen::SparseMatrix<double> create_shared_nearest_neighbors_matrix(const Eigen::MatrixXi& knn);

	Eigen::SparseMatrix<int> create_matrix_from_knn_index(const Eigen::MatrixXi& knn);

	QVector<double> empirical_cumulative_distribution(const QVector<double>& vec);

	// [low, high]
	QVector<int> sample_integer(int low, int high, int n, int seed = 1997, bool unique = true);

	// 0 : string, 1 : numeric, 2 : integer
	int detect_type(const QStringList& vec);

	bool is_integer(const char* str);
	bool is_integer(const QString& str);
	bool is_integer(const QStringList& vec);

	bool is_numeric(const char* str);
	bool is_numeric(const QString& str);
	bool is_numeric(const QStringList& vec);

	int line_length(const char* str);

	std::vector<int> generate_random_permutation(std::size_t n_element, std::default_random_engine& re);

	QVector<int> integer_divide(const QVector<int>& dividend, const int divisor);

	// only positive integer
	int atoi_specialized(const char* src);
	int atoi_specialized(const char** src);

	Eigen::ArrayXXd row_compressed(const Eigen::ArrayXXd& mat, int rate);

	Eigen::ArrayXXd column_compressed(const Eigen::ArrayXXd& mat, int rate);

	Eigen::ArrayXXd column_compressed(const Eigen::ArrayXXd& mat, const QVector<int>& segment_width);

	Eigen::ArrayXXd row_compressed(const Eigen::ArrayXXd& mat, const QVector<int>& segment_width);

	QColor random_color();

	QStringList multi_split(const QString& src, const QList<QChar>& splits);

	QString string_next_line(const QFontMetrics& font, const QString& text, int label_size, bool is_English = true);

	Eigen::SparseMatrix<double> eliminate_less_than(const Eigen::SparseMatrix<double>& mat, double threshold);
	Eigen::SparseMatrix<double> eliminate_zero(const Eigen::SparseMatrix<double>& mat);

	QVector<double> linspaced(int length, double begin, double end);

	QVector<int> integer_linspaced(int length, int begin, int end);

	QChar detect_delimiter(const QString& line, const QString& line2);

	QStringList digest(const QString& line, const QChar delimiter = ',');

	QStringList digest(const QString& line, const QString& delimiter);

	void digest(
		const std::string& line,
		const char delimiter,
		std::size_t start,
		std::size_t end,
		std::vector<std::pair<std::size_t, std::size_t>>& line_loc);

	void digest(
		const QString& line, 
		const QChar delimiter, 
		qsizetype start, 
		qsizetype end, 
		std::vector<std::pair<qsizetype, qsizetype>>& line_loc);

	QString merge_to_string(const QStringList& vec, const QChar delimiter, bool quote = true);
	QString merge_to_string(const QStringList& vec, const QString& delimiter);

	QStringList to_upper(const QStringList& list);

	double signal_to_noise(const Eigen::ArrayXd& X, const Eigen::ArrayXd& Y);

	double max(const Eigen::SparseMatrix<double>& mat);

	QStringList make_unique(const QStringList& list);

	Eigen::ArrayXd set_upper_bound(const Eigen::ArrayXd& vec, double bound);

	Eigen::ArrayXd adjust_p_value(const Eigen::ArrayXd& pval, const QString& method);

	QVector<double> adjust_p_value(const QVector<double>& pval, const QString& method);

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_lloyd(
		const Eigen::MatrixXd& mat, 
		int k, 
		bool remove_empty_cluster = false,
		int max_iteration = 20,
		int random_state = 1997
	);

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_lloyd(
		const Eigen::MatrixXd& mat, 
		Eigen::MatrixXd centers, 
		bool remove_empty_cluster = false,
		int max_iteration = 20);

	// k >= 2
	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_hartigan_wong_mt(const Eigen::MatrixXd& mat, int k, int max_iteration = 20, int random_state = 1997, int n_start = 1);

	std::pair<Eigen::MatrixXd, Eigen::ArrayXi> kmeans_hartigan_wong_once(const Eigen::MatrixXd& mat, int k, int max_iteration, int random_state);

	double calculate_total_withinss(const Eigen::MatrixXd& mat, const Eigen::MatrixXd& centers, const Eigen::ArrayXi& cluster);

	std::pair<bool, double> threshold_minimum(const Eigen::ArrayXd& arr);

	Eigen::ArrayXi uniform_filter_1d(const Eigen::ArrayXi& arr, int size); // mode : reflect

	Eigen::ArrayXd uniform_filter_1d(const Eigen::ArrayXd& arr, int size); // mode : reflect

	QStringList split_lines(const QString& s);

	double linear_percentile(const Eigen::ArrayXd& arr, double p);

	double linear_percentile_no_sort(const Eigen::ArrayXd& arr, double p);

	Eigen::ArrayXd quantile(const Eigen::ArrayXd& arr, const Eigen::ArrayXd& p);
	double iqr(const Eigen::ArrayXd& arr);

	double trimean(const Eigen::ArrayXd& arr);

	Eigen::ArrayXd tricubic_weighting(const Eigen::ArrayXd& x, double loc, double maximum_distance);	

	Eigen::MatrixXd pseudo_inverse(const Eigen::MatrixXd& mat);

	Eigen::ArrayXd least_square(const Eigen::ArrayXd& y, const Eigen::ArrayXd& x, const int degree);

	Eigen::ArrayXd weighted_least_square(const Eigen::ArrayXd& y, const Eigen::ArrayXd& x, const Eigen::ArrayXd& weight, const int degree);

	QString capitalize_first(const QString& str);

	Eigen::ArrayXd loess_mt(const Eigen::ArrayXd& _y, const Eigen::ArrayXd& _x, const int degree, const double span, const int delta);

	// note : scale by row
	Eigen::MatrixXd row_scale_mt(const Eigen::SparseMatrix<double>& mat);

	void scale_in_place(Eigen::MatrixXd& mat, bool by_column = true);

	void quiver_scale(const Eigen::MatrixXd& x, Eigen::MatrixXd& v, int factor = 40);

	Eigen::SparseMatrix<double> row_normalize(const Eigen::SparseMatrix<double>& mat, double scale_factor);

	Eigen::SparseMatrix<double> row_normalize2(const Eigen::SparseMatrix<double>& mat, double scale_factor);

	Eigen::SparseMatrix<double> normalize(const Eigen::SparseMatrix<double>& mat, double scale_factor);

	Eigen::SparseMatrix<double> normalize(const Eigen::SparseMatrix<int>& mat, double scale_factor);

	void normalize_in_place(Eigen::SparseMatrix<double>& mat, double scale_factor = -1.0);

	// rowsum / colsum == scale_factor
	Eigen::MatrixXd normalize(const Eigen::MatrixXd& mat, double scale_factor = -1.0, bool by_column = true);

	void normalize_in_place(Eigen::MatrixXd& mat, double scale_factor = -1.0, bool by_column = true);

	// arma::normalise(mat, 2, )
	Eigen::MatrixXd norm2(const Eigen::MatrixXd& mat, bool by_column = true);

	Eigen::ArrayXd cumsum(const Eigen::ArrayXd& vec);
}

