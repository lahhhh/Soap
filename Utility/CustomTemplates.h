#pragma once

#ifndef _Cs
#define _Cs ::custom::
#endif // !_Cs

#include "Identifier.h"

#include "CustomTypeTraits.h"
#include "CustomFunctions.h"

#include "annoylib.h"
#include "kissrandom.h"
#include "mman.h"

namespace custom {

	template<typename Dst, typename Src>
		requires _Cs not_container<Src>&& _Cs not_same<Dst, QString>
	SOAP_INLINE Dst cast(Src&& src) {
		return static_cast<Dst>(std::forward<Src>(src));
	}


	//--------- cast to qstring -----------//


	template<typename Dst, typename Src>
		requires std::is_arithmetic_v<std::decay_t<Src>>&& std::same_as<Dst, QString>
	SOAP_INLINE Dst cast(Src src) {
		return QString::number(src);
	}


	//--------- cast to qstring -----------//

	//--------- cast from qstring -----------//


	template<typename Dst>
	SOAP_INLINE Dst cast(const QString& src) {
		return static_cast<Dst>(src);
	}

	template<>
	SOAP_INLINE double cast<double>(const QString& src) {
		return src.toDouble();
	}

	template<>
	SOAP_INLINE int cast<int>(const QString& src) {
		return src.toInt();
	}


	//--------- cast from qstring -----------//


	template<
		template<typename...>typename Outer,
		typename SrcCon,
		typename Ele = _Cs get_element_raw_type<SrcCon>,
		typename ReCon = Outer<Ele>
	>
		requires _Cs construct_from_size<ReCon>
	ReCon cast(SrcCon&& src) {

		auto size = std::size(src);

		ReCon ret(size);

		auto iter = std::cbegin(src);
		auto end = std::cend(src);

		auto to_iter = std::begin(ret);

		for (; iter != end; ++iter, ++to_iter) {
			*to_iter = *iter;
		}

		return ret;
	}

	template<
		template<typename, int, int...>typename Outer,
		typename SrcCon,
		typename Ele = _Cs get_element_raw_type<SrcCon>,
		int Int, int... Ints
	>
	auto cast(SrcCon&& src) {

		using ReCon = Outer<Ele, Int, Ints...>;

		auto size = std::size(src);

		ReCon ret(size);

		auto iter = std::cbegin(src);
		auto end = std::cend(src);

		auto to_iter = std::begin(ret);

		for (; iter != end; ++iter, ++to_iter) {
			*to_iter = *iter;
		}

		return ret;
	}

	template<
		typename Ele,
		typename SrcCon,
		typename DstCon = change_element_type<std::decay_t<SrcCon>, Ele>
	>
		requires _Cs is_container<DstCon>
	&& _Cs is_container<SrcCon>
		&& _Cs construct_from_size<DstCon>
		DstCon cast(SrcCon&& src) {
		auto size = std::size(src);
		DstCon ret(size);

		auto iter = std::cbegin(src);
		auto end = std::cend(src);
		auto to_iter = std::begin(ret);

		for (; iter != end; ++iter, ++to_iter) {
			*to_iter = _Cs cast<Ele>(*iter);
		}

		return ret;
	}

	template<typename T>
	bool is_same(const T& con) {
		if (std::size(con) == 0) {
			return false;
		}

		auto first = std::cbegin(con);
		auto iter = std::next(first, 1);
		auto end = std::cend(con);

		for (; iter != end; ++iter) {
			if (*iter != *first) {
				return false;
			}
		}

		return true;
	}

	template<typename Con, typename Ele = get_element_raw_type<Con>>
	int find_most_frequent_element(const Con& con) {
		std::unordered_map<Ele, int> count_map;

		auto iter = std::cbegin(con);
		auto end = std::cend(con);

		for (; iter != end; ++iter) {
			++count_map[*iter];
		}

		Ele max_element = *std::cbegin(con);
		int max_count = 0;

		for (auto&& [ele, count] : count_map) {
			if (count > max_count) {
				max_count = count;
				max_element = ele;
			}
		}

		return max_element;
	}

	template<typename Mat>
	Mat cbind(const Mat& mat1, const Mat& mat2) {

		int ncol1 = mat1.cols();
		int ncol2 = mat2.cols();
		int nrow = mat1.rows();

		Mat m(nrow, ncol1 + ncol2);

		m.block(0, 0, nrow, ncol1) = mat1;
		m.block(0, ncol1, nrow, ncol2) = mat2;

		return m;
	}

	template<typename Mat>
	Mat cbind(const Mat& mat1, const Mat& mat2, const Mat& mat3) {

		int ncol1 = mat1.cols();
		int ncol2 = mat2.cols();
		int ncol3 = mat3.cols();
		int nrow = mat1.rows();

		Mat m(nrow, ncol1 + ncol2 + ncol3);

		m.block(0, 0, nrow, ncol1) = mat1;
		m.block(0, ncol1, nrow, ncol2) = mat2;
		m.block(0, ncol2 + ncol1, nrow, ncol3) = mat3;

		return m;
	}

	template<typename Mat>
	Mat rbind(const Mat& mat1, const Mat& mat2) {

		int nrow1 = mat1.rows();
		int nrow2 = mat2.rows();
		int ncol = mat1.cols();

		Mat m(nrow1 + nrow2, ncol);

		m.block(0, 0, nrow1, ncol) = mat1;
		m.block(nrow1, 0, nrow2, ncol) = mat2;

		return m;
	}

	template <typename Con>	requires is_container<Con>&& construct_from_size<Con>
	Con concatenated(
		const Con& container1,
		const Con& container2
	) {
		auto size1 = std::size(container1);
		auto size2 = std::size(container2);

		auto size = size1 + size2;

		Con ret(size);
		auto iter = std::begin(ret);

		for (auto&& ele1 : container1) {
			*(iter++) = ele1;
		}

		for (auto&& ele2 : container2) {
			*(iter++) = ele2;
		}

		return ret;
	}

	template <typename Con, typename Ele>
		requires is_container<Con>
	SOAP_INLINE void extend_minmax(const Con& container, Ele& min_val, Ele& max_val) {

		Ele min{ min_val }, max{ max_val };

		for (auto&& ele : container) {

			if (ele < min) {
				min = ele;
			}

			if (max < ele) {
				max = ele;
			}
		}

		min_val = min;
		max_val = max;
	}

	template <typename Con, typename Ele>
		requires is_container<Con>
	SOAP_INLINE bool all(const Con& container, Ele val) {

		for (auto&& ele : container) {
			if (ele != val) {
				return false;
			}
		}

		return true;
	}

	template <typename Con>
		requires is_specific_container<Con, double>
	SOAP_INLINE void remove_na(Con& container, double replace = 0.0) {

		for (auto&& d : container) {
			if (std::isnan(d)) {
				d = replace;
			}
		}
	}

	template <
		typename Con,
		typename Ele,
		typename Project = std::identity
	>
		requires is_container<Con>
	bool any(const Con& container, Ele val, Project pj = {}) {

		for (auto&& ele : container) {
			if (pj(ele) == val) {
				return true;
			}
		}

		return false;
	}

	template <typename Key, typename Value>
	SOAP_INLINE QList<Key> keys(const std::map<Key, Value>& map) {

		QList<Key> ret;

		std::ranges::for_each(map, [&ret](auto&& p) {ret << p.first; });

		return ret;
	}

	template <typename Key, typename Value>
	SOAP_INLINE QList<Value> values(const std::map<Key, Value>& map) {

		QList<Value> ret;

		std::ranges::for_each(map, [&ret](auto&& p) {ret << p.second; });

		return ret;
	}

	template <
		typename Con,
		typename ReCon = std::decay_t<Con>
	>
		requires _Cs is_arithmetic_container<Con>
	ReCon abs(Con&& container) {

		ReCon ret(std::forward<Con>(container));

		for (auto& val : ret) {
			val = std::abs(val);
		}

		return ret;
	}

	template <typename S, typename _S = std::decay_t<S>>
		requires _Cs is_slice_container<_S>
	_S flip(S&& slice) {

		_S flipped(std::forward<S>(slice));

		for (auto& bool_val : flipped) {
			bool_val = !bool_val;
		}

		return flipped;
	}

	template <
		typename Con,
		typename Ele = get_element_raw_type<Con>
	>
		requires less_comparable<Ele>
	&& is_random_access_container<Con>
		Eigen::ArrayXi order(const Con& ct, bool decrease = false) {

		const std::size_t size = std::ranges::size(ct);

		Eigen::ArrayXi ind = Eigen::ArrayXi::LinSpaced(size, 0, size - 1);

		if (decrease) {
			auto rule = [&ct](int i, int j)->bool {
				return ct[j] < ct[i];
			};
			std::sort(ind.begin(), ind.end(), rule);
		}
		else {
			auto rule = [&ct](int i, int j)->bool {
				return ct[i] < ct[j];
			};
			std::sort(ind.begin(), ind.end(), rule);
		}

		return ind;
	};

	template<typename Val>
	auto SOAP_INLINE constexpr min(const Val& val)
	{
		return val;
	}

	template<typename Val, typename Val2, typename... Vals>
	auto SOAP_INLINE constexpr min(const Val& val1, const Val2& val2, const Vals&... vals)
	{
		return val1 < val2 ? _Cs min(val1, vals...) : _Cs min(val2, vals...);
	}

	template<typename Con>
	int argmin(const Con& ct) {

		auto iter = std::cbegin(ct);
		auto end = std::cend(ct);

		int ind{ 0 }, min_loc{ -1 };

		if (iter == end) {
			return min_loc;
		}

		auto min = iter;
		min_loc = 0;

		while (++iter != end) {

			++ind;

			if (*iter < *min) {
				min = iter;
				min_loc = ind;
			}
		}

		return min_loc;
	}

	template<typename Val>
	auto SOAP_INLINE constexpr max(const Val& val)
	{
		return val;
	}

	template<typename Val, typename Val2, typename... Vals>
	auto SOAP_INLINE constexpr max(const Val& val1, const Val2& val2, const Vals&... vals)
	{
		return val1 < val2 ? _Cs max(val2, vals...) : _Cs max(val1, vals...);
	}

	template<typename Con>
	int argmax(const Con& ct) {

		auto iter = std::cbegin(ct);
		auto end = std::cend(ct);

		int ind{ 0 }, max_loc{ -1 };

		if (iter == end) {
			return max_loc;
		}

		auto max = iter;
		max_loc = 0;

		while (++iter != end) {

			++ind;

			if (*iter > *max) {
				max = iter;
				max_loc = ind;
			}
		}

		return max_loc;
	}

	template <
		typename Con,
		typename S
	>
		requires  is_container<Con>
	&& is_slice_container<S>
		&& construct_from_size<Con>
		[[nodiscard]] auto sliced(const Con& ct, const S& slice) {

		const auto count = std::ranges::count(slice, true);
		const auto size = std::ranges::size(slice);

		Con res(count);

		auto slice_iter = std::cbegin(slice);
		auto slice_end = std::cend(slice);

		auto from_iter = std::cbegin(ct);
		auto to_iter = std::begin(res);

		for (; slice_iter != slice_end; ++slice_iter, ++from_iter) {
			if (*slice_iter) {
				*to_iter = *from_iter;
				++to_iter;
			}
		}

		return res;
	};

	template <
		typename Con,
		typename S
	>
		requires  is_container<Con>
	&& is_slice_container<S>
		&& construct_from_size<Con>
		[[nodiscard]] auto sliced_ptr(const Con& ct, const S& slice) {

		const auto count = std::ranges::count(slice, true);
		const auto size = std::ranges::size(slice);

		Con* res = new Con(count);

		auto slice_iter = std::cbegin(slice);
		auto slice_end = std::cend(slice);

		auto from_iter = std::cbegin(ct);
		auto to_iter = std::begin(*res);

		for (; slice_iter != slice_end; ++slice_iter, ++from_iter) {
			if (*slice_iter) {
				*to_iter = *from_iter;
				++to_iter;
			}
		}

		return res;
	};

	template<typename S>
		requires is_slice_container<S>
	QVector<int> which(const S& slice) {

		QVector<int> index;

		int i{ 0 };

		for (auto&& ele : slice) {

			if (ele) {
				index << i;
			}

			++i;
		}

		return index;
	}

	template <typename Ele>
	QVector<int> get_column_index(const Eigen::SparseMatrix<Ele>& sp_mat, int ind) {
		const int start = sp_mat.m_outerIndex[ind], end = sp_mat.m_outerIndex[ind + 1], length = end - start;

		QVector<int> index(length);
		for (int i = start; i < end; ++i) {
			index[i - start] = sp_mat.m_data.m_indices[i];
		}

		return index;
	}

	template <typename Con, typename Ele = get_element_raw_type<Con> >
		requires  is_container<Con>
	SOAP_INLINE auto sum(const Con& ct)
	{
		return std::accumulate(std::cbegin(ct), std::cend(ct), static_cast<Ele>(0));
	}

	template <typename Con, typename OrderCon>
		requires is_random_access_container<Con>
	&& is_order_container<OrderCon>
		&& construct_from_size<Con>
		[[nodiscard]] Con reordered(const Con& ct, const OrderCon& order) {

		const std::size_t size = std::ranges::size(order);

		Con res(size);

		auto iter = std::cbegin(order);
		auto end = std::cend(order);

		auto to_iter = std::begin(res);

		for (; iter != end; ++iter, ++to_iter) {
			*to_iter = ct[*iter];
		}

		return res;
	};

	template <typename Con, typename OrderCon>
		requires is_random_access_container<Con>
	&& is_order_container<OrderCon>
		&& construct_from_size<Con>
		[[nodiscard]] Con* reordered_ptr(const Con& ct, const OrderCon& order) {

		const std::size_t size = std::ranges::size(order);

		Con* res = new Con(size);

		auto iter = std::cbegin(order);
		auto end = std::cend(order);

		auto to_iter = std::begin(*res);

		for (; iter != end; ++iter, ++to_iter) {
			*to_iter = ct[*iter];
		}

		return res;
	};

	template <
		typename Ele,
		typename Function
	>
	SOAP_INLINE void elementwise_in_place(
		Eigen::SparseMatrix<Ele>& sp_mat,
		const Function& fun) {

		for (int k = 0; k < sp_mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(sp_mat, k); it; ++it) {
				it.valueRef() = fun(it);
			}
		}
	}

	template <
		typename Ele,
		typename Function
	>
	[[nodiscard]] SOAP_INLINE Eigen::SparseMatrix<Ele>
		elementwise(
			const Eigen::SparseMatrix<Ele>& sp_mat,
			const Function& fun)
	{

		std::vector<Eigen::Triplet<Ele>> triplets;
		const auto nnz = sp_mat.nonZeros();
		triplets.reserve(nnz);

		for (int k = 0; k < sp_mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(sp_mat, k); it; ++it) {
				triplets.emplace_back(it.row(), k, fun(it));
			}
		}

		const auto n_row = sp_mat.rows(), n_col = sp_mat.cols();

		Eigen::SparseMatrix<Ele> ret;
		ret.resize(n_row, n_col);
		ret.setFromTriplets(triplets.cbegin(), triplets.cend());

		return ret;
	}

	// sp_mat should be square matrix
	template <typename Ele>
	[[nodiscard]] Eigen::SparseMatrix<Ele>
		set_diagonaled(const Eigen::SparseMatrix<Ele>& sp_mat, Ele val) {

		const int n_col = sp_mat.cols();

		std::vector<Eigen::Triplet<Ele> > triplets;
		auto nnz = sp_mat.nonZeros();
		triplets.reserve(nnz + n_col);

		for (int k = 0; k < sp_mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(sp_mat, k); it; ++it) {
				if (it.row() != k) {
					triplets.emplace_back(it.row(), k, it.value());
				}
			}
		}

		for (int i = 0; i < n_col; ++i) {
			triplets.emplace_back(i, i, val);
		}

		Eigen::SparseMatrix<Ele> ret;
		ret.resize(n_col, n_col);
		ret.setFromTriplets(triplets.cbegin(), triplets.cend());

		return ret;
	}

	// sp_mat should be square matrix
	template <typename Ele>
	[[nodiscard]] Eigen::SparseMatrix<Ele>
		set_diagonaled(const Eigen::SparseMatrix<Ele>& sp_mat, const Eigen::ArrayX<Ele>& arr) {

		const int n_col = sp_mat.cols();

		std::vector<Eigen::Triplet<Ele> > triplets;
		auto nnz = sp_mat.nonZeros();
		triplets.reserve(nnz + n_col);

		for (int k = 0; k < sp_mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(sp_mat, k); it; ++it) {
				if (it.row() != k) {
					triplets.emplace_back(it.row(), k, it.value());
				}
			}
		}

		for (int i = 0; i < n_col; ++i) {
			triplets.emplace_back(i, i, arr[i]);
		}

		Eigen::SparseMatrix<Ele> ret;
		ret.resize(n_col, n_col);
		ret.setFromTriplets(triplets.cbegin(), triplets.cend());

		return ret;
	}

	// sp_mat should be square matrix
	template <typename Ele>
	void set_diagonal(Eigen::SparseMatrix<Ele>& sp_mat, Ele val) {

		const int n_col = sp_mat.cols();

		std::vector<Eigen::Triplet<Ele> > triplets;
		auto nnz = sp_mat.nonZeros();
		triplets.reserve(nnz + n_col);

		for (int k = 0; k < sp_mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(sp_mat, k); it; ++it) {
				if (it.row() != k) {
					triplets.emplace_back(it.row(), k, it.value());
				}
			}
		}

		for (int i = 0; i < n_col; ++i) {
			triplets.emplace_back(i, i, val);
		}

		sp_mat.setFromTriplets(triplets.cbegin(), triplets.cend());
	}

	// sp_mat should be square matrix
	template <typename Ele>
	void set_diagonal(Eigen::SparseMatrix<Ele>& sp_mat, const Eigen::ArrayX<Ele>& arr) {

		const int n_col = sp_mat.cols();

		std::vector<Eigen::Triplet<Ele> > triplets;
		auto nnz = sp_mat.nonZeros();
		triplets.reserve(nnz + n_col);

		for (int k = 0; k < sp_mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(sp_mat, k); it; ++it) {
				if (it.row() != k) {
					triplets.emplace_back(it.row(), k, it.value());
				}
			}
		}

		for (int i = 0; i < n_col; ++i) {
			triplets.emplace_back(i, i, arr[i]);
		}

		sp_mat.setFromTriplets(triplets.cbegin(), triplets.cend());
	}

	template <typename Ele>
	Eigen::MatrixX<Ele> sign(const Eigen::MatrixX<Ele>& con) {

		Eigen::MatrixX<Ele> ret(con);

		const int nrow = ret.rows(), ncol = ret.cols();

		for (int j = 0; j < ncol; ++j) {
			for (int i = 0; i < nrow; ++i) {
				Ele val = ret(i, j);
				ret(i, j) = val > static_cast<Ele>(0) ?
					static_cast<Ele>(1) : (val < static_cast<Ele>(0) ?
						static_cast<Ele>(-1) : static_cast<Ele>(0));
			}
		}

		return ret;
	}

	template <typename Con>
		requires is_random_access_container<Con>
	std::tuple<Eigen::ArrayXd, bool, double> rank_average(const Con& ct) {

		const auto size = std::ranges::size(ct);

		Eigen::ArrayXi index = _Cs order(ct);
		Eigen::ArrayXd ret(size);

		int count{ 1 };
		bool tie{ false };
		double n_tie{ 0.0 };
		double current_value{ ct[index[0]] };

		for (std::size_t i = 1; i < size; ++i) {

			double this_value = ct[index[i]];

			if (current_value != this_value) {
				if (count > 1) {
					tie = true;
					n_tie += (count * count * count - count);
					double rank_mean = i - (count + 1) / 2.0;
					for (int j = 0; j < count; ++j) {
						ret[index[i - j - 1]] = rank_mean;
					}
				}
				else {
					ret[index[i - 1]] = i - 1;
				}
				count = 1;
			}
			else {
				++count;
			}
			current_value = this_value;
		}

		if (count > 1) {
			tie = true;
			n_tie += (count * count * count - count);
			double rank_mean = size - (count + 1) / 2.0;
			for (int j = 0; j < count; ++j) {
				ret[index[size - j - 1]] = rank_mean;
			}
		}
		else {
			ret[index[size - 1]] = size - 1;
		}

		return std::make_tuple(ret + 1, tie, n_tie);
	};

	// note : the return index include itself at 1st column
	template <typename KnnMethod = Euclidean, bool ReDistance = false>
	auto get_knn_mt(const Eigen::MatrixXd& mat, const int n_neighbors, const int n_trees = 50) {
		const int search_k = 2 * n_neighbors * n_trees;
		const int ncol = mat.cols(), nrow = mat.rows();

		AnnoyIndex<int, double, KnnMethod, Kiss64Random, AnnoyIndexSingleThreadedBuildPolicy> ann(ncol);

		for (Eigen::Index i = 0; i < nrow; ++i) {
			std::vector<double> trans(mat.row(i).cbegin(), mat.row(i).cend());
			ann.add_item(i, trans.data());// row.data() will cause error
		}

		ann.build(n_trees);

		Eigen::MatrixXi index(nrow, n_neighbors);
		Eigen::MatrixXd distances;

		if constexpr (ReDistance) {
			distances.resize(nrow, n_neighbors);
		}

	#pragma omp parallel for
		for (int i = 0; i < nrow; ++i) {
			std::vector<double> quest(mat.row(i).cbegin(), mat.row(i).cend());
			std::vector<int> neighbors;
			std::vector<double> distance;

			ann.get_nns_by_vector(quest.data(), n_neighbors, search_k, &neighbors, &distance);
			for (int j = 0; j < n_neighbors; ++j) {
				index(i, j) = neighbors[j];
				if constexpr (ReDistance) {
					distances(i, j) = distance[j];
				}
			}
		}

		ann.unload();

		if constexpr (ReDistance) {
			return std::make_pair(index, distances);
		}
		else {
			return index;
		}
	};


	template <typename KnnMethod, bool ReDistance = false>
	auto get_knn_mt(const Eigen::MatrixXd& base, const Eigen::MatrixXd& query, const int n_neighbors, const int n_trees = 50) {
		const int search_k = 2 * n_neighbors * n_trees;
		const int ncol = base.cols(), nrow = base.rows();

		AnnoyIndex<int, double, KnnMethod, Kiss64Random, AnnoyIndexSingleThreadedBuildPolicy> ann(ncol);

		for (Eigen::Index i = 0; i < nrow; ++i) {
			std::vector<double> trans(base.row(i).cbegin(), base.row(i).cend());
			ann.add_item(i, trans.data());// row.data() will cause error
		}

		ann.build(n_trees);

		const int nrow_query = query.rows();
		Eigen::MatrixXi index(nrow_query, n_neighbors);
		Eigen::MatrixXd distances;

		if constexpr (ReDistance) {
			distances.resize(nrow_query, n_neighbors);
		}

	#pragma omp parallel for
		for (int i = 0; i < nrow_query; ++i) {
			std::vector<double> quest(query.row(i).cbegin(), query.row(i).cend());
			std::vector<int> neighbors;
			std::vector<double> distance;

			ann.get_nns_by_vector(quest.data(), n_neighbors, search_k, &neighbors, &distance);
			for (int j = 0; j < n_neighbors; ++j) {
				index(i, j) = neighbors[j];
				if constexpr (ReDistance) {
					distances(i, j) = distance[j];
				}
			}
		}

		ann.unload();

		if constexpr (ReDistance) {
			return std::make_pair(index, distances);
		}
		else {
			return index;
		}
	};

	template <typename Con, typename Value, typename Op>
		requires std::is_arithmetic_v<Value>&& std::same_as<get_element_raw_type<Con>, Value>&& construct_from_size<Con>
	[[nodiscard]] Con calculate(const Con& ct, Value val, Op op) {

		const std::size_t size = std::ranges::size(ct);

		Con ret(size);

		auto from_iter = std::cbegin(ct);
		auto from_end = std::cend(ct);

		auto to_iter = std::begin(ret);

		for (; from_iter != from_end; ++from_iter, ++to_iter) {
			*to_iter = op(*from_iter, val);
		}

		return ret;
	}

	template <typename Con, typename Op>
		requires is_arithmetic_container<Con> && construct_from_size<Con>
		[[nodiscard]] Con calculate(const Con& ct1, const Con& ct2, Op op) {

		const std::size_t size = std::ranges::size(ct1);

		Con ret(size);

		auto iter1 = std::cbegin(ct1);
		auto end = std::cend(ct1);

		auto iter2 = std::cbegin(ct2);

		auto ret_iter = std::begin(ret);

		for (; iter1 != end; ++iter1, ++iter2, ++ret_iter) {
			*ret_iter = op(*iter1, *iter2);
		}

		return ret;
	}

	template <typename Con, typename Value>
		requires std::is_arithmetic_v<Value>&& std::same_as<get_element_raw_type<Con>, Value>&& construct_from_size<Con>
	[[nodiscard]] Con add(const Con& ct, Value val) {

		return _Cs calculate(ct, val, std::plus{});
	}

	template <typename Con>
		requires is_arithmetic_container<Con>&& construct_from_size<Con>
	[[nodiscard]] Con add(const Con& ct1, const Con& ct2) {

		return _Cs calculate(ct1, ct2, std::plus{});
	}

	template <typename Con, typename Value>
		requires std::is_arithmetic_v<Value>&& std::same_as<get_element_raw_type<Con>, Value>&& construct_from_size<Con>
	[[nodiscard]] Con minus(const Con& ct, Value val) {

		return _Cs calculate(ct, val, std::minus{});
	}

	template <typename Con>
		requires is_arithmetic_container<Con>&& construct_from_size<Con>
	[[nodiscard]] Con minus(const Con& ct1, const Con& ct2) {

		return _Cs calculate(ct1, ct2, std::minus{});
	}

	template <typename Con, typename Value>
		requires std::is_arithmetic_v<Value>&& std::same_as<get_element_raw_type<Con>, Value>&& construct_from_size<Con>
	[[nodiscard]] Con multiply(const Con& ct, Value val) {

		return _Cs calculate(ct, val, std::multiplies{});
	}

	template <typename Con>
		requires is_arithmetic_container<Con>&& construct_from_size<Con>
	[[nodiscard]] Con multiply(const Con& ct1, const Con& ct2) {

		return _Cs calculate(ct1, ct2, std::multiplies{});
	}

	template <
		typename Con1,
		typename Con2,
		bool Check = true,
		typename Ele1 = get_element_raw_type<Con1>,
		typename Ele2 = get_element_raw_type<Con2>,
		typename ReEle = decltype(std::declval<Ele1>() / std::declval<Ele2>()),
		typename ReCon = change_element_type<Con1, ReEle>
	>
		requires is_arithmetic_container<Con1>
	&& is_arithmetic_container<Con2>
		&& construct_from_size<ReCon>
		[[nodiscard]]
	auto divide(
		const Con1& dividend,
		const Con2& divisor,
		ReEle zero_division_result = static_cast<ReEle>(0)
	) {

		const std::size_t size = std::ranges::size(dividend);

		ReCon ret(size);

		auto dividend_iter = std::cbegin(dividend);
		auto dividend_end = std::cend(dividend);

		auto divisor_iter = std::cbegin(divisor);

		auto ret_iter = std::begin(ret);


		for (; dividend_iter != dividend_end; ++dividend_iter, ++divisor_iter, ++ret_iter) {

			if constexpr (Check) {
				if (*divisor_iter == static_cast<Ele2>(0)) { // check
					*ret_iter = zero_division_result;
				}
				else {
					*ret_iter = static_cast<ReEle>(*dividend_iter) / static_cast<ReEle>(*divisor_iter);
				}
			}
			else {
				*ret_iter = static_cast<ReEle>(*dividend_iter) / static_cast<ReEle>(*divisor_iter);
			}
		}

		return ret;
	}

	template <
		typename ReEle,
		typename Con1,
		typename Con2,
		bool Check = true,
		typename Ele1 = get_element_raw_type<Con1>,
		typename Ele2 = get_element_raw_type<Con2>,
		typename ReCon = change_element_type<Con1, ReEle>
	>
		requires is_arithmetic_container<Con1>
	&& is_arithmetic_container<Con2>
		&& construct_from_size<ReCon>
		[[nodiscard]]
	ReCon  partial_divide(
		const Con1& dividend,
		const Con2& divisor,
		ReEle zero_division_result = static_cast<ReEle>(0)) {

		const std::size_t size = std::ranges::size(dividend);

		ReCon ret(size);

		auto dividend_iter = std::cbegin(dividend);
		auto dividend_end = std::cend(dividend);

		auto divisor_iter = std::cbegin(divisor);

		auto ret_iter = std::begin(ret);

		for (; dividend_iter != dividend_end; ++dividend_iter, ++divisor_iter, ++ret_iter) {

			if constexpr (Check) {
				if (*divisor_iter == static_cast<Ele2>(0)) {
					*ret_iter = zero_division_result;
				}
				else {
					*ret_iter = static_cast<ReEle>(*dividend_iter) / static_cast<ReEle>(*divisor_iter);
				}
			}
			else {
				*ret_iter = static_cast<ReEle>(*dividend_iter) / static_cast<ReEle>(*divisor_iter);
			}
		}

		return ret;
	}

	template <
		typename Con,
		typename Value,
		typename Ele = get_element_raw_type<Con>,
		typename ReEle = decltype(std::declval<Ele>() / std::declval<Value>()),
		typename ReCon = change_element_type<Con, ReEle>
	>
		requires is_arithmetic_container<Con>
	&& std::is_arithmetic_v<Value>
		&& construct_from_size<ReCon>
		ReCon
		divide(
			const Con& dividend,
			const Value divisor)
	{
		const std::size_t size = std::ranges::size(dividend);

		ReCon ret(size);

		auto from_iter = std::cbegin(dividend);
		auto from_end = std::cend(dividend);

		auto to_iter = std::begin(ret);

		for (; from_iter != from_end; ++from_iter, ++to_iter) {
			*to_iter = *from_iter / divisor;
		}

		return ret;
	}

	template <
		typename ReEle,
		typename Con,
		typename Value,
		typename ReCon = change_element_type<Con, ReEle>
	>
		requires  is_arithmetic_container<Con>
	&& std::is_arithmetic_v<Value>
		&& construct_from_size<ReCon>
		ReCon
		partial_divide(
			const Con& dividend,
			const Value divisor)
	{
		const std::size_t size = std::ranges::size(dividend);

		ReCon ret(size);

		auto from_iter = std::cbegin(dividend);
		auto from_end = std::cend(dividend);

		auto to_iter = std::begin(ret);

		for (; from_iter != from_end; ++from_iter, ++to_iter) {
			*to_iter = static_cast<ReEle>(*from_iter) / static_cast<ReEle>(divisor);
		}

		return ret;
	}

	template <typename Ele, typename Order>
		requires std::is_arithmetic_v<Ele>&& is_order_container<Order>
	[[nodiscard]]
	Eigen::ArrayX<Ele> column_reorder_and_row_sum(const Eigen::SparseMatrix<Ele>& mat, const Order& order) {

		const int ncol = mat.cols(), nrow = mat.rows();

		Eigen::ArrayX<Ele> ret = Eigen::ArrayX<Ele>::Zero(nrow);

		for (auto&& ele : order) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, ele); it; ++it) {
				ret[it.row()] += it.value();
			}
		}

		return ret;
	}

	template <typename Ele, typename S>
		requires std::is_arithmetic_v<Ele>&& is_slice_container<S>
	[[nodiscard]]
	Eigen::ArrayX<Ele> column_slice_and_row_sum(const Eigen::SparseMatrix<Ele>& mat, const S& slice) {

		const int ncol = mat.cols(), nrow = mat.rows();

		Eigen::ArrayX<Ele> ret = Eigen::ArrayX<Ele>::Zero(nrow);

		for (int i = 0; i < ncol; ++i) {
			if (slice[i]) {
				for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, i); it; ++it) {
					ret[it.row()] += it.value();
				}
			}
		}

		return ret;
	}

	template <typename Ele, typename S>
		requires std::is_arithmetic_v<Ele>&& is_slice_container<S>
	[[nodiscard]]
	Eigen::ArrayXd column_slice_and_row_mean(const Eigen::SparseMatrix<Ele>& mat, const S& slice) {

		auto count = (double)std::ranges::count(slice, true);

		auto sum = column_slice_and_row_sum(mat, slice);

		Eigen::ArrayXd s = sum.cast<double>();

		if (count > 0.0) {
			s /= count;
		}

		return s;
	}

	template <typename Ele, typename Order>
		requires std::is_arithmetic_v<Ele>&& is_order_container<Order>
	[[nodiscard]]
	Eigen::ArrayX<Ele> row_reorder_and_column_sum(const Eigen::SparseMatrix<Ele>& mat, const Order& order) {

		const int ncol = mat.cols();

		Eigen::ArrayX<Ele> ret = Eigen::ArrayX<Ele>::Zero(ncol);

		for (std::size_t i = 0; i < ncol; ++i) {
			ret[i] = _Cs sum(_Cs reordered(Eigen::ArrayX<Ele>(mat.col(i)), order));
		}

		return ret;
	}

	template <typename Ele, typename Order>
		requires std::is_arithmetic_v<Ele>&& is_order_container<Order>
	Eigen::ArrayXd row_reorder_and_column_mean(const Eigen::SparseMatrix<Ele>& mat, const Order& order) {

		const std::size_t order_size = std::ranges::size(order);

		Eigen::ArrayXd ret = _Cs row_reorder_and_column_sum(mat, order).cast<double>();

		ret /= (double)order_size;

		return ret;
	}

	template <typename Ele>
		requires std::is_arithmetic_v<Ele>
	QVector<int> get_row_index(const Eigen::SparseMatrix<Ele>& sp_mat, int ind) {

		Eigen::ArrayX<Ele> row = sp_mat.row(ind);

		return _Cs which(row != static_cast<Ele>(0));
	}

	template <typename Ele>
		requires std::is_arithmetic_v<Ele>
	Eigen::ArrayX<Ele> get_row_data(const Eigen::SparseMatrix<Ele>& sp_mat, int ind) {

		Eigen::ArrayX<Ele> row = sp_mat.row(ind);

		return _Cs sliced(row, row != static_cast<Ele>(0));
	}

	template<typename Con, typename Ele = get_element_raw_type<Con>>
		requires is_random_access_container<Con>
	Ele median(const Con& con) {

		auto order = _Cs order(con);

		const std::size_t size = std::ranges::size(con);

		if (size % 2 == 1) {
			const std::size_t median_loc = size / 2;
			return con[order[median_loc]];
		}
		else {
			const std::size_t median_loc1 = size / 2;
			const std::size_t median_loc2 = median_loc1 - 1;

			return (con[order[median_loc1]] + con[order[median_loc2]]) / static_cast<Ele>(2);
		}
	}

	template <typename Con>
		requires is_arithmetic_container<Con>
	SOAP_INLINE double mean(const Con& ct) {

		const std::size_t size = std::ranges::size(ct);

		return static_cast<double>(_Cs sum(ct)) / size;
	}

	template <typename Con>
		requires is_arithmetic_container<Con>
	double var(const Con& ct) {

		std::size_t size = std::ranges::size(ct);

		if (size <= 1) {
			return std::nan(0);
		}

		double mean = _Cs mean(ct), variance{ 0.0 };

		for (auto&& ele : ct) {

			variance += (ele - mean) * (ele - mean);
		}

		return variance / (size - 1);
	}

	template <typename Con>
		requires is_arithmetic_container<Con>
	SOAP_INLINE double sd(const Con& ct) {

		return std::sqrt(_Cs var(ct));
	}

	template <typename Con, typename Function, typename Re = get_return_type<Function, get_element_arg_type<Con>>>
		requires is_container<Con> && std::same_as<Re, bool>
		[[nodiscard]] QVector<int> match(Con&& ct, Function&& fun)
	{
		QVector<int> index;

		int i{ 0 };

		for (auto&& ele : ct) {
			if (fun(ele)) {
				index << i;
			}

			++i;
		}

		return index;
	}

	template <typename Con, typename Function, typename Re = sapply_return_type<Function, get_element_arg_type<Con>>>
		Re sapply(Con&& ct, Function&& fun)
	{
		constexpr bool has_return_value = !std::same_as<Re, void>;

		if constexpr (has_return_value) {

			const std::size_t size = std::ranges::size(ct);

			Re ret(size);

			auto iter = std::begin(ct);
			auto end = std::end(ct);

			auto ret_iter = std::begin(ret);

			for (; iter != end; ++iter, ++ret_iter) {

				*ret_iter = fun(*iter);
			}

			return ret;
		}
		else {

			for (auto&& ele : ct) {

				fun(ele);
			}
		}
	}

	template < typename Outer, typename Inner = get_element_raw_type<Outer>, typename Ele = get_element_raw_type<Inner>>
		requires is_container<Outer>&& is_container<Inner>
	QVector<Ele> unroll(const Outer& ct)
	{
		const std::size_t size = _Cs sum(_Cs sapply(ct, [](auto&& ele) { return std::size(ele); }));

		QVector<Ele> ret;
		ret.reserve(size);

		for (auto&& ele : ct) {

			for (auto&& inner_ele : ele) {
				ret.append(inner_ele);
			}
		}

		return ret;
	}

	template <typename T>
		requires std::is_arithmetic_v<T>
	std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXd>
		histogram(const Eigen::ArrayX<T>& orig, int bins = 64, double zero_space = 0.01) {
		Eigen::ArrayXd arr = orig.template cast<double>();

		std::size_t original_size = orig.size();

		double amin = arr.minCoeff(), amax = arr.maxCoeff();

		if (amin == amax) {

			Eigen::ArrayXi counts = Eigen::ArrayXi::Zero(bins);
			counts[bins / 2] = original_size;

			Eigen::ArrayXd edges = Eigen::ArrayXd::LinSpaced(bins + 1, amin - bins * zero_space / 2, amin + bins * zero_space / 2);

			Eigen::ArrayXd centers = Eigen::ArrayXd::LinSpaced(bins, amin - (bins - 1) * zero_space / 2, amin + (bins - 1) * zero_space / 2);

			return std::make_tuple(counts, edges, centers);
		}

		Eigen::ArrayXd edges = Eigen::ArrayXd::LinSpaced(bins + 1, amin, amax);

		double d = (amax - amin) / bins;

		Eigen::ArrayXi counts = Eigen::ArrayXi::Zero(bins + 1);

		for (std::size_t i = 0; i < original_size; ++i) {
			++counts[floor((arr[i] - amin) / d)];
		}

		counts[bins - 1] += counts[bins];

		Eigen::ArrayXd centers = Eigen::ArrayXd::LinSpaced(bins, amin + d / 2, amax - d / 2);

		return std::make_tuple(counts.segment(0, bins), edges, centers);
	};

	template <typename Con>
	SOAP_INLINE Con reversed(const Con& ct) {

		return Con(std::ranges::crbegin(ct), std::ranges::crend(ct));
	}

	template <typename Con, typename Ele = get_element_raw_type<Con>>
		requires less_comparable<Ele>
	SOAP_INLINE void sort(Con& ct, bool decrease = false) {

		if (decrease) {
			std::ranges::sort(ct, std::greater<Ele>());
		}
		else {
			std::ranges::sort(ct, std::less<Ele>());
		}
	};

	template <typename Con, typename ReCon = std::decay_t<Con>, typename Ele = get_element_raw_type<ReCon>>
		requires less_comparable<Ele>&& is_random_access_container<ReCon>
	[[nodiscard]] SOAP_INLINE auto sorted(Con&& ct, bool decrease = false) {

		ReCon res(std::forward<Con>(ct));

		_Cs sort(res, decrease);

		return res;
	};

	template <typename Con, typename Con2, typename Ele = get_element_raw_type<Con>>
		requires less_comparable<Ele>&& is_random_access_container<Con>
	&& is_random_access_container<Con2>
		void sort_by_first(Con& first, Con2& second, bool decrease = false) {

		Eigen::ArrayXi ind = _Cs order(first, decrease);

		first = _Cs reordered(first, ind);

		second = _Cs reordered(second, ind);
	};

	template <typename Con, typename Ele = get_element_raw_type<Con>>
		requires less_comparable<Ele>&& is_random_access_container<Con>
	Eigen::ArrayXi stable_order(const Con& ct, bool decrease = false) {

		const std::size_t size = std::ranges::size(ct);

		Eigen::ArrayXi ind = Eigen::ArrayXi::LinSpaced(size, 0, size - 1);

		if (decrease) {
			auto rule = [&ct](int i, int j)->bool {
				return ct[j] < ct[i];
			};
			std::stable_sort(ind.begin(), ind.end(), rule);
		}
		else {
			auto rule = [&ct](int i, int j)->bool {
				return ct[i] < ct[j];
			};
			std::stable_sort(ind.begin(), ind.end(), rule);
		}

		return ind;
	};

	// length should be less than size
	template <typename Con>
	Con sample(const Con& ct, std::size_t length, unsigned int random_state) {

		Con res(length);

		std::default_random_engine re(random_state);

		std::ranges::sample(ct, std::begin(res), length, re);

		return res;
	};

	template <typename Con>
		requires  is_slice_container<Con>
	long long get_first_index(const Con& container) {

		auto iter = std::cbegin(container);
		auto end = std::cend(container);

		for (long long i = 0; iter != end; ++iter, ++i) {
			if (*iter) {
				return i;
			}
		}

		return -1;
	}

	// y.index_of(x)
	template <typename Con1, typename Con2>
		requires same_element<Con1, Con2>&& is_random_access_container<Con1>&& is_random_access_container<Con2>
	QVector<int> index_of(const Con1& x, const Con2& y) {

		Eigen::ArrayXi order_x = _Cs stable_order(x), order_y = _Cs stable_order(y);

		const std::size_t size_x = std::ranges::size(x), size_y = std::ranges::size(y);

		QVector<int> ret(size_x, -1);

		std::size_t rank_y = 0, rank_x = 0;

		while (rank_y < size_y && rank_x < size_x) {

			int loc_x = order_x[rank_x], loc_y = order_y[rank_y];

			if (x[loc_x] == y[loc_y]) {
				ret[loc_x] = loc_y;
				++rank_x;
			}
			else if (x[loc_x] < y[loc_y]) {
				++rank_x;
			}
			else {
				++rank_y;
			}
		}

		return ret;
	}

	template <typename Con1, typename Con2>
		requires same_element<Con1, Con2>&& is_random_access_container<Con1>&& is_random_access_container<Con2>
	QVector<int> last_index_of(const Con1& x, const Con2& y) {

		Eigen::ArrayXi order_x = _Cs stable_order(x), order_y = _Cs stable_order(y);

		const std::size_t size_x = std::ranges::size(x), size_y = std::ranges::size(y);

		QVector<int> ret(size_x, -1);

		long long rank_y = size_y - 1, rank_x = size_x - 1;

		while (rank_y >= 0 && rank_x >= 0) {

			int loc_x = order_x[rank_x], loc_y = order_y[rank_y];

			if (x[loc_x] == y[loc_y]) {
				ret[loc_x] = loc_y;
				--rank_x;
			}
			else if (x[loc_x] > y[loc_y]) {
				--rank_x;
			}
			else {
				--rank_y;
			}
		}

		return ret;
	}

	// should be unique elements
	template <typename Con1, typename Con2, typename _Con1 = std::decay_t<Con1>, typename _Con2 = std::decay_t<Con2>>
		requires same_element<_Con1, _Con2>
	std::size_t intersect_length(Con1&& ct1, Con2&& ct2) {

		_Con1 x(std::forward<Con1>(ct1));
		_Con2 y(std::forward<Con2>(ct2));

		std::ranges::sort(x);
		std::ranges::sort(y);

		auto x_last = std::unique(std::begin(x), std::end(x));
		auto x_first = std::begin(x);

		auto y_last = std::unique(std::begin(y), std::end(y));
		auto y_first = std::begin(y);

		std::size_t intersect = 0;

		while (x_first != x_last && y_first != y_last) {
			if (*x_first < *y_first) {
				++x_first;
			}
			else if (*y_first < *x_first) {
				++y_first;
			}
			else {
				++x_first;
				++y_first;
				++intersect;
			}
		}

		return intersect;
	}

	// number of elements only in ct1
	template <typename Con1, typename Con2, typename _Con1 = std::decay_t<Con1>, typename _Con2 = std::decay_t<Con2>>
		requires same_element<_Con1, _Con2>
	std::size_t set_difference_length(Con1&& ct1, Con2&& ct2) {

		_Con1 x(std::forward<Con1>(ct1));
		_Con2 y(std::forward<Con2>(ct2));

		std::ranges::sort(x);
		std::ranges::sort(y);

		auto x_last = std::unique(std::begin(x), std::end(x));
		auto x_first = std::begin(x);

		auto y_last = std::unique(std::begin(y), std::end(y));
		auto y_first = std::begin(y);


		std::size_t length = 0;

		while (x_first != x_last && y_first != y_last) {
			if (*x_first < *y_first) {
				++x_first;
				++length;
			}
			else if (*y_first < *x_first) {
				++y_first;
			}
			else {
				++x_first;
				++y_first;
			}
		}

		if (x_first != x_last) {
			length += x_last - x_first;
		}

		return length;
	}

	// set_difference(x, y) - remove y from x
	template <typename Con1, typename Con2, typename ReCon = std::decay_t<Con1>>
		requires std::same_as<std::decay_t<Con1>, std::decay_t<Con2>>&& is_container<ReCon>
	ReCon set_difference(Con1&& ct1, Con2&& ct2) {

		ReCon x(std::forward<Con1>(ct1));
		ReCon y(std::forward<Con2>(ct2));

		std::ranges::sort(x);
		std::ranges::sort(y);

		auto x_last = std::unique(std::begin(x), std::end(x));
		auto x_first = std::begin(x);

		auto y_last = std::unique(std::begin(y), std::end(y));
		auto y_first = std::begin(y);

		ReCon diff;
		auto bi = std::back_inserter(diff);

		while (x_first != x_last && y_first != y_last) {
			if (*x_first < *y_first) {
				bi = *x_first;
				++x_first;
			}
			else if (*y_first < *x_first) {
				++y_first;
			}
			else {
				++x_first;
				++y_first;
			}
		}

		while (x_first != x_last) {
			bi = *x_first;
			++x_first;
		}

		return diff;
	}

	template <typename Con1, typename Con2>
		requires same_element<Con1, Con2>&& is_random_access_container<Con1>&& is_random_access_container<Con2>
	QVector<int> SOAP_INLINE valid_index_of(const Con1& query, const Con2& data) {

		QVector<int> ret = _Cs index_of(query, data);

		ret.removeAll(-1);

		return ret;
	};

	template <typename To, typename From, typename S>
		requires same_element<To, From>&& is_slice_container<S>
	void assign(To& to_container, const From& from_container, const S& slice) {

		auto from_iter = std::cbegin(from_container);
		auto from_end = std::cend(from_container);

		auto to_iter = std::begin(to_container);

		auto slice_iter = std::cbegin(slice);

		for (; from_iter != from_end; ++from_iter, ++to_iter, ++slice_iter) {
			if (*slice_iter) {
				*to_iter = *from_iter;
			}
		}
	}

	template <typename Con, typename Ele, typename S, typename ToEle = get_element_raw_type<Con>>
		requires  is_container<Con>&& is_slice_container<S>&& std::same_as<ToEle, Ele>
	void assign(Con& to_container, const Ele& from, const S& slice) {

		auto to_iter = std::begin(to_container);
		auto to_end = std::end(to_container);

		auto slice_iter = std::cbegin(slice);

		for (; to_iter != to_end; ++to_iter, ++slice_iter) {
			if (*slice_iter) {
				*to_iter = from;
			}
		}
	}

	template <typename Con, typename Ele, typename Order, typename ToEle = get_element_raw_type<Con>>
		requires is_container<Con>
	&& is_order_container<Order>
		&& std::same_as<ToEle, Ele>
		void assign(Con& to_container, const Ele& from, const Order& order) {

		for (auto&& ele : order) {

			to_container[ele] = from;
		}
	}

	template <typename To, typename From, typename Order>
		requires same_element<To, From>&& is_order_container<Order>
	void assign(To& to_container, const From& from_container, const Order& order) {

		auto from_iter = std::cbegin(from_container);
		auto from_end = std::cend(from_container);

		auto to_iter = std::begin(to_container);

		auto order_iter = std::cbegin(order);

		for (; from_iter != from_end; ++from_iter, ++to_iter, ++order_iter) {
			to_container[*order_iter] = *from_iter;
		}
	}

	template <typename To, typename From, typename Order, typename Re = std::decay_t<To>>
		requires same_element<To, From>&& is_order_container<Order>
	[[nodiscard]] Re assigned(To&& to_container, const From& from_container, const Order& order) {

		Re ret(std::forward<To>(to_container));

		auto from_iter = std::cbegin(from_container);
		auto from_end = std::cend(from_container);

		auto order_iter = std::cbegin(order);

		for (; from_iter != from_end; ++from_iter, ++order_iter) {

			ret[*order_iter] = *from_iter;
		}

		return ret;
	}

	template <typename To, typename From, typename ToOrder, typename FromOrder, typename Re = std::decay_t<To>>
		requires same_element<To, From>&& is_order_container<ToOrder>&& is_order_container<FromOrder>
	[[nodiscard]] Re assigned(
		To&& to_container, const ToOrder& to_order, const From& from_container, const FromOrder& from_order) {

		Re ret(std::forward<To>(to_container));

		auto from_order_iter = std::cbegin(from_order);
		auto from_order_end = std::cend(from_order);

		auto to_order_iter = std::cbegin(to_order);

		for (; from_order_iter != from_order_end; ++from_order_iter, ++to_order_iter) {

			ret[*to_order_iter] = from_container[*from_order_iter];
		}

		return ret;
	}

	// return sorted elements
	template <typename Con, typename ReCon = std::decay_t<Con>, typename Ele = get_element_raw_type<ReCon>>
		requires less_comparable<Ele>&& is_container<Con>&& construct_from_size<ReCon>
	SOAP_INLINE ReCon unique(Con&& ct) {

		ReCon tmp(std::forward<Con>(ct));

		std::ranges::sort(tmp);

		auto unique_iter = std::unique(std::begin(tmp), std::end(tmp));

		const auto unique_number = unique_iter - std::begin(tmp);

		ReCon ret(unique_number);

		auto to_iter = std::begin(ret);
		auto from_iter = std::begin(tmp);

		for (; from_iter != unique_iter; ++to_iter, ++from_iter) {
			*to_iter = *from_iter;
		}

		return ret;
	}

	template <typename Con, typename _Con = std::decay_t<Con>, typename Ele = get_element_raw_type<_Con>>
		requires less_comparable<Ele>&& is_container<_Con>
	SOAP_INLINE auto unique_element_number(Con&& ct) {

		_Con tmp(std::forward<Con>(ct));

		std::ranges::sort(tmp);

		auto unique_iter = std::unique(std::begin(tmp), std::end(tmp));

		return unique_iter - std::begin(tmp);
	}

	template <typename Con, typename Ele = get_element_raw_type<Con> >
	Con stable_unique(const Con& ct) {

		Con ret;

		auto bi = std::back_inserter(ret);

		std::unordered_set<Ele> value_set;

		for (auto&& ele : ct) {
			if (!value_set.contains(ele)) {
				value_set.insert(ele);
				bi = ele;
			}
		}

		return ret;
	}

	template <typename T>
	Eigen::ArrayX<T> stable_unique(const Eigen::ArrayX<T>& ct) {

		return _Cs cast<Eigen::ArrayX>(_Cs stable_unique(_Cs cast<QVector>(ct)));
	}

	template <typename Con, typename Ele, typename Comp = std::less<Ele>>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	Eigen::ArrayX<bool> compare(const Con& ct, const Ele& val, Comp cp = {}) {
		const std::size_t size = std::ranges::size(ct);

		Eigen::ArrayX<bool> ret(size);

		auto iter = std::cbegin(ct);
		auto end = std::cend(ct);

		auto predict = std::begin(ret);

		for (; iter != end; ++iter, ++predict) {
			*predict = cp(*iter, val);
		}

		return ret;
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE Eigen::ArrayX<bool> equal(const Con& ct, const Ele& val) {

		return _Cs compare(ct, val, std::equal_to<Ele>{});
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE Eigen::ArrayX<bool> not_equal(const Con& ct, const Ele& val) {

		return _Cs compare(ct, val, std::not_equal_to<Ele>{});
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE Eigen::ArrayX<bool> greater_than(const Con& ct, const Ele& val) {

		return _Cs compare(ct, val, std::greater<Ele>{});
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE Eigen::ArrayX<bool> greater_equal(const Con& ct, const Ele& val) {

		return _Cs compare(ct, val, std::greater_equal<Ele>{});
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE Eigen::ArrayX<bool> less_than(const Con& ct, const Ele& val) {

		return _Cs compare(ct, val, std::less<Ele>{});
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE Eigen::ArrayX<bool> less_equal(const Con& ct, const Ele& val) {

		return _Cs compare(ct, val, std::less_equal<Ele>{});
	}

	template <typename Con, typename Ele>
		requires std::same_as<Ele, get_element_raw_type<Con>>
	SOAP_INLINE QVector<int> match(const Con& ct, const Ele& val) {

		return match(ct, [&val](auto&& ele) {return val == ele; });
	}

	template<typename Con, typename Func>
		requires random_accessible<Con>
	Eigen::ArrayX<bool> test(const Con& ct, const Func& fun) {

		const std::size_t size = std::ranges::size(ct);

		Eigen::ArrayX<bool> ret(size);

		for (std::size_t i = 0; i < size; ++i) {
			ret[i] = fun(ct[i]);
		}

		return ret;
	};

	template<typename Re, typename Condition, typename Con1, typename Con2>
		requires is_bool_container<Condition>&& is_random_access_container<Condition>
	&& same_element<Con1, Con2>&& construct_from_size<Re>&& is_random_access_container<Con1>
		&& is_random_access_container<Con2>&& same_element<Con1, Re>
		Re ifelse(const Condition& condition, const Con1& yes_value, const Con2& no_value) {

		const std::size_t size = std::ranges::size(condition);

		Re ret(size);

		for (std::size_t i = 0; i < size; ++i) {
			ret[i] = condition[i] ? yes_value[i] : no_value[i];
		}

		return ret;
	}

	template<typename Re, typename Condition, typename Con, typename Value>
		requires is_bool_container<Condition>
	&& is_random_access_container<Condition>
		&& std::same_as<get_element_raw_type<Con>, Value>
		&& construct_from_size<Re>
		&& is_random_access_container<Con>
		Re ifelse(const Condition& condition, const Con& yes_value, const Value& no_value) {

		const std::size_t size = std::ranges::size(condition);

		Re ret(size);

		for (std::size_t i = 0; i < size; ++i) {
			ret[i] = condition[i] ? yes_value[i] : no_value;
		}

		return ret;
	}

	template<typename Re, typename Condition, typename Con, typename Value>
		requires is_bool_container<Condition>&& is_random_access_container<Condition>
	&& std::same_as<get_element_raw_type<Con>, Value>&& construct_from_size<Re>
		&& is_random_access_container<Con>
		Re ifelse(const Condition& condition, const Value& yes_value, const Con& no_value) {

		const std::size_t size = std::ranges::size(condition);

		Re ret(size);

		for (std::size_t i = 0; i < size; ++i) {
			ret[i] = condition[i] ? yes_value : no_value[i];
		}

		return ret;
	}

	template<typename Re, typename Value, typename Condition>
		requires is_bool_container<Condition>&& is_random_access_container<Condition>
	&& construct_from_size<Re>&& std::same_as<get_element_raw_type<Re>, Value>
		Re ifelse(const Condition& condition, const Value& yes_value, const Value& no_value) {
		const std::size_t size = std::ranges::size(condition);

		Re ret(size);
		for (std::size_t i = 0; i < size; ++i) {
			ret[i] = condition[i] ? yes_value : no_value;
		}

		return ret;
	}

	template <typename T>
	Eigen::ArrayX<T> col_sum(const Eigen::SparseMatrix<T>& mat) {
		int ncol = mat.cols();
		Eigen::ArrayX<T> res(ncol);
		for (int k = 0; k < mat.cols(); ++k) {
			T colsum = 0;
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
				colsum += it.value();
			}
			res[k] = colsum;
		}
		return res;
	}

	template <typename T>
	Eigen::ArrayX<T> col_sum_mt(const Eigen::SparseMatrix<T>& mat) {
		int ncol = mat.cols();
		Eigen::ArrayX<T> res(ncol);
	#pragma omp parallel for
		for (int i = 0; i < ncol; ++i) {
			T colsum = 0;
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it) {
				colsum += it.value();
			}
			res[i] = colsum;
		}
		return res;
	}

	template <typename T>
	Eigen::ArrayX<T> row_max(const Eigen::SparseMatrix<T>& mat) {
		Eigen::ArrayX<T> res = Eigen::ArrayX<T>::Zero(mat.rows());
		for (int k = 0; k < mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
				if (res[it.row()] < it.value()) {
					res[it.row()] = it.value();
				}
			}
		}
		return res;
	}

	template <typename T>
	Eigen::ArrayX<T> col_max(const Eigen::SparseMatrix<T>& mat) {
		Eigen::ArrayX<T> res = Eigen::ArrayX<T>::Zero(mat.cols());
		for (int k = 0; k < mat.outerSize(); ++k) {
			T max_val{ 0 };
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
				if (max_val < it.value()) {
					max_val = it.value();
				}
			}
			res[k] = max_val;
		}
		return res;
	}

	template <typename T>
	Eigen::ArrayX<T> row_sum(const Eigen::SparseMatrix<T>& mat) {
		Eigen::ArrayX<T> res = Eigen::ArrayX<T>::Zero(mat.rows());
		for (int k = 0; k < mat.outerSize(); ++k) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
				res[it.row()] += it.value();
			}
		}
		return res;
	}

	template <typename T>
		requires std::is_arithmetic_v<T>
	Eigen::ArrayX<T> row_sum_abs(const Eigen::SparseMatrix<T>& mat) {

		Eigen::ArrayX<T> res = Eigen::ArrayX<T>::Zero(mat.rows());

		for (int k = 0; k < mat.outerSize(); ++k) {

			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {

				res[it.row()] += std::abs(it.value());
			}
		}

		return res;
	}

	template <typename T>
	Eigen::ArrayX<T> col_sliced_row_sum(const Eigen::SparseMatrix<T>& mat, const Eigen::ArrayX<bool>& slice) {
		Eigen::ArrayX<T> res = Eigen::ArrayX<T>::Zero(mat.rows());
		const int ncol = mat.cols();
		for (int k = 0; k < ncol; ++k) {
			if (!slice[k]) {
				continue;
			}
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
				res[it.row()] += it.value();
			}
		}
		return res;
	}

	template <typename T>
	Eigen::ArrayXd col_sliced_row_mean(const Eigen::SparseMatrix<T>& mat, const Eigen::ArrayX<bool>& slice) {
		Eigen::ArrayXd res = Eigen::ArrayXd::Zero(mat.rows());
		const int ncol = mat.cols();
		for (int k = 0; k < ncol; ++k) {
			if (!slice[k]) {
				continue;
			}
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
				res[it.row()] += it.value();
			}
		}
		return res / slice.count();
	}

	template <typename T>
	Eigen::ArrayXd row_mean(const Eigen::SparseMatrix<T>& mat) {

		return row_sum(mat).template cast<double>() / mat.cols();
	}

	template <typename T>
	Eigen::ArrayXd col_mean(const Eigen::SparseMatrix<T>& mat) {

		return col_sum(mat).template cast<double>() / mat.rows();
	}

	template <typename T>
	Eigen::ArrayXd col_mean_mt(const Eigen::SparseMatrix<T>& mat) {

		return col_sum_mt(mat).template cast<double>() / mat.rows();
	}

	template <typename T>
	Eigen::ArrayXd row_var(const Eigen::SparseMatrix<T>& mat) {
		const int nrow = mat.rows();
		const int ncol = mat.cols();

		Eigen::ArrayXd vars = Eigen::ArrayXd::Zero(nrow), means = row_mean(mat);
		Eigen::ArrayXi n_zero = Eigen::ArrayXi::Constant(nrow, ncol);

		for (int i = 0; i < ncol; ++i) {
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it)
			{
				int row = it.row();
				double val = it.value(), mean = means[row];

				--n_zero[row];
				vars[row] += (val - mean) * (val - mean);
			}
		}

		vars += means * means * n_zero.cast<double>();
		vars /= (ncol - 1);

		return vars;
	}

	template <typename T>
	Eigen::ArrayXd col_var_mt(const Eigen::SparseMatrix<T>& mat) {
		const int nrow = mat.rows();
		const int ncol = mat.cols();
		Eigen::ArrayXd vars(ncol);
	#pragma omp parallel for
		for (int i = 0; i < ncol; ++i) {
			double col_mean = 0;
			double col_var = 0;
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it)
			{
				col_mean += it.value();
			}
			col_mean /= ncol;
			int n_not_zero = 0;

			for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it)
			{
				++n_not_zero;
				col_var += std::pow((it.value() - col_mean), 2);
			}
			col_var += std::pow(col_mean, 2) * (nrow - n_not_zero);
			col_var /= (nrow - 1);
			vars[i] = col_var;
		}

		return vars;
	}

	// default : count n_elements >= threshold per row
	template <typename T, bool MoreThan = true, bool Equal = true>
	Eigen::ArrayXi row_count(const Eigen::SparseMatrix<T>& mat, double threshold) {
		Eigen::ArrayXi res = Eigen::ArrayXi::Zero(mat.rows());
		if constexpr (Equal) {
			for (int k = 0; k < mat.outerSize(); ++k) {
				for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
					if ((it.value() >= threshold && MoreThan) || (it.value() <= threshold && !MoreThan)) {
						++res[it.row()];
					}
				}
			}
		}
		else {
			for (int k = 0; k < mat.outerSize(); ++k) {
				for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
					if ((it.value() > threshold && MoreThan) || (it.value() < threshold && !MoreThan)) {
						++res[it.row()];
					}
				}
			}
		}
		return res;
	}

	template<typename Con>
		requires is_container<Con>
	bool is_unique(const Con& ct) {

		return _Cs unique_element_number(ct) == std::ranges::size(ct);
	}

	template <typename Ele, typename S>
		requires is_slice_container<S>
	[[nodiscard]] Eigen::MatrixX<Ele>
		row_sliced(const Eigen::MatrixX<Ele>& mat, const S& slice) {

		const int res_row = std::ranges::count(slice, true);
		const int ncol = mat.cols();

		Eigen::MatrixX<Ele> res(res_row, ncol);

		int i{ 0 }, index{ 0 };

		for (auto&& ele : slice) {

			if (ele) {
				res.row(index++) = mat.row(i);
			}

			++i;
		}

		return res;
	}

	template <typename Ele, typename S>
		requires  is_slice_container<S>
	[[nodiscard]] Eigen::MatrixX<Ele>
		col_sliced(const Eigen::MatrixX<Ele>& mat, const S& slice) {

		const int res_col = std::ranges::count(slice, true);
		const int nrow = mat.rows();

		Eigen::MatrixX<Ele> res(nrow, res_col);

		int i{ 0 }, index{ 0 };

		for (auto&& ele : slice) {

			if (ele) {
				res.col(index++) = mat.col(i);
			}

			++i;
		}

		return res;
	}

	template <typename Ele, typename S>
		requires  is_slice_container<S>
	[[nodiscard]] Eigen::ArrayXX<Ele>
		row_sliced(const Eigen::ArrayXX<Ele>& mat, const S& slice) {

		const int res_row = std::ranges::count(slice, true);
		const int ncol = mat.cols();

		Eigen::ArrayXX<Ele> res(res_row, ncol);

		int i{ 0 }, index{ 0 };

		for (auto&& ele : slice) {

			if (ele) {
				res.row(index++) = mat.row(i);
			}

			++i;
		}

		return res;
	}

	template <typename Ele, typename S>
		requires  is_slice_container<S>
	[[nodiscard]] Eigen::ArrayXX<Ele>
		col_sliced(const Eigen::ArrayXX<Ele>& mat, const S& slice) {

		const int res_col = std::ranges::count(slice, true);
		const int nrow = mat.rows();

		Eigen::ArrayXX<Ele> res(nrow, res_col);

		int i{ 0 }, index{ 0 };

		for (auto&& ele : slice) {

			if (ele) {
				res.col(index++) = mat.col(i);
			}

			++i;
		}

		return res;
	}

	template <typename Ele, typename S1, typename S2>
		requires  is_slice_container<S1>&& is_slice_container<S2>
	[[nodiscard]] Eigen::ArrayXX<Ele>
		sliced(const Eigen::ArrayXX<Ele>& mat, const S1& row_slice, const S2& col_slice) {

		return _Cs row_sliced(_Cs col_sliced(mat, col_slice), row_slice);
	}

	template <typename Ele, typename S1, typename S2>
		requires  is_slice_container<S1>&& is_slice_container<S2>
	[[nodiscard]] Eigen::MatrixX<Ele>
		sliced(const Eigen::MatrixX<Ele>& mat, const S1& row_slice, const S2& col_slice) {

		return _Cs row_sliced(_Cs col_sliced(mat, col_slice), row_slice);
	}

	template <typename Ele, typename S>
		requires  is_slice_container<S>
	[[nodiscard]]
	Eigen::SparseMatrix<Ele>
		col_sliced(const Eigen::SparseMatrix<Ele>& mat, const S& slice)
	{

		const std::size_t res_col = std::ranges::count(slice, true);
		const std::size_t nrow = mat.rows();

		Eigen::SparseMatrix<Ele> res(nrow, res_col);
		std::vector<Eigen::Triplet<Ele> > triplets;

		auto iter = std::cbegin(slice);
		auto end = std::cend(slice);

		std::size_t target_col = 0, index_col = 0;

		for (; iter != end; ++iter, ++index_col) {
			if (*iter) {
				for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, index_col); it; ++it) {
					triplets.emplace_back(it.row(), target_col, it.value());
				}
				++target_col;
			}
		}

		res.setFromTriplets(triplets.cbegin(), triplets.cend());

		return res;
	};

	template <typename Ele, typename Order>
		requires  is_order_container<Order>
	[[nodiscard]] Eigen::SparseMatrix<Ele>
		col_reordered(const Eigen::SparseMatrix<Ele>& mat, const Order& order) {

		const std::size_t ncol = std::ranges::size(order);
		const std::size_t nrow = mat.rows();

		Eigen::SparseMatrix<Ele> res(nrow, ncol);
		std::vector<Eigen::Triplet<Ele> > triplets;

		auto iter = std::cbegin(order);
		auto end = std::cend(order);

		for (std::size_t i = 0; iter != end; ++iter, ++i) {

			std::size_t col = *iter;

			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, col); it; ++it) {
				triplets.emplace_back(it.row(), i, it.value());
			}
		}
		res.setFromTriplets(triplets.cbegin(), triplets.cend());

		return res;
	};

	template <typename Ele, typename S>
		requires  is_slice_container<S>
	[[nodiscard]] Eigen::SparseMatrix<Ele>
		row_sliced(const Eigen::SparseMatrix<Ele>& mat, const S& slice) {

		const int res_row = std::ranges::count(slice, true);
		const int nrow = mat.rows();
		const int ncol = mat.cols();

		std::vector<signed long long> row_map(nrow, -1);

		int index = 0;

		auto iter = std::cbegin(slice);
		auto end = std::cend(slice);

		auto index_iter = row_map.begin();

		for (; iter != end; ++iter, ++index_iter) {
			if (*iter) {
				*index_iter = index;
				++index;
			}
		}

		Eigen::SparseMatrix<Ele> res(res_row, ncol);
		std::vector<Eigen::Triplet<Ele> > triplets;

		for (int i = 0; i < ncol; ++i) {
			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, i); it; ++it) {
				signed long long target_row = row_map[it.row()];
				if (target_row >= 0) {
					triplets.emplace_back(target_row, i, it.value());
				}
			}
		}

		res.setFromTriplets(triplets.cbegin(), triplets.cend());

		return res;
	};

	template <typename Ele, typename Order>
		requires  is_order_container<Order>
	[[nodiscard]] Eigen::SparseMatrix<Ele>
		row_reordered(const Eigen::SparseMatrix<Ele>& mat, const Order& order) {

		const std::size_t res_row = std::ranges::size(order);
		const std::size_t ncol = mat.cols();
		const std::size_t nrow = mat.rows();

		Eigen::SparseMatrix<Ele> res(res_row, ncol);
		std::vector<Eigen::Triplet<Ele> > triplets;

		std::vector<signed long long> row_map(nrow, -1);

		auto iter = std::cbegin(order);
		auto end = std::cend(order);

		for (std::size_t i = 0; iter != end; ++iter, ++i) {
			row_map[*iter] = i;
		}

		for (std::size_t i = 0; i < ncol; ++i) {

			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, i); it; ++it) {
				signed long long target_row = row_map[it.row()];
				if (target_row >= 0) {
					triplets.emplace_back(target_row, i, it.value());
				}
			}
		}

		res.setFromTriplets(triplets.cbegin(), triplets.cend());

		return res;
	};

	template <typename Ele, typename S1, typename S2>
		requires  is_slice_container<S1>&& is_slice_container<S2>
	[[nodiscard]] Eigen::SparseMatrix<Ele> sliced(
		const Eigen::SparseMatrix<Ele>& mat,
		const S1& row_slice,
		const S2& col_slice)
	{
		const std::size_t res_row = std::ranges::count(row_slice, true);
		const std::size_t res_col = std::ranges::count(col_slice, true);

		const std::size_t nrow = mat.rows();
		const std::size_t ncol = mat.cols();

		Eigen::SparseMatrix<Ele> res(res_row, res_col);

		std::vector<signed long long> row_map(nrow, -1);

		std::size_t row_index = 0;

		auto row_iter = std::cbegin(row_slice);
		auto row_end = std::cend(row_slice);

		auto row_index_iter = row_map.begin();

		for (; row_iter != row_end; ++row_iter, ++row_index_iter) {
			if (*row_iter) {
				*row_index_iter = row_index;
				++row_index;
			}
		}

		std::vector<Eigen::Triplet<Ele> > triplets;

		auto col_iter = std::cbegin(col_slice);
		auto col_end = std::cend(col_slice);

		std::size_t target_col = 0;

		for (std::size_t i = 0; col_iter != col_end; ++col_iter, ++i) {

			if (*col_iter) {
				for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, i); it; ++it) {
					signed long long target_row = row_map[it.row()];
					if (target_row >= 0) {
						triplets.emplace_back(target_row, target_col, it.value());
					}
				}
				++target_col;
			}
		}

		res.setFromTriplets(triplets.cbegin(), triplets.cend());

		return res;
	};

	template <typename Ele, typename Order1, typename Order2>
		requires is_order_container<Order1>&& is_order_container<Order2>
	[[nodiscard]] Eigen::SparseMatrix<Ele> reordered(
		const Eigen::SparseMatrix<Ele>& mat,
		const Order1& row_order,
		const Order2& col_order)
	{
		const std::size_t res_row = std::ranges::size(row_order);
		const std::size_t res_col = std::ranges::size(col_order);

		const std::size_t nrow = mat.rows();
		const std::size_t ncol = mat.cols();

		Eigen::SparseMatrix<Ele> res(res_row, res_col);

		std::vector<signed long long> row_map(nrow, -1);

		std::size_t row_index = 0;

		auto row_iter = std::cbegin(row_order);
		auto row_end = std::cend(row_order);

		auto row_index_iter = row_map.begin();

		for (std::size_t i = 0; row_iter != row_end; ++row_iter, ++i) {
			row_map[*row_iter] = i;
		}

		std::vector<Eigen::Triplet<Ele> > triplets;

		auto col_iter = std::cbegin(col_order);
		auto col_end = std::cend(col_order);

		std::size_t target_col = 0;

		for (std::size_t i = 0; col_iter != col_end; ++col_iter, ++i) {

			std::size_t from_col = *col_iter;

			for (typename Eigen::SparseMatrix<Ele>::InnerIterator it(mat, from_col); it; ++it) {
				signed long long target_row = row_map[it.row()];
				if (target_row >= 0) {
					triplets.emplace_back(target_row, i, it.value());
				}
			}
		}

		res.setFromTriplets(triplets.cbegin(), triplets.cend());

		return res;
	}

	// x %in% y
	template<typename Con1, typename Con2>
		requires same_element<Con1, Con2>
	[[nodiscard]] Eigen::ArrayX<bool> in(Con1&& x, Con2&& y) {

		using ReCon1 = std::decay_t<Con1>;
		using ReCon2 = std::decay_t<Con2>;

		ReCon1 _x(std::forward<Con1>(x));
		ReCon2 _y(std::forward<Con2>(y));

		auto order = _Cs order(_x);
		_x = _Cs reordered(_x, order);
		_Cs sort(_y);

		auto iter_x = std::cbegin(_x);
		auto end_x = std::cend(_x);

		auto iter_y = std::cbegin(_y);
		auto end_y = std::cend(_y);

		const std::size_t size_x = std::ranges::size(_x);
		Eigen::ArrayX<bool> ret = Eigen::ArrayX<bool>::Constant(size_x, false);

		std::size_t index_x = 0;
		while (iter_x != end_x && iter_y != end_y) {
			if (*iter_x < *iter_y) {
				++iter_x;
				++index_x;
			}
			else if (*iter_x == *iter_y) {
				ret[order[index_x]] = true;
				++iter_x;
				++index_x;
			}
			else {
				++iter_y;
			}
		}

		return ret;
	}

	// only used when elements are unique
	template <typename Con>
		requires is_random_access_container<Con>
	[[nodiscard]] Con intersect(const Con& ct1, const Con& ct2) {

		auto order1 = _Cs order(ct1);
		auto order2 = _Cs order(ct2);

		const std::size_t size1 = std::ranges::size(ct1);
		const std::size_t size2 = std::ranges::size(ct2);

		Con ret;
		auto bi = std::back_inserter(ret);

		std::size_t iter1 = 0, iter2 = 0;
		while (iter1 < size1 && iter2 < size2) {
			if (ct1[order1[iter1]] < ct2[order2[iter2]]) {
				++iter1;
			}
			else if (ct1[order1[iter1]] > ct2[order2[iter2]]) {
				++iter2;
			}
			else
			{
				bi = ct1[order1[iter1]];
				++iter1;
				++iter2;
			}
		}

		return ret;
	}

	template <typename Con, typename Ele = get_element_raw_type<Con> >
		requires  is_container<Con>
	auto table(const Con& ct) {

		QMap<Ele, int> ret;

		for (auto&& ele : ct) {
			++ret[ele];
		}

		return ret;
	}

	template <typename Con, typename Ele = get_element_raw_type<Con> >
		requires  is_container<Con>&& less_comparable<Ele>
	Con cummin(const Con& ct) {

		const std::size_t size = std::ranges::size(ct);

		if (size < 2) {
			return ct;
		}

		Con ret(size);

		auto iter = std::cbegin(ct);

		auto lag = std::begin(ret);
		*lag = *iter;
		++iter;

		auto lead = std::next(lag);
		auto end = std::end(ret);


		for (; lead != end; ++lead, ++lag, ++iter) {
			if (*iter < *lag) {
				*lead = *iter;
			}
			else {
				*lead = *lag;
			}
		}

		return ret;
	}

	template <typename Con1, typename Con2>
		requires is_specific_container<Con1, double>&& is_specific_container<Con2, double>
	double correlation_pearson(const Con1& ct1, const Con2& ct2) {

		auto iter1 = std::cbegin(ct1);
		auto end1 = std::cend(ct1);

		auto iter2 = std::cbegin(ct2);
		auto end2 = std::cend(ct2);

		double sum_x = 0, sum_y = 0, sum_x_y = 0, sum_square_x = 0, sum_square_y = 0;

		for (; iter1 != end1; ++iter1, ++iter2) {
			double x = *iter1, y = *iter2;
			sum_x += x;
			sum_y += y;
			sum_x_y += x * y;
			sum_square_x += x * x;
			sum_square_y += y * y;
		}

		const std::size_t size = std::ranges::size(ct1);
		if (size * sum_square_x == sum_x * sum_x || size * sum_square_y == sum_y * sum_y)
			return std::nan(0);

		return (size * sum_x_y - sum_x * sum_y) / std::sqrt((size * sum_square_x - sum_x * sum_x) * (size * sum_square_y - sum_y * sum_y));
	}

	template <typename Con1, typename Con2>
		requires is_specific_container<Con1, double>&& is_specific_container<Con2, double>
	double correlation_spearman(const Con1& ct1, const Con2& ct2) {

		auto [rank_x, tie1, n_tie1] = _Cs rank_average(ct1);
		auto [rank_y, tie2, n_tie2] = _Cs rank_average(ct2);

		auto iter1 = rank_x.cbegin();
		auto end1 = rank_x.cend();

		auto iter2 = rank_y.cbegin();

		if (tie1 || tie2) {

			const double size = rank_x.size();

			double rank_x_average = _Cs sum(rank_x) / size;
			double rank_y_average = _Cs sum(rank_y) / size;

			double numerator{ 0.0 }, denominator1{ 0.0 }, denominator2{ 0.0 };

			for (; iter1 != end1; ++iter1, ++iter2) {
				numerator += (*iter1 - rank_x_average) * (*iter2 - rank_y_average);
				denominator1 += (*iter1 - rank_x_average) * (*iter1 - rank_x_average);
				denominator2 += (*iter2 - rank_y_average) * (*iter2 - rank_y_average);
			}
			if (denominator1 == 0 || denominator2 == 0) {
				return std::nan(0);
			}

			return numerator / std::sqrt(denominator1 * denominator2);
		}

		double di2{ 0.0 };

		for (; iter1 != end1; ++iter1, ++iter2) {
			double dis = (*iter1 - *iter2);
			di2 += dis * dis;
		}

		const double size = std::ranges::size(ct1);

		return 1 - (6 * di2) / (size * size * size - size);
	}
};

