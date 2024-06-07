#pragma once

#include "Identifier.h"

#include <fstream>
#include <QFont>

#include "BulkRna.h"
#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"
#include "DataFrame.h"
#include "StringVector.h"
#include "NumericMatrix.h"
#include "GraphSettings.h"
#include "Custom.h"


template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val);
template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val);

template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val);
template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QSet<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QSet<T>& val);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::SparseMatrix<T>& val);
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::SparseMatrix<T>& val);

template <typename... Types>
void sread(std::ifstream& in, std::tuple<Types...>& val);
template <typename... Types>
void swrite(std::ofstream& out, const std::tuple<Types...>& val);

template <typename T1, typename T2>
void sread(std::ifstream& in, std::pair<T1, T2>& val);
template <typename T1, typename T2>
void swrite(std::ofstream& out, const std::pair<T1, T2>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, std::vector<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const std::vector<T>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QVector<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QVector<T>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, RunLengthEncoding<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const RunLengthEncoding<T>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, SparseVector<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const SparseVector<T>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val);
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val);
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val);
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val);

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, std::map<T1, T2>& val);
template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const std::map<T1, T2>& val);

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, QMap<T1, T2>& val);
template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const QMap<T1, T2>& val);

void sread(std::ifstream& in, igraph_vector_int_t& val);
void swrite(std::ofstream& out, const igraph_vector_int_t& val);

void sread(std::ifstream& in, igraph_vector_t& val);
void swrite(std::ofstream& out, const igraph_vector_t& val);

void sread(std::ifstream& in, QColor& val);
void swrite(std::ofstream& out, const QColor& val);

void sread(std::ifstream& in, QString& val);
void swrite(std::ofstream& out, const QString& val);

void sread(std::ifstream& in, QFont& val);
void swrite(std::ofstream& out, const QFont& val);

void sread(std::ifstream& in, std::string& val);
void swrite(std::ofstream& out, const std::string& val);

void sread(std::ifstream& in, igraph_t& val);
void swrite(std::ofstream& out, const igraph_t& val);

void sread(std::ifstream& in, CustomMatrix& val);
void swrite(std::ofstream& out, const CustomMatrix& val);

void sread(std::ifstream& in, PatternWeightMatrix& val);
void swrite(std::ofstream& out, const PatternWeightMatrix& val);

void sread(std::ifstream& in, NumericMatrix& val);
void swrite(std::ofstream& out, const NumericMatrix& val);

void sread(std::ifstream& in, DataFrame& val);
void swrite(std::ofstream& out, const DataFrame& val);

void sread(std::ifstream& in, Enrichment& val);
void swrite(std::ofstream& out, const Enrichment& val);

void sread(std::ifstream& in, Metadata& val);
void swrite(std::ofstream& out, const Metadata& val);

void sread(std::ifstream& in, SparseInt& val);
void swrite(std::ofstream& out, const SparseInt& val);

void sread(std::ifstream& in, SparseDouble& val);
void swrite(std::ofstream& out, const SparseDouble& val);

void sread(std::ifstream& in, DenseInt& val);
void swrite(std::ofstream& out, const DenseInt& val);

void sread(std::ifstream& in, DenseDouble& val);
void swrite(std::ofstream& out, const DenseDouble& val);

void sread(std::ifstream& in, GSEA& val);
void swrite(std::ofstream& out, const GSEA& val);

void sread(std::ifstream& in, Embedding& val);
void swrite(std::ofstream& out, const Embedding& val);

void sread(std::ifstream& in, DifferentialAnalysis& val);
void swrite(std::ofstream& out, const DifferentialAnalysis& val);

void sread(std::ifstream& in, GeneName& val);
void swrite(std::ofstream& out, const GeneName& val);

void sread(std::ifstream& in, Footprint& val);
void swrite(std::ofstream& out, const Footprint& val);

void sread(std::ifstream& in, ChromVAR& val);
void swrite(std::ofstream& out, const ChromVAR& val);

void sread(std::ifstream& in, MotifPosition::match& val);
void swrite(std::ofstream& out, const MotifPosition::match& val);

void sread(std::ifstream& in, MotifPosition& val);
void swrite(std::ofstream& out, const MotifPosition& val);

void sread(std::ifstream& in, StringVector& val);
void swrite(std::ofstream& out, const StringVector& val);

void sread(std::ifstream& in, CNV& val);
void swrite(std::ofstream& out, const CNV& val);

void sread(std::ifstream& in, CellChat& val);
void swrite(std::ofstream& out, const CellChat& val);

void sread(std::ifstream& in, IRange& val);
void swrite(std::ofstream& out, const IRange& val);

void sread(std::ifstream& in, SeqInfo& val);
void swrite(std::ofstream& out, const SeqInfo& val);

void sread(std::ifstream& in, Fragments& val);
void swrite(std::ofstream& out, const Fragments& val);

void sread(std::ifstream& in, GenomicRange& val);
void swrite(std::ofstream& out, const GenomicRange& val);

void sread(std::ifstream& in, CoverageTrack& val);
void swrite(std::ofstream& out, const CoverageTrack& val);

void sread(std::ifstream& in, VelocityEstimate& val);
void swrite(std::ofstream& out, const VelocityEstimate& val);

void sread(std::ifstream& in, ScveloEstimate& val);
void swrite(std::ofstream& out, const ScveloEstimate& val);

void sread(std::ifstream& in, VelocytoBase& val);
void swrite(std::ofstream& out, const VelocytoBase& val);

void sread(std::ifstream& in, Pando& val);
void swrite(std::ofstream& out, const Pando& val);

void sread(std::ifstream& in, Monocle3& val);
void swrite(std::ofstream& out, const Monocle3& val);

void sread(std::ifstream& in, Cicero& val);
void swrite(std::ofstream& out, const Cicero& val);

void sread(std::ifstream& in, GraphSettings& val);
void swrite(std::ofstream& out, const GraphSettings& val);

void sread(std::ifstream& in, DataField& val);
void swrite(std::ofstream& out, const DataField& val);

void sread(std::ifstream& in, BulkRna& val);
void swrite(std::ofstream& out, const BulkRna& val);

void sread(std::ifstream& in, SingleCellRna& val);
void swrite(std::ofstream& out, const SingleCellRna& val);

void sread(std::ifstream& in, SingleCellAtac& val);
void swrite(std::ofstream& out, const SingleCellAtac& val);

void sread(std::ifstream& in, SingleCellMultiome& val);
void swrite(std::ofstream& out, const SingleCellMultiome& val);

template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val) {
	in.read(reinterpret_cast<char*>(&val), sizeof(T));
}

template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val) {
	out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val) {
	int enum_val{ 0 };
	sread(in, enum_val);
	val = static_cast<T>(enum_val);
}

template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val) {
	int enum_val = static_cast<int>(val);
	swrite(out, enum_val);
}

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QSet<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	T ele;
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, ele);
		val.insert(ele);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QSet<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	for (auto&& ele : val) {
		swrite(out, ele);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::SparseMatrix<T>& val) {
	int rows{ 0 }, cols{ 0 };
	std::size_t nnz{ 0 };
	sread(in, rows);
	sread(in, cols);
	sread(in, nnz);

	val.resize(rows, cols);
	val.makeCompressed();
	val.resizeNonZeros(nnz);

	in.read(reinterpret_cast<char*>(val.valuePtr()), nnz * sizeof(T));
	in.read(reinterpret_cast<char*>(val.innerIndexPtr()), nnz * sizeof(int));
	in.read(reinterpret_cast<char*>(val.outerIndexPtr()), (cols + 1) * sizeof(int));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::SparseMatrix<T>& val) {
	int rows = val.rows(), cols = val.cols();
	std::size_t nnz = val.nonZeros();
	swrite(out, rows);
	swrite(out, cols);
	swrite(out, nnz);

	out.write(reinterpret_cast<const char*>(val.valuePtr()), nnz * sizeof(T));
	out.write(reinterpret_cast<const char*>(val.innerIndexPtr()), nnz * sizeof(int));
	out.write(reinterpret_cast<const char*>(val.outerIndexPtr()), (cols + 1) * sizeof(int));
};


template <typename... Types>
void sread(std::ifstream& in, std::tuple<Types...>& val)
{
	std::apply([&](auto &...element) { (sread(in, element), ...); }, val);
}

template <typename... Types>
void swrite(std::ofstream& out, const std::tuple<Types...>& val)
{
	std::apply([&](const auto &...element) { (swrite(out, element), ...); }, val);
}

template <typename T1, typename T2>
void sread(std::ifstream& in, std::pair<T1, T2>& val)
{
	sread(in, val.first);
	sread(in, val.second);
}

template <typename T1, typename T2>
void swrite(std::ofstream& out, const std::pair<T1, T2>& val)
{
	swrite(out, val.first);
	swrite(out, val.second);
}

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, std::vector<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, val[i]);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const std::vector<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, val[i]);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, std::vector<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const std::vector<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QVector<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, val[i]);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QVector<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, val[i]);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, QVector<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const QVector<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, RunLengthEncoding<T>& val) {
	sread(in, val.data_);
	sread(in, val.index_);
	sread(in, val.size_);
	sread(in, val.encoding_size_);
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const RunLengthEncoding<T>& val) {
	swrite(out, val.data_);
	swrite(out, val.index_);
	swrite(out, val.size_);
	swrite(out, val.encoding_size_);
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, SparseVector<T>& val) {
	sread(in, val.data_);
	sread(in, val.index_);
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const SparseVector<T>& val) {
	swrite(out, val.data_);
	swrite(out, val.index_);
};


template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);
	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			sread(in, val(i, j));
		}
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	nrow = val.rows();
	ncol = val.cols();
	size = nrow * ncol;
	swrite(out, size);
	swrite(out, nrow);
	swrite(out, ncol);

	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			swrite(out, val(i, j));
		}
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);

	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	nrow = val.rows();
	ncol = val.cols();
	size = nrow * ncol;
	swrite(out, size);
	swrite(out, nrow);
	swrite(out, ncol);

	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);
	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			sread(in, val(i, j));
		}
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	nrow = val.rows();
	ncol = val.cols();
	size = nrow * ncol;
	swrite(out, size);
	swrite(out, nrow);
	swrite(out, ncol);

	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			swrite(out, val(i, j));
		}
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);

	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	nrow = val.rows();
	ncol = val.cols();
	size = nrow * ncol;
	swrite(out, size);
	swrite(out, nrow);
	swrite(out, ncol);

	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, val[i]);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, val[i]);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val) {
	std::size_t size = val.size();
	swrite(out, size);
	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, std::map<T1, T2>& val) {
	std::size_t size{ 0 };
	sread(in, size);

	T1 key;
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, key);
		T2& value = val[key];
		sread(in, value);
	}
};

template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const std::map<T1, T2>& val) {
	std::size_t size = val.size();
	swrite(out, size);

	for (auto&& [key, value] : val) {
		swrite(out, static_cast<const T1&>(key));
		swrite(out, static_cast<const T2&>(value));
	}
};

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, QMap<T1, T2>& val) {
	std::size_t size{ 0 };
	sread(in, size);

	T1 key;
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, key);
		sread(in, val[key]);
	}
};

template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const QMap<T1, T2>& val) {
	std::size_t size = val.size();
	swrite(out, size);

	for (auto&& [key, value] : val.asKeyValueRange()) {
		swrite(out, key);
		swrite(out, value);
	}
};

