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
SOAP_INLINE void sread(std::ifstream& in, T& val, const QString& edition = {});
template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val, const QString& edition = {});

template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val, const QString& edition = {});
template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val, const QString& edition = {});

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QSet<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QSet<T>& val, const QString& edition);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::SparseMatrix<T>& val, const QString& edition = {});
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::SparseMatrix<T>& val, const QString& edition = {});

template <typename... Types>
void sread(std::ifstream& in, std::tuple<Types...>& val, const QString& edition);
template <typename... Types>
void swrite(std::ofstream& out, const std::tuple<Types...>& val, const QString& edition);

template <typename T1, typename T2>
void sread(std::ifstream& in, std::pair<T1, T2>& val, const QString& edition);
template <typename T1, typename T2>
void swrite(std::ofstream& out, const std::pair<T1, T2>& val, const QString& edition);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, std::vector<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const std::vector<T>& val, const QString& edition);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QVector<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QVector<T>& val, const QString& edition);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, RunLengthEncoding<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const RunLengthEncoding<T>& val, const QString& edition);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, SparseVector<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const SparseVector<T>& val, const QString& edition);

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val, const QString& edition);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val, const QString& edition = {});
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val, const QString& edition = {});

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val, const QString& edition);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val, const QString& edition = {});
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val, const QString& edition = {});

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val, const QString& edition);
template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val, const QString& edition);

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val, const QString& edition = {});
template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val, const QString& edition = {});

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, std::map<T1, T2>& val, const QString& edition);
template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const std::map<T1, T2>& val, const QString& edition);

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, QMap<T1, T2>& val, const QString& edition);
template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const QMap<T1, T2>& val, const QString& edition);

void sread(std::ifstream& in, igraph_vector_int_t& val, const QString& edition = {});
void swrite(std::ofstream& out, const igraph_vector_int_t& val, const QString& edition = {});

void sread(std::ifstream& in, igraph_vector_t& val, const QString& edition = {});
void swrite(std::ofstream& out, const igraph_vector_t& val, const QString& edition = {});

void sread(std::ifstream& in, QColor& val, const QString& edition = {});
void swrite(std::ofstream& out, const QColor& val, const QString& edition = {});

void sread(std::ifstream& in, QString& val, const QString& edition = {});
void swrite(std::ofstream& out, const QString& val, const QString& edition = {});

void sread(std::ifstream& in, QFont& val, const QString& edition);
void swrite(std::ofstream& out, const QFont& val, const QString& edition);

void sread(std::ifstream& in, std::string& val, const QString& edition = {});
void swrite(std::ofstream& out, const std::string& val, const QString& edition = {});

void sread(std::ifstream& in, igraph_t& val, const QString& edition = {});
void swrite(std::ofstream& out, const igraph_t& val, const QString& edition = {});

void sread(std::ifstream& in, CustomMatrix& val, const QString& edition);
void swrite(std::ofstream& out, const CustomMatrix& val, const QString& edition);

void sread(std::ifstream& in, PatternWeightMatrix& val, const QString& edition);
void swrite(std::ofstream& out, const PatternWeightMatrix& val, const QString& edition);

void sread(std::ifstream& in, NumericMatrix& val, const QString& edition);
void swrite(std::ofstream& out, const NumericMatrix& val, const QString& edition);

void sread(std::ifstream& in, DataFrame& val, const QString& edition);
void swrite(std::ofstream& out, const DataFrame& val, const QString& edition);

void sread(std::ifstream& in, Enrichment& val, const QString& edition);
void swrite(std::ofstream& out, const Enrichment& val, const QString& edition);

void sread(std::ifstream& in, Metadata& val, const QString& edition);
void swrite(std::ofstream& out, const Metadata& val, const QString& edition);

void sread(std::ifstream& in, SparseInt& val, const QString& edition);
void swrite(std::ofstream& out, const SparseInt& val, const QString& edition);

void sread(std::ifstream& in, SparseDouble& val, const QString& edition);
void swrite(std::ofstream& out, const SparseDouble& val, const QString& edition);

void sread(std::ifstream& in, DenseInt& val, const QString& edition);
void swrite(std::ofstream& out, const DenseInt& val, const QString& edition);

void sread(std::ifstream& in, DenseDouble& val, const QString& edition);
void swrite(std::ofstream& out, const DenseDouble& val, const QString& edition);

void sread(std::ifstream& in, GSEA& val, const QString& edition);
void swrite(std::ofstream& out, const GSEA& val, const QString& edition);

void sread(std::ifstream& in, Embedding& val, const QString& edition);
void swrite(std::ofstream& out, const Embedding& val, const QString& edition);

void sread(std::ifstream& in, DifferentialAnalysis& val, const QString& edition);
void swrite(std::ofstream& out, const DifferentialAnalysis& val, const QString& edition);

void sread(std::ifstream& in, GeneName& val, const QString& edition);
void swrite(std::ofstream& out, const GeneName& val, const QString& edition);

void sread(std::ifstream& in, Footprint& val, const QString& edition);
void swrite(std::ofstream& out, const Footprint& val, const QString& edition);

void sread(std::ifstream& in, ChromVAR& val, const QString& edition);
void swrite(std::ofstream& out, const ChromVAR& val, const QString& edition);

void sread(std::ifstream& in, MotifPosition::match& val, const QString& edition);
void swrite(std::ofstream& out, const MotifPosition::match& val, const QString& edition);

void sread(std::ifstream& in, MotifPosition& val, const QString& edition);
void swrite(std::ofstream& out, const MotifPosition& val, const QString& edition);

void sread(std::ifstream& in, StringVector& val, const QString& edition);
void swrite(std::ofstream& out, const StringVector& val, const QString& edition);

void sread(std::ifstream& in, CNV& val, const QString& edition);
void swrite(std::ofstream& out, const CNV& val, const QString& edition);

void sread(std::ifstream& in, CellChat& val, const QString& edition);
void swrite(std::ofstream& out, const CellChat& val, const QString& edition);

void sread(std::ifstream& in, IRange& val, const QString& edition);
void swrite(std::ofstream& out, const IRange& val, const QString& edition);

void sread(std::ifstream& in, SeqInfo& val, const QString& edition);
void swrite(std::ofstream& out, const SeqInfo& val, const QString& edition);

void sread(std::ifstream& in, Fragments& val, const QString& edition);
void swrite(std::ofstream& out, const Fragments& val, const QString& edition);

void sread(std::ifstream& in, GenomicRange& val, const QString& edition);
void swrite(std::ofstream& out, const GenomicRange& val, const QString& edition);

void sread(std::ifstream& in, CoverageTrack& val, const QString& edition);
void swrite(std::ofstream& out, const CoverageTrack& val, const QString& edition);

void sread(std::ifstream& in, VelocityEstimate& val, const QString& edition);
void swrite(std::ofstream& out, const VelocityEstimate& val, const QString& edition);

void sread(std::ifstream& in, ScveloEstimate& val, const QString& edition);
void swrite(std::ofstream& out, const ScveloEstimate& val, const QString& edition);

void sread(std::ifstream& in, VelocytoBase& val, const QString& edition);
void swrite(std::ofstream& out, const VelocytoBase& val, const QString& edition);

void sread(std::ifstream& in, Pando& val, const QString& edition);
void swrite(std::ofstream& out, const Pando& val, const QString& edition);

void sread(std::ifstream& in, Monocle3& val, const QString& edition);
void swrite(std::ofstream& out, const Monocle3& val, const QString& edition);

void sread(std::ifstream& in, Cicero& val, const QString& edition);
void swrite(std::ofstream& out, const Cicero& val, const QString& edition);

void sread(std::ifstream& in, GraphSettings& val, const QString& edition);
void swrite(std::ofstream& out, const GraphSettings& val, const QString& edition);

void sread(std::ifstream& in, DataField& val, const QString& edition);
void swrite(std::ofstream& out, const DataField& val, const QString& edition);

void sread(std::ifstream& in, BulkRna& val, const QString& edition);
void swrite(std::ofstream& out, const BulkRna& val, const QString& edition);

void sread(std::ifstream& in, SingleCellRna& val, const QString& edition);
void swrite(std::ofstream& out, const SingleCellRna& val, const QString& edition);

void sread(std::ifstream& in, SingleCellAtac& val, const QString& edition);
void swrite(std::ofstream& out, const SingleCellAtac& val, const QString& edition);

void sread(std::ifstream& in, SingleCellMultiome& val, const QString& edition);
void swrite(std::ofstream& out, const SingleCellMultiome& val, const QString& edition);

template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val, const QString& edition) {

	in.read(reinterpret_cast<char*>(&val), sizeof(T));
}

template <typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val, const QString& edition) {

	out.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void sread(std::ifstream& in, T& val, const QString& edition) {
	int enum_val{ 0 };
	sread(in, enum_val, edition);
	val = static_cast<T>(enum_val);
}

template <typename T>
	requires std::is_enum_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const T& val, const QString& edition) {
	int enum_val = static_cast<int>(val);
	swrite(out, enum_val, edition);
}

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QSet<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	T ele;
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, ele);
		val.insert(ele);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QSet<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	for (auto&& ele : val) {
		swrite(out, ele);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::SparseMatrix<T>& val, const QString& edition) {
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
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::SparseMatrix<T>& val, const QString& edition) {
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
void sread(std::ifstream& in, std::tuple<Types...>& val, const QString& edition)
{
	std::apply([&](auto &...element) { (sread(in, element, edition), ...); }, val);
}

template <typename... Types>
void swrite(std::ofstream& out, const std::tuple<Types...>& val, const QString& edition)
{
	std::apply([&](const auto &...element) { (swrite(out, element, edition), ...); }, val);
}

template <typename T1, typename T2>
void sread(std::ifstream& in, std::pair<T1, T2>& val, const QString& edition)
{
	sread(in, val.first, edition);
	sread(in, val.second, edition);
}

template <typename T1, typename T2>
void swrite(std::ofstream& out, const std::pair<T1, T2>& val, const QString& edition)
{
	swrite(out, val.first, edition);
	swrite(out, val.second, edition);
}

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, std::vector<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, val[i], edition);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const std::vector<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, val[i], edition);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, std::vector<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const std::vector<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, QVector<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, val[i], edition);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const QVector<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, val[i], edition);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, QVector<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const QVector<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, RunLengthEncoding<T>& val, const QString& edition) {
	sread(in, val.data_, edition);
	sread(in, val.index_, edition);
	sread(in, val.size_, edition);
	sread(in, val.encoding_size_, edition);
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const RunLengthEncoding<T>& val, const QString& edition) {
	swrite(out, val.data_, edition);
	swrite(out, val.index_, edition);
	swrite(out, val.size_, edition);
	swrite(out, val.encoding_size_, edition);
};

template<typename T>
SOAP_INLINE void sread(std::ifstream& in, SparseVector<T>& val, const QString& edition) {
	sread(in, val.data_, edition);
	sread(in, val.index_, edition);
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const SparseVector<T>& val, const QString& edition) {
	swrite(out, val.data_, edition);
	swrite(out, val.index_, edition);
};


template<typename T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val, const QString& edition) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);
	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			sread(in, val(i, j), edition);
		}
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val, const QString& edition) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	nrow = val.rows();
	ncol = val.cols();
	size = nrow * ncol;
	swrite(out, size);
	swrite(out, nrow);
	swrite(out, ncol);

	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			swrite(out, val(i, j), edition);
		}
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::MatrixX<T>& val, const QString& edition) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);

	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::MatrixX<T>& val, const QString& edition) {
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
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val, const QString& edition) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);
	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			sread(in, val(i, j), edition);
		}
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val, const QString& edition) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	nrow = val.rows();
	ncol = val.cols();
	size = nrow * ncol;
	swrite(out, size);
	swrite(out, nrow);
	swrite(out, ncol);

	for (std::size_t j = 0; j < ncol; ++j) {
		for (std::size_t i = 0; i < nrow; ++i) {
			swrite(out, val(i, j), edition);
		}
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayXX<T>& val, const QString& edition) {
	std::size_t nrow{ 0 }, ncol{ 0 }, size{ 0 };
	sread(in, size);
	sread(in, nrow);
	sread(in, ncol);
	val.resize(nrow, ncol);

	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayXX<T>& val, const QString& edition) {
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
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, val[i], edition);
	}
};

template<typename T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, val[i], edition);
	}
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void sread(std::ifstream& in, Eigen::ArrayX<T>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(T));
};

template<typename T>
	requires std::is_arithmetic_v<T>
SOAP_INLINE void swrite(std::ofstream& out, const Eigen::ArrayX<T>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);
	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(T));
};

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, std::map<T1, T2>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);

	T1 key;
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, key, edition);
		T2& value = val[key];
		sread(in, value, edition);
	}
};

template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const std::map<T1, T2>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);

	for (auto&& [key, value] : val) {
		swrite(out, static_cast<const T1&>(key), edition);
		swrite(out, static_cast<const T2&>(value), edition);
	}
};

template<typename T1, typename T2>
SOAP_INLINE void sread(std::ifstream& in, QMap<T1, T2>& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size);

	T1 key;
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, key, edition);
		sread(in, val[key], edition);
	}
};

template<typename T1, typename T2>
SOAP_INLINE void swrite(std::ofstream& out, const QMap<T1, T2>& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size);

	for (auto&& [key, value] : val.asKeyValueRange()) {
		swrite(out, key, edition);
		swrite(out, value, edition);
	}
};

