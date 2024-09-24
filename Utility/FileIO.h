#pragma once

#include "Identifier.h"

#include <fstream>

#include <QFile>
#include <QTextStream>

#include "SingleCellMultiome.h"

CustomMatrix* read_table(const QString& file_name, bool fast = false);

CustomMatrix* read_sv(const QString& file_name, const QChar& delimiter = ',');

// Note : do not support unicode, char only
CustomMatrix* read_sv_fast(const QString& file_name, const char delimiter = ',');

void write_sv(
	const CustomMatrix& mat, 
	QTextStream& out, 
	const QChar delimiter,
	bool quote,
	bool keep_rownames,
	bool keep_colnames);

void write_sv(
	const DenseDouble& dd,
	QTextStream& out, 
	const QChar delimiter, 
	bool quote,
	bool keep_rownames,
	bool keep_colnames);

void write_sv(
	const MotifPosition& mp,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const GenomicRange& gr,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const Eigen::MatrixXd& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const Eigen::ArrayXXd& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const Eigen::SparseMatrix<double>& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const Eigen::MatrixXi& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const Eigen::ArrayXXi& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const Eigen::SparseMatrix<int>& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const QStringList& s,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames = {},
	const QStringList& colnames = {});

void write_sv(
	const QStringList& s,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	bool keep_rownames,
	bool keep_colnames);


std::map<QString, PatternWeightMatrix> read_motif_database(const QString& file_name);

template <typename MatrixType, typename TextStream>
void write_csv(const MatrixType& mat, TextStream& out, const QStringList& settings) {

	bool quote = !settings.contains(CSV_NO_QUOTE);
	bool keep_rownames = !settings.contains(CSV_NO_ROWNAMES);
	bool keep_colnames = !settings.contains(CSV_NO_COLNAMES);

	write_sv(mat, out, ',', quote, keep_rownames, keep_colnames);
}

template <typename MatrixType, typename TextStream>
void write_tsv(const MatrixType& mat, TextStream& out, const QStringList& settings) {

	bool quote = !settings.contains(TSV_NO_QUOTE);
	bool keep_rownames = !settings.contains(TSV_NO_ROWNAMES);
	bool keep_colnames = !settings.contains(TSV_NO_COLNAMES);

	write_sv(mat, out, '\t', quote, keep_rownames, keep_colnames);
}

//SparseInt read_sparse_int(const QString& fileName);
