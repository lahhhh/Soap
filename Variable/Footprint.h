#pragma once

#include "Identifier.h"

#include "DenseInt.h"
#include "PatternWeightMatrix.h"
#include "GraphSettings.h"
#include "GenomicRange.h"
#include "GeneName.h"

class Footprint
{
public:

	enum class DataType : int {	Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(Footprint);

	Footprint(
		const PatternWeightMatrix& motif,
		const Eigen::MatrixXi& insertion_matrix,
		const QStringList& cell_names,
		const QStringList& base_position_names,
		const QVector<double>& expected_insertions,
		const GenomicRange& motif_location
	);

	G_SET_IDENTIFIER("Footprint");

	DataType data_type_{ DataType::Plain };

	PatternWeightMatrix motif_;

	DenseInt insertion_matrix_;

	QVector<double> expected_insertions_;

	GenomicRange motif_location_;

	SOAP_SUBMODULES(GeneName);

	bool is_empty() const;

	bool draw(
		const QString& file_name,
		const QString& factor_name,
		const QStringList& levels,
		const QStringList& factors,
		const GraphSettings& gs,
		const int height,
		const int width
	) const;
};

