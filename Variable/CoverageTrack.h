#pragma once

#include "Identifier.h"

#include <QString>
#include <QMap>
#include <QVector>
#include <vector>

#include "GenomicRange.h"
#include "SparseVector.h"

struct Location {
	QString sequence_name;
	int start = 0;
	int end = 0;
};

class CoverageTrack
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(CoverageTrack);

	DataType data_type_{ DataType::Plain };

	QMap<QString, std::vector< SparseVector<float> > > insertion_matrix_;

	QString level_name_;
	QStringList levels_;

	GenomicRange annotation_;

	G_SET_IDENTIFIER("Coverage Track");
};

