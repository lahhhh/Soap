#pragma once

#include "Identifier.h"

#include <QVector>
#include <QString>

class IRange
{
public:
	G_CLASS_FUNCTION_DEFAULT(IRange);

	IRange(const QVector<int>& start, const QVector<int>& width);

	QVector<int> start_;

	QVector<int> width_;

	void clear();

	void append(int start, int width);

	void append(const IRange& range);

	qsizetype size() const;

	G_SET_IDENTIFIER("IRange");

};

