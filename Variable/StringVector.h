#pragma once

#include "Identifier.h"

#include <QStringList>

#include "Custom.h"

class StringVector
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(StringVector);

	StringVector(const QStringList& str) : data_(str) {}

	DataType data_type_{ DataType::Plain };

	QStringList data_;

	G_SET_IDENTIFIER("StringVector");
};

