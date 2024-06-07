#pragma once

#include "Identifier.h"

#include <QStringList>

#include "Custom.h"
#include "Enrichment.h"

class GeneName
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(GeneName);

	explicit GeneName(const QStringList& gene_names) :
		data_(gene_names)
	{}

	DataType data_type_{ DataType::Plain };

	QStringList data_;

	SOAP_SUBMODULES(Enrichment);

	G_SET_IDENTIFIER("GeneName");
};

