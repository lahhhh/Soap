#pragma once

#include "Identifier.h"

#include "DenseInt.h"
#include "DenseDouble.h"
#include "Embedding.h"
#include "Metadata.h"
#include "DifferentialAnalysis.h"

class BulkRna
{
public:
	enum class DataType : int { Plain = 0, Integrated = 1 };

	G_CLASS_FUNCTION_DEFAULT(BulkRna);

	G_SET_IDENTIFIER("BulkRna");

	G_QUICK_ACCESS(Metadata, metadata);
	G_QUICK_ACCESS_TYPE(DenseInt, counts, Counts);
	G_QUICK_ACCESS_TYPE(DenseDouble, normalized, Normalized);
	G_QUICK_ACCESS_TYPE(Embedding, pca, Pca);
	G_QUICK_ACCESS_TYPE(Embedding, umap, Umap);
	G_QUICK_ACCESS_TYPE(Embedding, tsne, Tsne);

	/******       data        *****/

	DataType data_type_{ DataType::Plain };

	soap::Species species_ = soap::Species::Undefined;
	int random_state_ = 1997;

	std::map<QString, QString> string_information_;
	std::map<QString, int> integer_information_;
	std::map<QString, double> double_information_;

	std::map<QString, QStringList> string_vectors_;
	std::map<QString, QVector<int> > integer_vectors_;
	std::map<QString, QVector<double> > double_vectors_;

	SOAP_SUBMODULES(Metadata);
	SOAP_SUBMODULES(DenseInt);
	SOAP_SUBMODULES(DenseDouble);
	SOAP_SUBMODULES(Embedding);
	SOAP_SUBMODULES(DifferentialAnalysis);
	SOAP_SUBMODULES(GSEA);
};

