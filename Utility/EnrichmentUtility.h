#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"

QMap<QString, QMap<QString, QStringList>> get_pathway_information(soap::Species species);

CustomMatrix enrich(
    const QString& database_name,
    const QStringList& gene_names,
    const QString& ontology,
    soap::Species species,
    const QString& adjust_p_value_method,
    double p_threshold
);

