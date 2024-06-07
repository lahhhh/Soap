#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"

class EnrichmentDatabase
{

public:

    EnrichmentDatabase() = default;
    EnrichmentDatabase(const EnrichmentDatabase&) = delete;

    CustomMatrix enrich_go(
        const QStringList & gene_names, 
        const QString & ontology, 
        const soap::Species &species,
        const QString &adjust_p_value_method,
        double p_threshold
    );

    CustomMatrix enrich_kegg(
        const QStringList & gene_names, 
        const soap::Species &species, 
        const QString &adjust_p_value_method, 
        double p_threshold
    );

    QMap<QString, QMap<QString, QStringList>> get_pathway_information(soap::Species species);

private:

    QMap<QString, QMap<QString, QStringList>> symbol_to_pathway_;

    QMap<QString, QMap<QString, QStringList>> pathway_to_symbol_;

    QMap<QString, QMap<QString, QString>> pathway_to_ontology_;

    QMap<QString, QMap<QString, QString>> pathway_to_pathway_name_;

    QMap<QString, QMap<QString, int>> pathway_to_size_;

private:

    void load_go_human();

    void load_go_mouse();

    void load_go_pathway_name();

    void load_go_path_to_symbol_human();

    void load_go_symbol_to_path_human();

    void load_go_path_to_symbol_mouse();

    void load_go_symbol_to_path_mouse();

    void load_kegg_human();

    void load_kegg_mouse();

    void load_kegg_pathway_name_human();

    void load_kegg_pathway_name_mouse();

    void load_kegg_path_to_symbol_human();

    void load_kegg_symbol_to_path_human();

    void load_kegg_path_to_symbol_mouse();

    void load_kegg_symbol_to_path_mouse();

};