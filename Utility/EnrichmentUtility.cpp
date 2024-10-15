#include "EnrichmentUtility.h"

#include <QFile>
#include "Custom.h"
#include "Phyper.h"

std::tuple<QMap<QString, QString>, QMap<QString, int>, QMap<QString, QStringList>> 
load_path_to_symbol(const QString& database_name, soap::Species species) {

	QString filename;

	if (database_name == "GO") {
		if (species == soap::Species::Human) {
			filename = FILE_GO_PATH2SYMBOL_HUMAN;
		}
		else {
			filename = FILE_GO_PATH2SYMBOL_MOUSE;
		}
	}
	else {
		if (species == soap::Species::Human) {
			filename = FILE_KEGG_PATH2SYMBOL_HUMAN;
		}
		else {
			filename = FILE_KEGG_PATH2SYMBOL_MOUSE;
		}
	}

	QFile file(filename);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QString> pathway_to_ontology;
	QMap<QString, int> pathway_to_size;
	QMap<QString, QStringList> pathway_to_symbol;

	QString line = in.readLine();

	while (!line.isNull()) {

		QStringList tmp = line.split(' ');
		QString path = tmp[0];

		if (database_name == "GO") {

			pathway_to_ontology[path] = tmp[1];
			pathway_to_symbol[path] = tmp.sliced(2);
			pathway_to_size[path] = tmp.size() - 2;
		}
		else {

			pathway_to_symbol[path] = tmp.sliced(1);
			pathway_to_size[path] = tmp.size() - 1;
		}

		line = in.readLine();
	}

	return { pathway_to_ontology, pathway_to_size, pathway_to_symbol };
};

QMap<QString, QString> load_pathway_name(const QString& database_name, soap::Species species) {

	QString filename;

	if (database_name == "GO") {
			filename = FILE_GO_PATH2NAME;
	}
	else {
		if (species == soap::Species::Human) {
			filename = FILE_KEGG_PATH2NAME_HUMAN;
		}
		else {
			filename = FILE_KEGG_PATH2NAME_MOUSE;
		}
	}

	QFile file(filename);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QString> pathway_to_name;
	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index_of_space = line.indexOf(' ');
		pathway_to_name[line.sliced(0, index_of_space)] = line.sliced(index_of_space + 1);

		line = in.readLine();
	}

	return pathway_to_name;
};

QMap<QString, QStringList> load_symbol_to_path(const QString& database_name, soap::Species species) {

	QString filename;

	if (database_name == "GO") {
		if (species == soap::Species::Human) {
			filename = FILE_GO_SYMBOL2PATH_HUMAN;
		}
		else {
			filename = FILE_GO_SYMBOL2PATH_MOUSE;
		}
	}
	else {
		if (species == soap::Species::Human) {
			filename = FILE_KEGG_SYMBOL2PATH_HUMAN;
		}
		else {
			filename = FILE_KEGG_SYMBOL2PATH_MOUSE;
		}
	}

	QFile file(filename);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QStringList> symbol_to_pathway;

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');

		symbol_to_pathway[tmp[0]] = tmp.sliced(1);

		line = in.readLine();
	}

	return symbol_to_pathway;
};

CustomMatrix enrich(
	const QString& database_name,
    const QStringList& gene_names,
    const QString& ontology,
    soap::Species species,
    const QString& adjust_p_value_method,
    double p_threshold
) {
    
	int minimum_pathway_size{ 10 }, maximum_pathway_size{ 500 };

	QStringList ret_path, path_id, query_genes, pathway_genes;
	QVector<int> count;
	QMap<QString, QString> query_map;
	QStringList gene_ratio, background_ratio;
	QVector<double> p_value, p_adjusted;

	auto [pathway_to_ontology, pathway_to_size, pathway_to_symbol] = load_path_to_symbol(database_name, species);
	auto pathway_to_name = load_pathway_name(database_name, species);
	auto symbol_to_pathway = load_symbol_to_path(database_name, species);

	QStringList filtered_genes = custom::intersect(gene_names, symbol_to_pathway.keys());

	auto filtered_size = filtered_genes.size();

	// no enough gene found, return empty matrix
	if (filtered_size < 10) {
		return {};
	}

	QStringList path_names;
	for (const auto& gene_name : filtered_genes) {
		path_names << symbol_to_pathway[gene_name];
		auto& paths = symbol_to_pathway[gene_name];
		for (QString& path_name : paths) {
			query_map[path_name] += (gene_name + ";");
		}
	}


	QMap<QString, int> query;
	// calculate query times for pathways
	for (auto& path_name : path_names) {
		++query[path_name];
	}

	QStringList unique_paths = custom::unique(path_names);
	QStringList filtered_paths;

	if (ontology == "ALL") {
		filtered_paths = unique_paths;
	}
	else {
		for (auto& path_name : unique_paths) {
			if (pathway_to_ontology[path_name] == ontology) {
				filtered_paths << path_name;
			}
		}
	}


	double database_size = symbol_to_pathway.size();

	for (auto& path_name : filtered_paths) {

		// check pathway size
		int pathway_size = pathway_to_size[path_name];
		if (pathway_size < minimum_pathway_size || pathway_size > maximum_pathway_size) {
			continue;
		}

		double p = phyper(
			query[path_name] - 1.0,
			(double)pathway_to_size[path_name],
			database_size - pathway_to_size[path_name],
			(double)filtered_size
		);

		path_id << path_name;
		ret_path << pathway_to_name[path_name];
		p_value << p;
		query_genes << query_map[path_name];
		pathway_genes << pathway_to_symbol[path_name].join(",");
		count << query[path_name];
		gene_ratio << QString::number(query[path_name]) + "/" + QString::number(filtered_size);
		background_ratio << QString::number(pathway_to_size[path_name]) + "/" + QString::number(database_size);
	}

	if (p_value.isEmpty()) {
		return {};
	}

	p_adjusted = custom::adjust_p_value(p_value, adjust_p_value_method);
	auto index = custom::order(p_adjusted);
	p_adjusted = custom::reordered(p_adjusted, index);
	path_id = custom::reordered(path_id, index);
	ret_path = custom::reordered(ret_path, index);
	p_value = custom::reordered(p_value, index);
	query_genes = custom::reordered(query_genes, index);
	pathway_genes = custom::reordered(pathway_genes, index);
	gene_ratio = custom::reordered(gene_ratio, index);
	background_ratio = custom::reordered(background_ratio, index);
	count = custom::reordered(count, index);

	CustomMatrix ret(path_id);

	ret.update(METADATA_ENRICHMENT_PATHWAY_NAMES, ret_path);
	ret.update(METADATA_ENRICHMENT_P_VALUE, p_value);
	ret.update(METADATA_ENRICHMENT_ADJUSTED_P_VALUE, p_adjusted);
	ret.update(METADATA_ENRICHMENT_QUERY_GENES, query_genes);
	ret.update(METADATA_ENRICHMENT_GENE_RATIO, gene_ratio);
	ret.update(METADATA_ENRICHMENT_BACKGROUND_RATIO, background_ratio);
	ret.update(METADATA_ENRICHMENT_COUNT, count);
	ret.row_slice(custom::less_than(p_adjusted, p_threshold));

	return ret;
};

QMap<QString, QMap<QString, QStringList>> get_pathway_information(soap::Species species) {

	QMap<QString, QMap<QString, QStringList>> pathway_information;

	{
		auto [pathway_to_ontology, pathway_to_size, pathway_to_symbol] = load_path_to_symbol("GO", species);
		auto pathway_to_name = load_pathway_name("GO", species);

		auto pathways = pathway_to_symbol.keys();
		for (const auto& pathway_name : pathways) {
			pathway_information["GO"][pathway_to_name[pathway_name]] = pathway_to_symbol[pathway_name];
		}
	}

	{
		auto [pathway_to_ontology, pathway_to_size, pathway_to_symbol] = load_path_to_symbol("KEGG", species);
		auto pathway_to_name = load_pathway_name("KEGG", species);

		auto pathways = pathway_to_symbol.keys();
		for (const auto& pathway_name : pathways) {
			pathway_information["KEGG"][pathway_to_name[pathway_name]] = pathway_to_symbol[pathway_name];
		}
	}

	return pathway_information;
};