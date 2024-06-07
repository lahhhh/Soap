#include "EnrichmentDatabase.h"

#include <QFile>
#include <QTextStream>

#include "Phyper.h"
#include "Custom.h"


QMap<QString, QMap<QString, QStringList>> EnrichmentDatabase::get_pathway_information(soap::Species species) {
	QMap<QString, QMap<QString, QStringList>> pathway_information;

	if (species == soap::Species::Human) {
		this->load_go_human();
		this->load_kegg_human();

		const auto& pathway_to_symbol_kegg = this->pathway_to_symbol_["KEGG_HUMAN"];
		const auto& pathway_to_name_kegg = this->pathway_to_pathway_name_["KEGG_HUMAN"];

		for (const auto& pathway_name : pathway_to_symbol_kegg.keys()) {
			pathway_information["KEGG_HUMAN"][pathway_to_name_kegg[pathway_name]] = pathway_to_symbol_kegg[pathway_name];
		}

		const auto& pathway_to_symbol_go = this->pathway_to_symbol_["GO_HUMAN"];
		const auto& pathway_to_name_go = this->pathway_to_pathway_name_["GO"];

		for (const auto& pathway_name : pathway_to_symbol_go.keys()) {
			pathway_information["GO"][pathway_to_name_go[pathway_name]] = pathway_to_symbol_go[pathway_name];
		}

		return pathway_information;
	}
	else {
		load_go_mouse();
		load_kegg_mouse();

		const auto& pathway_to_symbol_kegg = this->pathway_to_symbol_["KEGG_HUMAN"];
		const auto& pathway_to_name_kegg = this->pathway_to_pathway_name_["KEGG_HUMAN"];

		for (const auto& pathway_name : pathway_to_symbol_kegg.keys()) {
			pathway_information["KEGG_HUMAN"][pathway_to_name_kegg[pathway_name]] = pathway_to_symbol_kegg[pathway_name];
		}

		const auto& pathway_to_symbol_go = this->pathway_to_symbol_["GO_MOUSE"];
		const auto& pathway_to_name_go = this->pathway_to_pathway_name_["GO"];

		for (const auto& pathway_name : pathway_to_symbol_go.keys()) {
			pathway_information["GO"][pathway_to_name_go[pathway_name]] = pathway_to_symbol_go[pathway_name];
		}

		return pathway_information;
	}
};

void EnrichmentDatabase::load_go_pathway_name() {

	QFile file(FILE_GO_PATH2NAME);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QString>& pathway_to_name = this->pathway_to_pathway_name_["GO"];
	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index_of_space = line.indexOf(' ');
		pathway_to_name[line.sliced(0, index_of_space)] = line.sliced(index_of_space + 1);

		line = in.readLine();
	}
};

void EnrichmentDatabase::load_kegg_pathway_name_human() {

	QFile file(FILE_KEGG_PATH2NAME_HUMAN);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QString>& pathway_to_name = this->pathway_to_pathway_name_["KEGG_HUMAN"];
	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index_of_space = line.indexOf(' ');
		pathway_to_name[line.sliced(0, index_of_space)] = line.sliced(index_of_space + 1);
		line = in.readLine();
	}
};

void EnrichmentDatabase::load_kegg_pathway_name_mouse() {

	QFile file(FILE_KEGG_PATH2NAME_MOUSE);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QMap<QString, QString>& pathway_to_name = this->pathway_to_pathway_name_["KEGG_MOUSE"];
	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index_of_space = line.indexOf(' ');
		pathway_to_name[line.sliced(0, index_of_space)] = line.sliced(index_of_space + 1);
		line = in.readLine();
	}

};

void EnrichmentDatabase::load_go_path_to_symbol_human() {

	QFile file(FILE_GO_PATH2SYMBOL_HUMAN);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QString>& pathway_to_ontology = this->pathway_to_ontology_["GO_HUMAN"];
	QMap<QString, int>& pathway_to_size = this->pathway_to_size_["GO_HUMAN"];
	QMap<QString, QStringList>& pathway_to_symbol = this->pathway_to_symbol_["GO_HUMAN"];

	QString line = in.readLine();

	while (!line.isNull()) {

		QStringList tmp = line.split(' ');
		QString path = tmp[0];

		pathway_to_ontology[path] = tmp[1];
		pathway_to_symbol[path] = tmp.sliced(2);
		pathway_to_size[path] = tmp.size() - 2;

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_kegg_path_to_symbol_human() {

	QFile file(FILE_KEGG_PATH2SYMBOL_HUMAN);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, int>& pathway_to_size = this->pathway_to_size_["KEGG_HUMAN"];
	QMap<QString, QStringList>& pathway_to_symbol = this->pathway_to_symbol_["KEGG_HUMAN"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');
		QString path = tmp[0];

		pathway_to_symbol[path] = tmp.sliced(1);
		pathway_to_size[path] = tmp.size() - 1;

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_go_symbol_to_path_human() {

	QFile file(FILE_GO_SYMBOL2PATH_HUMAN);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QStringList>& symbol_to_pathway = this->symbol_to_pathway_["GO_HUMAN"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');

		symbol_to_pathway[tmp[0]] = tmp.sliced(1);

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_kegg_symbol_to_path_human() {

	QFile file(FILE_KEGG_SYMBOL2PATH_HUMAN);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);
	QMap<QString, QStringList>& symbol_to_pathway = this->symbol_to_pathway_["KEGG_HUMAN"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');

		symbol_to_pathway[tmp[0]] = tmp.sliced(1);

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_go_path_to_symbol_mouse() {

	QFile file(FILE_GO_PATH2SYMBOL_MOUSE);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QString>& pathway_to_ontology = this->pathway_to_ontology_["GO_MOUSE"];
	QMap<QString, int>& pathway_to_size = this->pathway_to_size_["GO_MOUSE"];
	QMap<QString, QStringList>& pathway_to_symbol = this->pathway_to_symbol_["GO_MOUSE"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');

		QString path = tmp[0];

		pathway_to_ontology[path] = tmp[1];
		pathway_to_symbol[path] = tmp.sliced(2);
		pathway_to_size[path] = tmp.size() - 2;

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_kegg_path_to_symbol_mouse() {

	QFile file(FILE_KEGG_PATH2SYMBOL_MOUSE);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, int>& pathway_to_size = this->pathway_to_size_["KEGG_MOUSE"];
	QMap<QString, QStringList>& pathway_to_symbol = this->pathway_to_symbol_["KEGG_MOUSE"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');
		QString path = tmp[0];

		pathway_to_symbol[path] = tmp.sliced(1);
		pathway_to_size[path] = tmp.size() - 1;

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_go_symbol_to_path_mouse() {

	QFile file(FILE_GO_SYMBOL2PATH_MOUSE);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QStringList>& symbol_to_pathway = this->symbol_to_pathway_["GO_MOUSE"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');

		symbol_to_pathway[tmp[0]] = tmp.sliced(1);

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_kegg_symbol_to_path_mouse() {

	QFile file(FILE_KEGG_SYMBOL2PATH_MOUSE);
	file.open(QIODevice::ReadOnly | QIODevice::Text);

	QTextStream in(&file);

	QMap<QString, QStringList>& symbol_to_pathway = this->symbol_to_pathway_["KEGG_MOUSE"];

	QString line = in.readLine();

	while (!line.isNull()) {
		QStringList tmp = line.split(' ');

		symbol_to_pathway[tmp[0]] = tmp.sliced(1);

		line = in.readLine();
	}

};

void EnrichmentDatabase::load_go_human() {
	this->load_go_path_to_symbol_human();
	this->load_go_pathway_name();
	this->load_go_symbol_to_path_human();
};

void EnrichmentDatabase::load_go_mouse() {
	this->load_go_path_to_symbol_mouse();
	this->load_go_pathway_name();
	this->load_go_symbol_to_path_mouse();
};

void EnrichmentDatabase::load_kegg_human() {
	this->load_kegg_path_to_symbol_human();
	this->load_kegg_pathway_name_human();
	this->load_kegg_symbol_to_path_human();
};

void EnrichmentDatabase::load_kegg_mouse() {
	this->load_kegg_path_to_symbol_mouse();
	this->load_kegg_pathway_name_mouse();
	this->load_kegg_symbol_to_path_mouse();
};

CustomMatrix EnrichmentDatabase::enrich_go(
	const QStringList& gene_names, 
	const QString& ontology, 
	const soap::Species& species, 
	const QString& adjust_p_value_method, 
	double p_threshold
) {
	int minimum_pathway_size = 10, maximum_pathway_size = 500;
	QStringList ret_path, path_id, query_genes, pathway_genes;
	QVector<int> count;
	QMap<QString, QString> query_map;
	QStringList gene_ratio, background_ratio;
	QVector<double> p_value, p_adjusted;

	QMap<QString, QStringList>* symbol_to_pathway, * pathway_to_symbol;
	QMap<QString, QString>* pathway_to_ontology, * pathway_to_name;
	QMap<QString, int>* pathway_to_size;

	if (species == soap::Species::Human) {
		this->load_go_human();
		symbol_to_pathway = &this->symbol_to_pathway_["GO_HUMAN"];
		pathway_to_symbol = &this->pathway_to_symbol_["GO_HUMAN"];
		pathway_to_ontology = &this->pathway_to_ontology_["GO_HUMAN"];
		pathway_to_size = &this->pathway_to_size_["GO_HUMAN"];
	}
	else {
		this->load_go_mouse();
		symbol_to_pathway = &this->symbol_to_pathway_["GO_MOUSE"];
		pathway_to_symbol = &this->pathway_to_symbol_["GO_MOUSE"];
		pathway_to_ontology = &this->pathway_to_ontology_["GO_MOUSE"];
		pathway_to_size = &this->pathway_to_size_["GO_MOUSE"];
	}

	pathway_to_name = &this->pathway_to_pathway_name_["GO"];

	QStringList filtered_genes = _Cs intersect(gene_names, symbol_to_pathway->keys());

	qsizetype filtered_size = filtered_genes.size();

	// no enough gene found, return empty matrix
	if (filtered_size < 10) {
		return CustomMatrix();
	}

	QStringList path_names;

	for (const auto& gene_name : filtered_genes) {
		path_names << symbol_to_pathway->operator[](gene_name);
		auto& paths = symbol_to_pathway->operator[](gene_name);
		for (QString& path_name : paths) {
			query_map[path_name] += (gene_name + ";");
		}
	}


	QMap<QString, int> query;
	// calculate query times for pathways
	for (auto& path_name : path_names) {
		++query[path_name];
	}

	QStringList unique_paths = _Cs unique(path_names);
	QStringList filtered_paths;

	if (ontology == "ALL") {
		filtered_paths = unique_paths;
	}
	else {
		for (auto& path_name : unique_paths) {
			if (pathway_to_ontology->operator[](path_name) == ontology) {
				filtered_paths << path_name;
			}
		}
	}


	double database_size = symbol_to_pathway->size();

	for (auto& path_name : filtered_paths) {

		// check pathway size
		if (pathway_to_size->operator[](path_name) < minimum_pathway_size ||
			pathway_to_size->operator[](path_name) > maximum_pathway_size)
			continue;

		double p = phyper(
			query[path_name] - 1.0, 
			(double)pathway_to_size->operator[](path_name), 
			database_size - pathway_to_size->operator[](path_name), 
			(double)filtered_size
		);

		path_id << path_name;
		ret_path << pathway_to_name->operator[](path_name);
		p_value << p;
		query_genes << query_map[path_name];
		pathway_genes << pathway_to_symbol->operator[](path_name).join(",");
		count << query[path_name];
		gene_ratio << QString::number(query[path_name]) + "/" + QString::number(filtered_size);
		background_ratio << QString::number(pathway_to_size->operator[](path_name)) + "/" + QString::number(database_size);
	}

	p_adjusted = _Cs adjust_p_value(p_value, adjust_p_value_method);
	auto index = _Cs order(p_adjusted);
	p_adjusted = _Cs reordered(p_adjusted, index);
	path_id = _Cs reordered(path_id, index);
	ret_path = _Cs reordered(ret_path, index);
	p_value = _Cs reordered(p_value, index);
	query_genes = _Cs reordered(query_genes, index);
	pathway_genes = _Cs reordered(pathway_genes, index);
	gene_ratio = _Cs reordered(gene_ratio, index);
	background_ratio = _Cs reordered(background_ratio, index);
	count = _Cs reordered(count, index);

	CustomMatrix ret(path_id);

	ret.update(METADATA_ENRICHMENT_PATHWAY_NAMES, ret_path);
	ret.update(METADATA_ENRICHMENT_P_VALUE, p_value);
	ret.update(METADATA_ENRICHMENT_ADJUSTED_P_VALUE, p_adjusted);
	ret.update(METADATA_ENRICHMENT_QUERY_GENES, query_genes);
	ret.update(METADATA_ENRICHMENT_GENE_RATIO, gene_ratio);
	ret.update(METADATA_ENRICHMENT_BACKGROUND_RATIO, background_ratio);
	ret.update(METADATA_ENRICHMENT_COUNT, count);
	ret.row_slice(_Cs less_than(p_adjusted, p_threshold));

	return ret;
};

CustomMatrix EnrichmentDatabase::enrich_kegg(
	const QStringList& gene_names, 
	const soap::Species& species, 
	const QString& adjust_p_valueMethod, 
	double pThreshold
) {
	int minimum_pathway_size = 10, maximum_pathway_size = 500;
	QStringList ret_path, path_id, query_genes, pathway_genes;
	QVector<int> count;
	QMap<QString, QString> query_map;
	QStringList gene_ratio, background_ratio;
	QVector<double> p_value, p_adjusted;

	QMap<QString, QStringList>* symbol_to_pathway, * pathway_to_symbol;
	QMap<QString, QString>* pathway_to_name;
	QMap<QString, int>* pathway_to_size;

	if (species == soap::Species::Human) {
		load_kegg_human();
		symbol_to_pathway = &this->symbol_to_pathway_["KEGG_HUMAN"];
		pathway_to_symbol = &this->pathway_to_symbol_["KEGG_HUMAN"];
		pathway_to_size = &this->pathway_to_size_["KEGG_HUMAN"];
		pathway_to_name = &this->pathway_to_pathway_name_["KEGG_HUMAN"];
	}
	else {
		load_kegg_mouse();
		symbol_to_pathway = &this->symbol_to_pathway_["KEGG_MOUSE"];
		pathway_to_symbol = &this->pathway_to_symbol_["KEGG_MOUSE"];
		pathway_to_size = &this->pathway_to_size_["KEGG_MOUSE"];
		pathway_to_name = &this->pathway_to_pathway_name_["KEGG_MOUSE"];
	}

	QStringList filtered_genes = _Cs intersect(gene_names, symbol_to_pathway->keys());

	qsizetype filtered_size = filtered_genes.size();
	// not enough query gene
	if (filtered_size < 10) {
		return CustomMatrix();
	}

	QStringList path_names;
	for (auto& gene_name : filtered_genes) {
		path_names << symbol_to_pathway->operator[](gene_name);
		auto& paths = symbol_to_pathway->operator[](gene_name);
		for (QString& path_name : paths) {
			query_map[path_name] += (gene_name + ";");
		}
	}

	QStringList unique_paths = _Cs unique(path_names);

	QMap<QString, int> query;
	for (auto& path_name : path_names) {
		++query[path_name];
	}

	double database_size = symbol_to_pathway->size();

	for (auto& path_name : unique_paths) {

		//check pathway size
		if (pathway_to_size->operator[](path_name) < minimum_pathway_size ||
			pathway_to_size->operator[](path_name) > maximum_pathway_size)
			continue;

		double p = phyper(
			query[path_name] - 1.0, 
			(double)pathway_to_size->operator[](path_name), 
			database_size - pathway_to_size->operator[](path_name), 
			(double)filtered_size
		);

		path_id << path_name;
		ret_path << pathway_to_name->operator[](path_name);
		p_value << p;
		query_genes << query_map[path_name];
		pathway_genes << pathway_to_symbol->operator[](path_name).join(",");
		count << query[path_name];
		gene_ratio << QString::number(query[path_name]) + "/" + QString::number(filtered_size);
		background_ratio << QString::number(pathway_to_size->operator[](path_name)) + "/" + QString::number(database_size);
	}

	p_adjusted = _Cs adjust_p_value(p_value, adjust_p_valueMethod);

	auto index = _Cs order(p_adjusted);

	p_adjusted = _Cs reordered(p_adjusted, index);
	path_id = _Cs reordered(path_id, index);
	ret_path = _Cs reordered(ret_path, index);
	p_value = _Cs reordered(p_value, index);
	query_genes = _Cs reordered(query_genes, index);
	pathway_genes = _Cs reordered(pathway_genes, index);
	gene_ratio = _Cs reordered(gene_ratio, index);
	background_ratio = _Cs reordered(background_ratio, index);
	count = _Cs reordered(count, index);

	CustomMatrix ret(path_id);
	ret.update(METADATA_ENRICHMENT_PATHWAY_NAMES, ret_path);
	ret.update(METADATA_ENRICHMENT_P_VALUE, p_value);
	ret.update(METADATA_ENRICHMENT_ADJUSTED_P_VALUE, p_adjusted);
	ret.update(METADATA_ENRICHMENT_QUERY_GENES, query_genes);
	ret.update(METADATA_ENRICHMENT_GENE_RATIO, gene_ratio);
	ret.update(METADATA_ENRICHMENT_BACKGROUND_RATIO, background_ratio);
	ret.update(METADATA_ENRICHMENT_COUNT, count);
	ret.row_slice(_Cs less_than(p_adjusted, pThreshold));

	return ret;
};
