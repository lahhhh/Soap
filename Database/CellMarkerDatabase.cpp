#include "CellMarkerDatabase.h"

#include <QFile>
#include <QTextStream>

QMap<QString, QMap<QString, QStringList>> load_cell_marker_mouse_database() {

	QMap<QString, QMap<QString, QStringList>> db;

	QFile file(FILE_MOUSE_CELLMARKER);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index1 = line.indexOf('$'), index2 = line.indexOf(':');

		QString tissue = line.sliced(0, index1);
		QString cell_type = line.sliced(index1 + 1, index2 - index1 - 1);
		QString markers = line.sliced(index2 + 1, line.size() - index2 - 1);

		db[tissue][cell_type] = markers.split(';');

		line = in.readLine();
	}

	return db;
};

QMap<QString, QMap<QString, QStringList>> load_cell_marker_human_database() {

	QMap<QString, QMap<QString, QStringList>> db;

	QFile file(FILE_HUMAN_CELLMARKER);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QString line = in.readLine();
	while (!line.isNull()) {
		qsizetype index1 = line.indexOf('$'), index2 = line.indexOf(':');

		QString tissue = line.sliced(0, index1);
		QString cell_type = line.sliced(index1 + 1, index2 - index1 - 1);
		QString markers = line.sliced(index2 + 1, line.size() - index2 - 1);

		db[tissue][cell_type] = markers.split(';');

		line = in.readLine();
	}

	return db;
};

QMap<QString, QMap<QString, QStringList>> load_cell_marker_database(soap::Species species) {

	if (species == soap::Species::Human) {
		return load_cell_marker_human_database();
	}
	else {
		return load_cell_marker_mouse_database();
	}
};

QMap<QString, QMap<QString, QStringList>> load_panglaodb_mouse_database() {

	QMap<QString, QMap<QString, QStringList>> db;

	QFile file(FILE_MOUSE_PANGLAODB);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index = line.indexOf(':');

		QString cell_type = line.sliced(0, index);
		QString markers = line.sliced(index + 1, line.size() - index - 1);

		db["All"][cell_type] = markers.split(';');

		line = in.readLine();
	}

	return db;
};

QMap<QString, QMap<QString, QStringList>> load_panglaodb_human_database() {

	QMap<QString, QMap<QString, QStringList>> db;

	QFile file(FILE_HUMAN_PANGLAODB);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QString line = in.readLine();

	while (!line.isNull()) {
		qsizetype index = line.indexOf(':');

		QString cell_type = line.sliced(0, index);
		QString markers = line.sliced(index + 1, line.size() - index - 1);

		db["All"][cell_type] = markers.split(';');

		line = in.readLine();
	}

	return db;
};

QMap<QString, QMap<QString, QStringList>> load_panglaodb_database(soap::Species species) {

	if (species == soap::Species::Human) {
		return load_panglaodb_human_database();
	}
	else if(species == soap::Species::Mouse){
		return load_panglaodb_mouse_database();
	}
};

QMap<QString, QMap<QString, QMap<QString, QStringList>>> CellMarkerDatabase::get_database(soap::Species species) {

	QMap<QString, QMap<QString, QMap<QString, QStringList>>> db;

	db["CellMarker"] = load_cell_marker_database(species);
	db["PanglaoDB"] = load_panglaodb_database(species);

	return db;
};
