#pragma once

#include "Identifier.h"

#include <QFile>
#include <QTextStream>
#include <iostream>
#include <fstream>
#include <cstdio>

#include "Serialization.h"

#include "windows.h"

class ItemDatabase
{
public:
	ItemDatabase() = default;

	ItemDatabase(const ItemDatabase&) = delete;

	template<typename ItemType>
	static bool read_item(const QString& item_name, ItemType& item) {

		std::ifstream ifs(item_name.toStdString(), std::ios::binary | std::ios::ate);

		if (!ifs.is_open()) {
			return false;
		}

		std::streamsize size = ifs.tellg();

		if (size < 16) {
			return false;
		}

		ifs.seekg(0, std::ios::beg);


		int id1{ 0 };
		int id2{ 0 };

		sread(ifs, id1);
		sread(ifs, id2);

		if (id1 != ITEM_IDENTIFIER_1 || id2 != ITEM_IDENTIFIER_2) {
			return false;
		}

		QString edition, name, type;

		sread(ifs, edition);

		if (edition != SOAP_EDITION) {
			return false;
		}

		sread(ifs, name);
		sread(ifs, type);

		if (type != ItemType::g_identifier()) {
			return false;
		}

		sread(ifs, item);

		return true;
	}

	template <typename ItemType>
	static bool write_item(const QString& item_name, const ItemType& item, const QString& file_name) {
		std::ofstream ofs(file_name.toStdString(), std::ios::binary);

		if (!ofs.is_open()) {
			return false;
		}

		int id1 = ITEM_IDENTIFIER_1;
		int id2 = ITEM_IDENTIFIER_2;
		QString edition = SOAP_EDITION;

		swrite(ofs, id1);
		swrite(ofs, id2);
		swrite(ofs, edition);
		swrite(ofs, item_name);
		swrite(ofs, ItemType::g_identifier());
		swrite(ofs, item);

		return true;
	}
};

