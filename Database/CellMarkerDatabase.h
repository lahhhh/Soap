#pragma once

#include "Identifier.h"

#include <QMap>


class CellMarkerDatabase
{
public:

    CellMarkerDatabase() = default;

    CellMarkerDatabase(const CellMarkerDatabase&) = delete;

    QMap<QString, QMap<QString, QMap<QString, QStringList>>> static get_database(soap::Species species);
};
