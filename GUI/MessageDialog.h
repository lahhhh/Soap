#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>

class MessageDialog 
    : public QDialog
{

public:

    MessageDialog(const QString& title, const QString& message);

    void static get_response(const QString& title, const QString& message);

};

