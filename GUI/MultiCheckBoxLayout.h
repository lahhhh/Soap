#pragma once

#include "Identifier.h"

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>

class MultiCheckBoxLayout :
    public QWidget
{
public:
    MultiCheckBoxLayout(
        const QStringList& item_names,
        QWidget* parent = nullptr
    );

    QString current_value();

private:

    QStringList item_names_;

    QVBoxLayout* main_layout_;

    QCheckBox* all_box_{ nullptr };

    QList<QCheckBox*> checkboxes_;

private slots:

    void s_set_all(int);

};

