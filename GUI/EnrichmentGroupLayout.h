#pragma once

#include "Identifier.h"

#include <QWidget>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLineEdit>

#include "FactorDoubleLineEditWithCompleterLayout.h"

class EnrichmentGroupLayout :
    public QWidget
{
public:
    EnrichmentGroupLayout(
        const QString& name, 
        const QMap<QString, QStringList>& map, 
        QWidget* parent = nullptr
    );

    QVBoxLayout* main_layout_;

    FactorDoubleLineEditWithCompleterLayout* select_layout_;

    QPushButton* legend_color_button_;

    QLabel* legend_color_label_;

    QStringList current_value();

private slots:

    void s_choose_color();
    
};

