#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

class SimpleChoiceWithLineEditLayout : public QWidget
{
    Q_OBJECT
public:
    SimpleChoiceWithLineEditLayout(
        const QString& name, 
        const QString& hint, 
        const QStringList& item_names, 
        QWidget* parent = nullptr
    );

    QString name_;

    QString hint_;

    QStringList item_names_;

    QVBoxLayout* main_layout_;

    QVBoxLayout* item_layout_;

    QPushButton* add_button_;

    QList<QHBoxLayout*> row_layouts_;

    QList<QComboBox*> combo_boxes_;

    QList<QLineEdit*> line_edits_;

    QList<QPushButton*> buttons_;

    QString current_value();

signals:

private slots:

    void s_delete_choice();

    void s_add_choice();
};
