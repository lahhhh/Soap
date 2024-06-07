#pragma once

#include "Identifier.h"

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>
#include <QCompleter>

class MultipleDoubleLineEditWithCompleterLayout : public QWidget
{
    Q_OBJECT
public:
    MultipleDoubleLineEditWithCompleterLayout(
        const QString& name, 
        const QString& hint, 
        const QStringList& item_names, 
        QWidget* parent = nullptr
    );

    QString name_;

    QString hint_;

    QStringList item_names_;

    QCompleter* completer_;

    QVBoxLayout* main_layout_;

    QVBoxLayout* item_layout_;

    QPushButton* add_button_;

    QList<QHBoxLayout*> row_layouts_;

    QList<QLineEdit*> line_edits1_;

    QList<QLineEdit*> line_edits2_;

    QList<QPushButton*> buttons_;

    QString current_value();

signals:

private slots:
    void s_delete_line();

    void s_add_line();

};
