#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

class MultiLineEditWithColorChoiceLayout : public QWidget
{
    Q_OBJECT
public:
    MultiLineEditWithColorChoiceLayout(
        const QString& name,
        QWidget* parent = nullptr
    );

    QString name_;

    QVBoxLayout* main_layout_;

    QVBoxLayout* item_layout_;

    QPushButton* add_button_;

    QList<QHBoxLayout*> row_layouts_;

    QList<QLineEdit*> line_edits_;

    QList<QPushButton*> color_buttons_;

    QList<QLabel*> color_labels_;

    QList<QPushButton*> buttons_;

    QString current_value();

private slots:

    void s_choose_color();

    void s_delete_choice();

    void s_add_choice();
};

