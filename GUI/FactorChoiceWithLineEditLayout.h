#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

class FactorChoiceWithLineEditLayout : public QWidget
{
    Q_OBJECT
public:
    explicit FactorChoiceWithLineEditLayout(const QString& name, const QString& hint, const QMap<QString, QStringList>& factor_map, QWidget* parent = nullptr);

    QString name_;

    QString hint_;

    QMap<QString, QStringList> factor_map_;

    QComboBox* factor_box_;

    QLineEdit* factor_line_edit_;

    QVBoxLayout* main_layout_;

    QVBoxLayout* item_layout_;

    QPushButton* add_button_;

    QList<QHBoxLayout*> row_layouts_;

    QList<QComboBox*> combo_boxes_;

    QList<QLineEdit*> line_edits_;

    QList<QPushButton*> buttons_;

    QString current_value();

private slots:

    void s_factor_changed();

    void s_delete_choice();

    void s_add_choice();

};
