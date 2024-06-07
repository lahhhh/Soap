#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>

class SimpleChoiceLayout : public QWidget
{
	Q_OBJECT
public:
	explicit SimpleChoiceLayout(const QString& name, const QStringList& item_names, QWidget* parent = nullptr);

	QString name_;

	QStringList item_names_;

	QVBoxLayout* main_layout_{ nullptr };
	QVBoxLayout* item_layout_{ nullptr };

	QPushButton* add_button_{ nullptr };
	QPushButton* expand_button_{ nullptr };

	QList<QHBoxLayout*> row_layouts_;
	QList<QComboBox*> combo_boxes_;
	QList<QPushButton*> buttons_;

	QString current_value();

signals:

private slots:

	void s_delete_choice();

	void s_add_choice();

	void s_expand();
};

