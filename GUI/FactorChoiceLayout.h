#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>

class FactorChoiceLayout : public QWidget
{
	Q_OBJECT
public:
	FactorChoiceLayout(
		const QString& name,
		const QMap<QString, QStringList>& factor_map_, 
		QWidget* parent
	);

	QString name_;

	QMap<QString, QStringList> factor_map_;

	QComboBox* factor_box_ = nullptr;

	QVBoxLayout* main_layout_ = nullptr;

	QVBoxLayout* item_layout_ = nullptr;

	QPushButton* add_button_ = nullptr;

	QPushButton* expand_button_ = nullptr;

	QList<QHBoxLayout*> row_layouts_;

	QList<QComboBox*> combo_boxes_;

	QList<QPushButton*> buttons_;

	QString current_value();

	void clear_choice();

private slots:

	void s_factor_changed();

	void s_expand_factor();

	void s_delete_choice();

	void s_add_choice();

};
