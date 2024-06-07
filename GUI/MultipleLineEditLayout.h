#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

#include "SignalEmitter.h"

class MultipleLineEditLayout : public QWidget
{
	Q_OBJECT
public:
	MultipleLineEditLayout(
		const QString& name,
		SignalEmitter* signal_emitter,
		const QStringList& item_names = {},
		QWidget* parent = nullptr
	);

	SignalEmitter* signal_emitter_{ nullptr };

	QString name_;

	QVBoxLayout* main_layout_{ nullptr };

	QVBoxLayout* item_layout_{ nullptr };

	QPushButton* add_button_{ nullptr };
	QPushButton* import_button_{ nullptr };

	QList<QHBoxLayout*> row_layouts_;

	QList<QLineEdit*> line_edits_;

	QList<QPushButton*> buttons_;

	QString current_value();

	void clear();
	void add_item(const QString& s);
	void add_items(const QStringList& items);
	void set_items(const QStringList& items);

private slots:

	void s_delete_line();

	void s_add_line();

	void s_import();

};
