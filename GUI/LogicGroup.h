#pragma once

#include "Identifier.h"

#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QCompleter>
#include <QPushButton>
#include "CustomMatrix.h"

#include "LogicHandler.h"

using SingleLogic = std::tuple<QString, QString, QString>;

class LogicGroup : public QWidget
{
	Q_OBJECT
public:
	LogicGroup(const CustomMatrix* metadata, const QStringList* gene_names, QWidget* parent);

	LogicGroup(LogicHandler* lh, QWidget* parent);

	LogicHandler* logic_handler_{ nullptr };

	const CustomMatrix* metadata_;
	const QStringList* gene_names_;

	QHBoxLayout* full_layout_;

	QLabel* feature_label_;
	QLineEdit* feature_name_;
	QCompleter* feature_completer_;
	QComboBox* compare_box_;
	QLineEdit* value_line_edit_;
	QComboBox* value_box_;
	QPushButton* delete_button_;

	bool factor_mode_ = false;

	void check_finished();

	SingleLogic get_filter();

	QString current_value();

	void to_factor_mode(const QStringList& levels);
	void to_numeric_mode();

signals:
	void x_delete_this(LogicGroup*);

private slots:

	void __s_delete_this();

	void s_check_feature();
};
