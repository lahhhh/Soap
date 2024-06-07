#pragma once

#include "Identifier.h"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include "CustomMatrix.h"
#include "subLogicLayout.h"

#include "LogicHandler.h"

using CompleteLogic = QList<OrLogic>;

class LogicLayout : public QWidget
{
	Q_OBJECT
public:

	LogicLayout(const CustomMatrix* metadata, const QStringList* gene_names, QWidget* parent);

	LogicLayout(LogicHandler* lh, QWidget* parent);

	const CustomMatrix* metadata_;
	const QStringList* gene_names_;

	LogicHandler* logic_handler_{ nullptr };

	QHBoxLayout* full_layout_{ nullptr };
	QHBoxLayout* main_layout_{ nullptr };

	QList<SubLogicLayout* > or_layouts_;
	QList<QLabel* > and_labels_;
	QPushButton* and_button_;

	CompleteLogic get_filters();

	QString current_value();

private slots:

	void s_delete_sublogic_layout(SubLogicLayout* sll);

	void s_and_clicked();
};
