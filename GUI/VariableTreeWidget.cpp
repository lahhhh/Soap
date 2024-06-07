#include "VariableTreeWidget.h"

#include "VariableItem.h"

VariableTreeWidget::VariableTreeWidget(
	QWidget* parent, 
	PlotsSuite* draw_suite, 
	InformationTextBrowser* information_area, 
	SignalEmitter* signal_emitter) : 
	QTreeWidget(parent),
	draw_suite_(draw_suite), 
	information_area_(information_area),
	signal_emitter_(signal_emitter)
{
	connect(this, &QTreeWidget::itemDoubleClicked, this, &VariableTreeWidget::s_item_double_clicked);
	this->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(this, &QWidget::customContextMenuRequested, this, &VariableTreeWidget::s_item_right_clicked);

	QTreeWidgetItem* item = new VariableItem();
	this->addTopLevelItem(item);
}

void VariableTreeWidget::s_item_double_clicked(QTreeWidgetItem* tree, int) {

	if (this->currentItem() == nullptr)return;

	VariableItem* item = (VariableItem*)this->currentItem();

	item->__show_this();
}

void VariableTreeWidget::s_item_right_clicked(const QPoint&) {

	if (this->currentItem() == nullptr)return;
	
	VariableItem* item = (VariableItem*)this->currentItem();

	if(item->menus_.contains("Root"))
		item->menus_["Root"]->exec(QCursor::pos());
}
