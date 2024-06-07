#include "EnrichmentMultiGroupPlotDialog.h"

#include "Identifier.h"
#include <QScrollArea>

#include "SoapGUI.h"

EnrichmentMultiGroupPlotDialog::EnrichmentMultiGroupPlotDialog(
	const QMap<QString, QStringList>& group_map, 
	QWidget* parent
):
	QDialog(parent), 
	group_map_(group_map)
{
	this->main_layout_ = new QVBoxLayout;
	this->all_layout_ = new QVBoxLayout;

	QHBoxLayout* row_layout = new QHBoxLayout;

	EnrichmentGroupLayout* enrichment_layout = new EnrichmentGroupLayout("Pathway", this->group_map_, this);
	this->layout_groups_ << enrichment_layout;

	G_SET_NEW_BUTTON(button, "Delete", soap::MiddleSize)

	connect(button, &QPushButton::clicked, this, &EnrichmentMultiGroupPlotDialog::s_delete_group);
	this->delete_buttons_ << button;

	row_layout->addWidget(enrichment_layout);
	row_layout->addWidget(button);

	this->row_layouts_ << row_layout;

	this->group_layout_ = new QVBoxLayout;

	this->group_layout_->addLayout(row_layout);

	this->main_layout_->addLayout(this->group_layout_);

	G_SET_BUTTON(this->add_button_, "+", soap::MiddleSize)
	connect(this->add_button_, &QPushButton::clicked, this, &EnrichmentMultiGroupPlotDialog::s_add_group);

	this->main_layout_->addWidget(this->add_button_);
	this->main_layout_->setSizeConstraint(QLayout::SetFixedSize);

	QWidget* main_widget = new QWidget(this);
	main_widget->setLayout(this->main_layout_);
	main_widget->setObjectName("MainWidget");
	main_widget->setStyleSheet("#MainWidget{background-color:#f4f4ff;}");

	QScrollArea* scroll_area = new QScrollArea(this);
	scroll_area->setWidget(main_widget);
	scroll_area->setObjectName("ScrollArea");
	scroll_area->setStyleSheet("#ScrollArea{background-color:#f4f4ff;}");

	this->all_layout_->addWidget(scroll_area);

	G_SET_FINISH_BUTTON;
	G_SET_CANCEL_BUTTON;
	G_DOUBLE_ITEM_ROWLAYOUT(row_layout, this->finish_button_, this->cancel_button_);
	this->all_layout_->addLayout(row_layout);	

	G_SET_LABEL(this->capitalize_label_, "Capitalize First Character", soap::MiddleSize);
	G_SET_SWITCH(this->capitalize_switch_, true, this->capitalize_label_, soap::MiddleSize);
	G_DOUBLE_ITEM_ROWLAYOUT(row_layout, this->capitalize_label_, this->capitalize_switch_);
	this->all_layout_->addLayout(row_layout);

	this->setLayout(this->all_layout_);

	connect(this->finish_button_, &QPushButton::clicked, this, &EnrichmentMultiGroupPlotDialog::accept);
	connect(this->cancel_button_, &QPushButton::clicked, this, &EnrichmentMultiGroupPlotDialog::reject);

	G_SET_ICON;
	this->setWindowTitle("Enrichment Multi Group Plot Settings");
	this->exec();
};

void EnrichmentMultiGroupPlotDialog::s_add_group() {

	QHBoxLayout* row_layout = new QHBoxLayout;

	EnrichmentGroupLayout* enrichment_layout = new EnrichmentGroupLayout("Pathway", this->group_map_, this);
	this->layout_groups_ << enrichment_layout;
	G_SET_NEW_BUTTON(button, "Delete", soap::MiddleSize)
	connect(button, &QPushButton::clicked, this, &EnrichmentMultiGroupPlotDialog::s_delete_group);
	this->delete_buttons_ << button;

	row_layout->addWidget(enrichment_layout);
	row_layout->addWidget(button);

	this->row_layouts_ << row_layout;
	this->group_layout_->addLayout(row_layout);
};

void EnrichmentMultiGroupPlotDialog::s_delete_group() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
	int index = this->delete_buttons_.indexOf(button);

	this->group_layout_->removeItem(this->row_layouts_[index]);
	
	delete this->layout_groups_[index];
	delete this->delete_buttons_[index];
	delete this->row_layouts_[index];
	
	this->layout_groups_.remove(index);
	this->delete_buttons_.remove(index);
	this->row_layouts_.remove(index);
};

void EnrichmentMultiGroupPlotDialog::accept() {
	this->is_accepted_ = true;
	QDialog::accept();
};

void EnrichmentMultiGroupPlotDialog::reject() {
	QDialog::reject();
};

QPair<bool, QList<QStringList>> EnrichmentMultiGroupPlotDialog::get_plot_settings(const QMap<QString, QStringList>& group_map, QWidget* parent) {
	EnrichmentMultiGroupPlotDialog dlg(group_map, parent);
	QPair<bool, QList<QStringList>> ret;
	
	if (!dlg.is_accepted_) {
		return ret;
	}

	ret.first = dlg.capitalize_switch_->value_;

	for (auto ptr : dlg.layout_groups_) {
		ret.second << ptr->current_value();
	}

	return ret;
};