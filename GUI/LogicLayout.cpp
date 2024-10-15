#include "LogicLayout.h"
#include "Identifier.h"

#include "SoapGUI.h"

LogicLayout::LogicLayout(const CustomMatrix* metadata, const QStringList* gene_names, QWidget* parent) :
	QWidget(parent),
	metadata_(metadata),
	gene_names_(gene_names)
{
	this->full_layout_ = new QHBoxLayout;
	this->main_layout_ = new QHBoxLayout;
	this->or_layouts_ << new SubLogicLayout(this->metadata_, this->gene_names_, this);

	connect(this->or_layouts_[0], &SubLogicLayout::x_delete_this, this, &LogicLayout::s_delete_sublogic_layout);

	this->main_layout_->addWidget(this->or_layouts_[0]);
	this->full_layout_->addLayout(this->main_layout_);

	G_SET_BUTTON(this->and_button_, "and", QSize(50, 30));
	this->full_layout_->addWidget(this->and_button_);

	setLayout(this->full_layout_);

	connect(this->and_button_, &QPushButton::clicked, this, &LogicLayout::s_and_clicked);
}

LogicLayout::LogicLayout(LogicHandler* lh, QWidget* parent) :
	QWidget(parent),
	logic_handler_(lh)
{
	this->full_layout_ = new QHBoxLayout;
	this->main_layout_ = new QHBoxLayout;
	this->or_layouts_ << new SubLogicLayout(lh, this);

	connect(this->or_layouts_[0], &SubLogicLayout::x_delete_this, this, &LogicLayout::s_delete_sublogic_layout);

	this->main_layout_->addWidget(this->or_layouts_[0]);
	this->full_layout_->addLayout(this->main_layout_);

	G_SET_BUTTON(this->and_button_, "and", QSize(50, 30));
	this->full_layout_->addWidget(this->and_button_);

	setLayout(this->full_layout_);

	connect(this->and_button_, &QPushButton::clicked, this, &LogicLayout::s_and_clicked);
}

void LogicLayout::s_delete_sublogic_layout(SubLogicLayout* sll) {
	if (this->or_layouts_.size() == 1)return;

	int index = this->or_layouts_.indexOf(sll);

	if (index == 0) {
		this->full_layout_->removeWidget(this->and_labels_[0]);
		delete this->and_labels_[0];
		this->and_labels_.removeAt(0);
	}
	else {
		this->full_layout_->removeWidget(this->and_labels_[index - 1]);
		delete this->and_labels_[index - 1];
		this->and_labels_.removeAt(index - 1);
	}
	this->full_layout_->removeWidget(sll);
	delete sll;
	this->or_layouts_.removeAt(index);
};

void LogicLayout::s_and_clicked() {
	QLabel* andLabel;
	G_SET_LABEL(andLabel, "and", QSize(50, 30));
	this->and_labels_ << andLabel;
	this->or_layouts_ << new SubLogicLayout(this->logic_handler_, this);

	connect(this->or_layouts_[this->or_layouts_.size() - 1], &SubLogicLayout::x_delete_this, this, &LogicLayout::s_delete_sublogic_layout);

	this->main_layout_->addWidget(this->and_labels_[this->and_labels_.size() - 1]);
	this->main_layout_->addWidget(this->or_layouts_[this->or_layouts_.size() - 1]);
}

CompleteLogic LogicLayout::get_filters() {
	CompleteLogic ret;
	for (auto orLayout : this->or_layouts_) {
		ret << orLayout->get_filters();
	}
	return ret;
};

QString LogicLayout::current_value() {
	QStringList val;

	for (auto orLayout : this->or_layouts_) {
		val << orLayout->current_value();
	}

	return custom::merge_to_string(val, SOAP_DELIMITER3);
};
