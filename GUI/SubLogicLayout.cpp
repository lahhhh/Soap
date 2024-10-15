#include "subLogicLayout.h"

#include "SoapGUI.h"

SubLogicLayout::SubLogicLayout(const CustomMatrix* metadata, const QStringList* gene_names, QWidget* parent) :
	QWidget(parent), 
	metadata_(metadata), 
	gene_names_(gene_names)
{
	this->full_layout_ = new QVBoxLayout;
	this->main_layout_ = new QVBoxLayout;

	this->logic_groups_ << new LogicGroup(metadata, gene_names, this);

	connect(this->logic_groups_[0], &LogicGroup::x_delete_this, this, &SubLogicLayout::s_delete_logic_group);

	G_SET_BUTTON(this->or_button_, "or", QSize(50, 30));

	this->main_layout_->addWidget(this->logic_groups_[0]);

	this->full_layout_->addLayout(this->main_layout_);

	G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT(this->full_layout_, row_layout, this->or_button_);

	setLayout(this->full_layout_);
	connect(this->or_button_, &QPushButton::clicked, this, &SubLogicLayout::s_or_clicked);
}

SubLogicLayout::SubLogicLayout(LogicHandler* lh, QWidget* parent)
	: QWidget(parent), logic_handler_(lh)
{
	this->full_layout_ = new QVBoxLayout;
	this->main_layout_ = new QVBoxLayout;

	this->logic_groups_ << new LogicGroup(lh, this);

	connect(this->logic_groups_[0], &LogicGroup::x_delete_this, this, &SubLogicLayout::s_delete_logic_group);

	G_SET_BUTTON(this->or_button_, "or", QSize(50, 30));

	this->main_layout_->addWidget(this->logic_groups_[0]);

	this->full_layout_->addLayout(this->main_layout_);

	G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT(this->full_layout_, row_layout, this->or_button_);

	setLayout(this->full_layout_);
	connect(this->or_button_, &QPushButton::clicked, this, &SubLogicLayout::s_or_clicked);
}

void SubLogicLayout::s_or_clicked() {
	G_SET_NEW_LABEL(or_label, "or", QSize(50, 30));
	this->or_labels_ << or_label;

	this->logic_groups_ << new LogicGroup(this->logic_handler_, this);
	connect(this->logic_groups_[this->logic_groups_.size() - 1], &LogicGroup::x_delete_this, this, &SubLogicLayout::s_delete_logic_group);

	G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, or_label);
	this->main_layout_->addWidget(this->logic_groups_[this->logic_groups_.size() - 1]);
}

void SubLogicLayout::s_delete_logic_group(LogicGroup* lg) {
	if (this->logic_groups_.size() == 1) {
		emit x_delete_this(this);
	}
	else {
		int index = this->logic_groups_.indexOf(lg);
		this->full_layout_->removeWidget(lg);

		if (index == 0) {
			this->full_layout_->removeWidget(this->or_labels_[0]);
			delete this->or_labels_[0];
			this->or_labels_.removeAt(0);
		}
		else {
			this->full_layout_->removeWidget(this->or_labels_[index - 1]);
			delete this->or_labels_[index - 1];
			this->or_labels_.removeAt(index - 1);
		}
		delete lg;
		this->logic_groups_.removeAt(index);
	}
};

OrLogic SubLogicLayout::get_filters() {
	OrLogic ret;
	for (auto lg : this->logic_groups_) {
		ret << lg->get_filter();
	}
	return ret;
};

QString SubLogicLayout::current_value() {

	QStringList val;

	for (auto lg : this->logic_groups_) {
		val << lg->current_value();
	}

	return custom::merge_to_string(val, SOAP_DELIMITER2);
};
