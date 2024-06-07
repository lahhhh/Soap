#include "FactorChoiceLayout.h"

#include "Identifier.h"

#include "SoapGUI.h"

FactorChoiceLayout::FactorChoiceLayout(
	const QString& name, 
	const QMap<QString, QStringList>& factor_map,
	QWidget* parent) :
	QWidget(parent),
	name_(name),
	factor_map_(factor_map)
{
	this->main_layout_ = new QVBoxLayout;
	this->setLayout(this->main_layout_);

	G_SET_COMBOBOX(this->factor_box_, this->factor_map_.keys(), 30);
	connect(this->factor_box_, &QComboBox::currentIndexChanged, this, &FactorChoiceLayout::s_factor_changed);
	G_SET_BUTTON(this->expand_button_, "+", QSize(30, 30));
	connect(this->expand_button_, &QPushButton::clicked, this, &FactorChoiceLayout::s_expand_factor);
	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, row_layout, this->factor_box_, this->expand_button_);

	this->item_layout_ = new QVBoxLayout;
	this->main_layout_->addLayout(this->item_layout_);

	G_SET_BUTTON(this->add_button_, "+ " + this->name_, soap::MiddleSize);
	connect(this->add_button_, &QPushButton::clicked, this, &FactorChoiceLayout::s_add_choice);
	G_ADD_SINGLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_);
}

void FactorChoiceLayout::s_expand_factor() {

	this->clear_choice();

	QString factor = this->factor_box_->currentText();
	const auto& levels = this->factor_map_[factor];

	for (int i = 0; i < levels.size() && i < 50; ++i) {
		G_SET_NEW_COMBOBOX(box, levels, 30);
		box->setCurrentIndex(i);
		this->combo_boxes_ << box;

		G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
		this->buttons_ << button;
		connect(button, &QPushButton::clicked, this, &FactorChoiceLayout::s_delete_choice);

		G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, button);
		this->row_layouts_ << row_layout;
	}
};

void FactorChoiceLayout::clear_choice() {

	const qsizetype size = this->buttons_.size();

	if (size < 1)return;

	for (qsizetype i = 0; i < size; ++i) {		

		this->item_layout_->removeItem(this->row_layouts_[i]);

		delete this->combo_boxes_[i];
		delete this->buttons_[i];
		delete this->row_layouts_[i];
	}

	this->combo_boxes_.clear();
	this->buttons_.clear();
	this->row_layouts_.clear();
};

void FactorChoiceLayout::s_factor_changed() {

	this->clear_choice();
};

void FactorChoiceLayout::s_delete_choice() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
	int index = this->buttons_.indexOf(button);

	this->item_layout_->removeItem(this->row_layouts_[index]);

	delete this->combo_boxes_[index];
	delete this->buttons_[index];
	delete this->row_layouts_[index];

	this->combo_boxes_.remove(index);
	this->buttons_.remove(index);
	this->row_layouts_.remove(index);

};

void FactorChoiceLayout::s_add_choice() {

	G_SET_NEW_COMBOBOX(box, this->factor_map_[this->factor_box_->currentText()], 30);
	this->combo_boxes_ << box;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;
	connect(button, &QPushButton::clicked, this, &FactorChoiceLayout::s_delete_choice);

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, button);
	this->row_layouts_ << row_layout;
};

QString FactorChoiceLayout::current_value() {
	QString ret = this->factor_box_->currentText();

	for (const auto ptr : this->combo_boxes_) {
		ret += SOAP_DELIMITER + ptr->currentText();
	}

	return ret;
};
