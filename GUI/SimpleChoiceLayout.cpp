#include "SimpleChoiceLayout.h"

#include "SoapGUI.h"

SimpleChoiceLayout::SimpleChoiceLayout(const QString& name, const QStringList& item_names, QWidget* parent) :
	QWidget(parent),
	name_(name),
	item_names_(item_names)
{
	this->main_layout_ = new QVBoxLayout;
	setLayout(this->main_layout_);

	this->item_layout_ = new QVBoxLayout;

	G_SET_NEW_COMBOBOX(box, this->item_names_, 30);
	this->combo_boxes_ << box;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;

	connect(button, &QPushButton::clicked, this, &SimpleChoiceLayout::s_delete_choice);

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, button);

	this->row_layouts_ << row_layout;

	this->main_layout_->addLayout(this->item_layout_);

	G_SET_BUTTON(this->add_button_, "+ " + this->name_, soap::MiddleSize);
	connect(this->add_button_, &QPushButton::clicked, this, &SimpleChoiceLayout::s_add_choice);

	G_SET_BUTTON(this->expand_button_, "++", QSize(40, 30));
	connect(this->expand_button_, &QPushButton::clicked, this, &SimpleChoiceLayout::s_expand);

	G_ADD_DOUBLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_, this->expand_button_);
}

void SimpleChoiceLayout::s_expand() {

	int n_layout = this->row_layouts_.size();

	for (int i = 0; i < n_layout; ++i) {

		this->item_layout_->removeItem(this->row_layouts_[i]);

		delete this->combo_boxes_[i];
		delete this->buttons_[i];
		delete this->row_layouts_[i];

	}

	this->combo_boxes_.clear();
	this->buttons_.clear();
	this->row_layouts_.clear();

	int n_item = this->item_names_.size();

	for (int i = 0; i < n_item; ++i) {
		G_SET_NEW_COMBOBOX(box, this->item_names_, 30);
		this->combo_boxes_ << box;

		box->setCurrentIndex(i);

		G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
		this->buttons_ << button;
		connect(button, &QPushButton::clicked, this, &SimpleChoiceLayout::s_delete_choice);

		G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, button);
		this->row_layouts_ << row_layout;
	}
};

void SimpleChoiceLayout::s_delete_choice() {
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

void SimpleChoiceLayout::s_add_choice() {
	G_SET_NEW_COMBOBOX(box, this->item_names_, 30);
	this->combo_boxes_ << box;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;
	connect(button, &QPushButton::clicked, this, &SimpleChoiceLayout::s_delete_choice);

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, button);
	this->row_layouts_ << row_layout;
};

QString SimpleChoiceLayout::current_value() {
	QString ret;
	for (const auto ptr : this->combo_boxes_) {
		ret += SOAP_DELIMITER + ptr->currentText();
	}
	return ret;
};
