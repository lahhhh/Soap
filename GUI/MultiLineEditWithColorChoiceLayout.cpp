#include "MultiLineEditWithColorChoiceLayout.h"

#include <QColorDialog>
#include "Identifier.h"
#include "Custom.h"

#include "SoapGUI.h"

MultiLineEditWithColorChoiceLayout::MultiLineEditWithColorChoiceLayout(
    const QString& name,
    QWidget* parent
):
    name_(name)
{

	this->main_layout_ = new QVBoxLayout;
	setLayout(this->main_layout_);

	this->item_layout_ = new QVBoxLayout;

	QLineEdit* line = new QLineEdit(this);
	line->setFixedSize(120, 30);
	this->line_edits_ << line;

	G_SET_NEW_BUTTON(color_button, "Choose Color", soap::MiddleSize);
		this->color_buttons_ << color_button;

	connect(color_button, &QPushButton::clicked, this, &MultiLineEditWithColorChoiceLayout::s_choose_color);

	QColor color = custom::random_color();
	QLabel* color_label = new QLabel(color.name(), this);
	color_label->setStyleSheet("QLabel{color:" + color.name() + "}");
	color_label->setFixedSize(soap::MiddleSize);
	this->color_labels_ << color_label;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;

	connect(button, &QPushButton::clicked, this, &MultiLineEditWithColorChoiceLayout::s_delete_choice);

	G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line, color_button, color_label, button);
	this->row_layouts_ << row_layout;

	this->main_layout_->addLayout(this->item_layout_);

	G_SET_BUTTON(this->add_button_, "+ " + this->name_, soap::MiddleSize);
	connect(this->add_button_, &QPushButton::clicked, this, &MultiLineEditWithColorChoiceLayout::s_add_choice);

	G_ADD_SINGLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_);
};

void MultiLineEditWithColorChoiceLayout::s_choose_color() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
	int index = this->color_buttons_.indexOf(button);

	QColor color = QColorDialog::getColor(QColor(this->color_labels_[index]->text()), this, "Choose Color for Legend");

	if (color.isValid()) {
		this->color_labels_[index]->setText(color.name());
		this->color_labels_[index]->setStyleSheet("QLabel{color:" + color.name() + "}");
	}
};

void MultiLineEditWithColorChoiceLayout::s_delete_choice() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
	int index = this->buttons_.indexOf(button);
	this->item_layout_->removeItem(this->row_layouts_[index]);

	delete this->line_edits_[index];
	delete this->color_buttons_[index];
	delete this->color_labels_[index];
	delete this->buttons_[index];
	delete this->row_layouts_[index];

	this->line_edits_.remove(index);
	this->color_buttons_.remove(index);
	this->color_labels_.remove(index);
	this->buttons_.remove(index);
	this->row_layouts_.remove(index);
};

void MultiLineEditWithColorChoiceLayout::s_add_choice() {

	QLineEdit* line = new QLineEdit(this);
	line->setFixedSize(120, 30);
	this->line_edits_ << line;

	G_SET_NEW_BUTTON(color_button, "Choose Color", soap::MiddleSize);
	this->color_buttons_ << color_button;

	connect(color_button, &QPushButton::clicked, this, &MultiLineEditWithColorChoiceLayout::s_choose_color);

	QColor color = custom::random_color();
	QLabel* color_label = new QLabel(color.name(), this);
	color_label->setStyleSheet("QLabel{color:" + color.name() + "}");
	this->color_labels_ << color_label;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;
	connect(button, &QPushButton::clicked, this, &MultiLineEditWithColorChoiceLayout::s_delete_choice);

	G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line, color_button, color_label, button);
	this->row_layouts_ << row_layout;
};

QString MultiLineEditWithColorChoiceLayout::current_value() {
	QString ret;
	int size = this->buttons_.size();
	for (int i = 0; i < size; ++i) {
		ret += SOAP_DELIMITER + this->line_edits_[i]->text() + SOAP_DELIMITER + this->color_labels_[i]->text();
	}
	return ret;
};