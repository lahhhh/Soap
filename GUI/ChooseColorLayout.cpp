#include "ChooseColorLayout.h"

#include <QColorDialog>

ChooseColorLayout::ChooseColorLayout(
	const QString& default_color_name,
	QWidget* parent
):
	QWidget(parent)
{
	QVBoxLayout* main_layout = new QVBoxLayout;
	
	this->button_ = new QPushButton(this);

	this->button_->setFixedSize(soap::MiddleSize);

	QColor color = QColor::fromString(default_color_name);
	if (!color.isValid()) {
		color = Qt::black;
	}

	this->color_name_ = color.name();

	this->button_->setText(this->color_name_);
	this->button_->setStyleSheet("QPushButton{color:" + this->color_name_ + "}");

	connect(this->button_, &QPushButton::clicked, this, &ChooseColorLayout::s_choose_color);

	main_layout->addWidget(this->button_);
	this->setLayout(main_layout);
};

QString ChooseColorLayout::current_value() {

	return this->color_name_;
};

void ChooseColorLayout::s_choose_color() {

	QColor color = QColorDialog::getColor(QColor(this->color_name_), this, "Choose Color");

	if (color.isValid()) {

		this->color_name_ = color.name();

		this->button_->setText(this->color_name_);
		this->button_->setStyleSheet("QPushButton{color:" + this->color_name_ + "}");
	}
};