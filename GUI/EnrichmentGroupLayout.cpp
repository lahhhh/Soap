#include "EnrichmentGroupLayout.h"

#include <QColorDialog>

#include "Custom.h"

EnrichmentGroupLayout::EnrichmentGroupLayout(
	const QString& name, 
	const QMap<QString, QStringList>& map, 
	QWidget* parent
): 
	QWidget(parent) 
{
	this->main_layout_ = new QVBoxLayout;

	this->select_layout_ = new FactorDoubleLineEditWithCompleterLayout(name, "Custom Label", map, this);
	this->main_layout_->addWidget(this->select_layout_);

	QHBoxLayout* row_layout = new QHBoxLayout;

	this->legend_color_button_ = new QPushButton("Choose Color", this);
	QColor color = custom::random_color();
	this->legend_color_label_ = new QLabel(color.name(), this);
	this->legend_color_label_->setStyleSheet("QLabel{color:" + color.name() + "}");

	connect(this->legend_color_button_, &QPushButton::clicked, this, &EnrichmentGroupLayout::s_choose_color);

	row_layout->addWidget(this->legend_color_button_);
	row_layout->addWidget(this->legend_color_label_);
	this->main_layout_->addLayout(row_layout);

	setLayout(this->main_layout_);
};

void EnrichmentGroupLayout::s_choose_color() {
	QColor color = QColorDialog::getColor(QColor(this->legend_color_label_->text()), this, "Choose Color for Legend");
	if (color.isValid()) {
		this->legend_color_label_->setText(color.name());
		this->legend_color_label_->setStyleSheet("QLabel{color:" + color.name() + "}");
	}
};

QStringList EnrichmentGroupLayout::current_value() {
	return QStringList() << this->select_layout_->current_value() << this->legend_color_label_->text();
};