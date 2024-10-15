#include "LogicGroup.h"
#include "Custom.h"
#include "Identifier.h"

#include "CustomPlot.h"

#include "SoapGUI.h"

LogicGroup::LogicGroup(LogicHandler* lh, QWidget* parent) :
	QWidget(parent),
	logic_handler_(lh)
{
	this->full_layout_ = new QHBoxLayout;
	this->full_layout_->setContentsMargins(QMargins(0, 0, 0, 0));

	G_SET_LABEL(this->feature_label_, "Feature", QSize(60, 30));

	G_SET_LINEEDIT(this->feature_name_, "", soap::MiddleSize);

	this->feature_completer_ = new QCompleter(lh->data_names(), this);
	this->feature_name_->setCompleter(this->feature_completer_);

	G_SET_COMBOBOX_FIXED_WIDTH(this->compare_box_, QStringList() << QString(u'=') << QString(u'<') << QString(u'≤') << QString(u'>') << QString(u'≥') << QString(u'≠'), 50, 30);
	
	G_SET_LINEEDIT(this->value_line_edit_, "", soap::MiddleSize);
	this->value_line_edit_->setValidator(new QDoubleValidator(this));

	G_SET_COMBOBOX_FIXED_WIDTH(this->value_box_, QStringList(), 0, 30);
	this->value_box_->setVisible(false);

	G_SET_BUTTON(this->delete_button_, "-", QSize(20, 30));
	connect(this->delete_button_, &QPushButton::clicked, this, &LogicGroup::__s_delete_this);

	this->full_layout_->addWidget(this->feature_label_);
	this->full_layout_->addWidget(this->feature_name_);
	this->full_layout_->addWidget(this->compare_box_);
	this->full_layout_->addWidget(this->value_line_edit_);
	this->full_layout_->addWidget(this->value_box_);
	this->full_layout_->addWidget(this->delete_button_);

	setLayout(this->full_layout_);

	connect(this->feature_name_, &QLineEdit::editingFinished, this, &LogicGroup::s_check_feature);
}

LogicGroup::LogicGroup(const CustomMatrix* metadata, const QStringList* gene_names, QWidget* parent) :
	QWidget(parent),
	metadata_(metadata),
	gene_names_(gene_names)
{
	this->full_layout_ = new QHBoxLayout;
	this->full_layout_->setContentsMargins(QMargins(0, 0, 0, 0));

	G_SET_LABEL(this->feature_label_, "Feature", QSize(60, 30));

	G_SET_LINEEDIT(this->feature_name_, "", soap::MiddleSize);

	this->feature_completer_ = new QCompleter(QStringList() << *this->gene_names_ << this->metadata_->colnames_, this);
	this->feature_name_->setCompleter(this->feature_completer_);

	G_SET_COMBOBOX_FIXED_WIDTH(this->compare_box_, QStringList() << QString(u'=') << QString(u'<') << QString(u'≤') << QString(u'>') << QString(u'≥') << QString(u'≠'), 50, 30);

	G_SET_LINEEDIT(this->value_line_edit_, "", soap::MiddleSize);
	this->value_line_edit_->setValidator(new QDoubleValidator(this));

	G_SET_COMBOBOX_FIXED_WIDTH(this->value_box_, QStringList(), 0, 30);
	this->value_box_->setVisible(false);

	G_SET_BUTTON(this->delete_button_, "-", QSize(20, 30));
	connect(this->delete_button_, &QPushButton::clicked, this, &LogicGroup::__s_delete_this);

	this->full_layout_->addWidget(this->feature_label_);
	this->full_layout_->addWidget(this->feature_name_);
	this->full_layout_->addWidget(this->compare_box_);
	this->full_layout_->addWidget(this->value_line_edit_);
	this->full_layout_->addWidget(this->value_box_);
	this->full_layout_->addWidget(this->delete_button_);

	setLayout(this->full_layout_);

	connect(this->feature_name_, &QLineEdit::editingFinished, this, &LogicGroup::s_check_feature);
}

void LogicGroup::__s_delete_this() {
	emit x_delete_this(this);
};

void LogicGroup::to_factor_mode(const QStringList& levels) {

	this->compare_box_->clear();
	this->compare_box_->addItems(QStringList() << QString(u'=') << QString(u'≠'));
	this->value_line_edit_->clear();
	this->value_line_edit_->setVisible(false);
	this->value_box_->clear();
	auto width = custom_plot::utility::get_max_text_width(levels, QFont("Arial", 15));
	this->value_box_->setFixedWidth(width + 20);
	this->value_box_->setVisible(true);
	this->value_box_->addItems(levels);

	this->factor_mode_ = true;
};

void LogicGroup::to_numeric_mode() {

	this->compare_box_->clear();
	this->compare_box_->addItems(QStringList() << QString(u'=') << QString(u'<') << QString(u'≤') << QString(u'>') << QString(u'≥') << QString(u'≠') );
	this->value_line_edit_->clear();
	this->value_line_edit_->setVisible(true);
	this->value_box_->clear();
	this->value_box_->setVisible(false);

	this->factor_mode_ = false;
};

void LogicGroup::s_check_feature() {

	QString name = this->feature_name_->text();
	if (name.isEmpty()) {
		return;
	}

	auto [valid, type, levels] = this->logic_handler_->get_content(name);

	if (!valid) {
		this->to_numeric_mode(); // maybe peak name in single cell multiome
		this->feature_name_->setStyleSheet("QLineEdit{border:2px solid red}");
		return;
	}

	if (levels.isEmpty()) {
		this->to_numeric_mode();
	}
	else {
		this->to_factor_mode(levels);
	}

	this->feature_name_->setStyleSheet("QLineEdit{border:2px solid green}");
};

SingleLogic LogicGroup::get_filter() {

	if (this->factor_mode_) {
		return std::make_tuple(this->feature_name_->text(), this->compare_box_->currentText(), this->value_box_->currentText());
	}
	else {
		return std::make_tuple(this->feature_name_->text(), this->compare_box_->currentText(), this->value_line_edit_->text());
	}
};

QString LogicGroup::current_value() {

	if (this->factor_mode_) {
		return this->feature_name_->text() + SOAP_DELIMITER + this->compare_box_->currentText() + SOAP_DELIMITER + this->value_box_->currentText();
	}
	else {
		return this->feature_name_->text() + SOAP_DELIMITER + this->compare_box_->currentText() + SOAP_DELIMITER + this->value_line_edit_->text();
	}
};
