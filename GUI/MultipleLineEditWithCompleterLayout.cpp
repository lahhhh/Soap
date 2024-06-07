#include "MultipleLineEditWithCompleterLayout.h"
#include <QCompleter>

#include "CommonDialog.h"
#include "StringVector.h"

#include "SoapGUI.h"

MultipleLineEditWithCompleterLayout::MultipleLineEditWithCompleterLayout(
	const QString& name,
	const QStringList& completer,
	SignalEmitter* signal_emitter,
	QWidget* parent
) :
	QWidget(parent),
	name_(name),
	completer_(completer),
	signal_emitter_(signal_emitter)
{
	this->main_layout_ = new QVBoxLayout;
	setLayout(this->main_layout_);

	this->item_layout_ = new QVBoxLayout;

	G_SET_NEW_LINEEDIT_WITH_COMPLETER(line, "", completer, soap::MiddleHeight);
	this->line_edits_ << line;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;

	connect(button, &QPushButton::clicked, this, &MultipleLineEditWithCompleterLayout::s_delete_line);

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line, button);

	this->row_layouts_ << row_layout;

	this->main_layout_->addLayout(this->item_layout_);

	G_SET_BUTTON(this->add_button_, "+ " + this->name_, QSize(100, 30));
	connect(this->add_button_, &QPushButton::clicked, this, &MultipleLineEditWithCompleterLayout::s_add_line);

	G_SET_BUTTON(this->import_button_, "import", QSize(50, 30));
	connect(this->import_button_, &QPushButton::clicked, this, &MultipleLineEditWithCompleterLayout::s_import);

	G_ADD_DOUBLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_, this->import_button_);

}

void MultipleLineEditWithCompleterLayout::s_import() {
	if (this->signal_emitter_ == nullptr) {
		return;
	}

	auto info = this->signal_emitter_->get_type_variable(soap::VariableType::StringVector);
	if (info.isEmpty()) {
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Import",
		{ "Data" },
		{soap::InputStyle::ComboBox},
		{info.keys()}
	);

	if (settings.isEmpty()) {
		return;
	}

	StringVector* sv = static_cast<StringVector*>(info[settings[0]]);
	QStringList& imports = sv->data_;

	for (auto&& str : imports) {
		G_SET_NEW_LINEEDIT_WITH_COMPLETER(line, str, this->completer_, soap::MiddleHeight);
		this->line_edits_ << line;

		G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
		this->buttons_ << button;
		connect(button, &QPushButton::clicked, this, &MultipleLineEditWithCompleterLayout::s_delete_line);

		G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line, button);
		this->row_layouts_ << row_layout;
	}
};

void MultipleLineEditWithCompleterLayout::s_delete_line() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
	int index = this->buttons_.indexOf(button);

	this->item_layout_->removeItem(this->row_layouts_[index]);

	delete this->line_edits_[index];
	delete this->buttons_[index];
	delete this->row_layouts_[index];

	this->line_edits_.remove(index);
	this->buttons_.remove(index);
	this->row_layouts_.remove(index);
};

void MultipleLineEditWithCompleterLayout::s_add_line() {

	G_SET_NEW_LINEEDIT_WITH_COMPLETER(line, "", this->completer_, soap::MiddleHeight);
	this->line_edits_ << line;

	G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
	this->buttons_ << button;
	connect(button, &QPushButton::clicked, this, &MultipleLineEditWithCompleterLayout::s_delete_line);

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line, button);
	this->row_layouts_ << row_layout;
};

QString MultipleLineEditWithCompleterLayout::current_value() {
	QString ret;
	for (const auto ptr : this->line_edits_) {
		ret += SOAP_DELIMITER + ptr->text();
	}
	return ret;
};

