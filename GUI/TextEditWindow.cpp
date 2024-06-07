#include "TextEditWindow.h"

#include "YesOrNoDialog.h"

#include "SoapGUI.h"

TextEditWindow::TextEditWindow(
	QString* note,
	void* from,
	SignalEmitter* signal_emitter) 
	:
	QMainWindow(signal_emitter->widget_),
	mode_(WorkMode::Note),
	note_(note),
	from_(from),
	signal_emitter_(signal_emitter)
{
	//connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &TextEditWindow::s_check_data);
	//connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &TextEditWindow::s_check_data);	

	this->set_layout();

	this->text_edit_->setHtml(*this->note_);

	this->set_property();

	connect(this->text_edit_, &QTextEdit::textChanged, this, &TextEditWindow::s_edit);
};

TextEditWindow::TextEditWindow(
	QStringList* string_vector_data,
	void* from,
	SignalEmitter* signal_emitter) :
	QMainWindow(signal_emitter->widget_),
	mode_(WorkMode::StringVectorEdit),
	string_vector_data_(string_vector_data),
	from_(from),
	signal_emitter_(signal_emitter)
{
	//connect(this->signal_emitter_, &SignalEmitter::x_data_delete_soon, this, &TextEditWindow::s_check_data);
	//connect(this->signal_emitter_, &SignalEmitter::x_data_edit_soon, this, &TextEditWindow::s_check_data);

	this->set_layout();

	for (auto&& str : *string_vector_data) {
		this->text_edit_->append(str);
	}

	this->set_property();

	connect(this->text_edit_, &QTextEdit::textChanged, this, &TextEditWindow::s_edit);
};

void TextEditWindow::set_layout() {

	this->main_layout_ = new QHBoxLayout;

	this->text_edit_ = new QTextEdit();
	this->main_layout_->addWidget(this->text_edit_);

	this->main_interface_ = new QWidget(this);
	this->main_interface_->setLayout(this->main_layout_);
	this->setCentralWidget(this->main_interface_);
};

void TextEditWindow::set_property() {

	this->setAttribute(Qt::WA_DeleteOnClose);
	G_SET_ICON;
	this->resize(600, 600);
	this->setWindowTitle("Edit");

	this->show();
};

void TextEditWindow::s_edit() {
	this->edited_ = true;
};

void TextEditWindow::s_check_data(void* data, soap::VariableType type, void* item) {
	if (this->from_ == data) {

		this->valid_ = false;

		this->close();
	}
};

void TextEditWindow::view(
	QString* note,
	void* from,
	SignalEmitter* signal_emitter) {
	TextEditWindow* window = new TextEditWindow(note, from, signal_emitter);
};

void TextEditWindow::view(
	QStringList* string_vector_data,
	void* from,
	SignalEmitter* signal_emitter) {
	TextEditWindow* window = new TextEditWindow(string_vector_data, from, signal_emitter);
};

void TextEditWindow::closeEvent(QCloseEvent* e) {

	if (this->edited_ && this->valid_) {

		bool save = YesOrNoDialog::get_response("Close", "Save change?");

		if (save) {
			if (this->mode_ == WorkMode::Note) {
				*this->note_ = this->text_edit_->toHtml();
			}
			else if (this->mode_ == WorkMode::StringVectorEdit) {
				*this->string_vector_data_ = _Cs split_lines(this->text_edit_->toPlainText());
				this->signal_emitter_->x_update_interface();
			}
		}
	}
};