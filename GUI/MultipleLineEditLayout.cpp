#include "MultipleLineEditLayout.h"

#include "CommonDialog.h"
#include "StringVector.h"

#include "SoapGUI.h"

MultipleLineEditLayout::MultipleLineEditLayout(
    const QString & name, 
    SignalEmitter* signal_emitter,
    const QStringList& item_names,
    QWidget *parent
) : 
    QWidget(parent), 
    name_(name),
    signal_emitter_(signal_emitter)    
{
    this->main_layout_ = new QVBoxLayout;
    setLayout(this->main_layout_);

    this->item_layout_ = new QVBoxLayout;

    this->main_layout_->addLayout(this->item_layout_);

    G_SET_BUTTON(this->add_button_, "+ " + this->name_, QSize(100, 30));
    connect(this->add_button_, &QPushButton::clicked, this, &MultipleLineEditLayout::s_add_line);

    if (signal_emitter != nullptr) {
        G_SET_BUTTON(this->import_button_, "import", QSize(50, 30));
        connect(this->import_button_, &QPushButton::clicked, this, &MultipleLineEditLayout::s_import);

        G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_, this->import_button_);
    }
    else {
        G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_);
    }

    if (item_names.isEmpty()) {
        this->add_item("");
    }
    else {
        this->add_items(item_names);
    }

}

void MultipleLineEditLayout::s_delete_line(){
    QPushButton * button = dynamic_cast<QPushButton *>(QObject::sender());
    int index = this->buttons_.indexOf(button);

    this->item_layout_->removeItem(this->row_layouts_[index]);

    delete this->line_edits_[index];
    delete this->buttons_[index];
    delete this->row_layouts_[index];

    this->line_edits_.remove(index);
    this->buttons_.remove(index);
    this->row_layouts_.remove(index);
};

void MultipleLineEditLayout::s_import() {

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

    this->add_items(imports);
};

void MultipleLineEditLayout::clear() {

    for (auto* layout : this->row_layouts_) {
        this->item_layout_->removeItem(layout);
    }

    for (auto* d : this->line_edits_) {
        delete d;
    }
    for (auto* d : this->buttons_) {
        delete d;
    }
    for (auto* d : this->row_layouts_) {
        delete d;
    }

    this->line_edits_.clear();
    this->buttons_.clear();
    this->row_layouts_.clear();
};

void MultipleLineEditLayout::set_items(const QStringList& items) {

    this->clear();

    this->add_items(items);
};

void MultipleLineEditLayout::add_items(const QStringList& items) {

    for (auto&& s : items) {
        this->add_item(s);
    }
};

void MultipleLineEditLayout::add_item(const QString& s) {
    G_SET_NEW_LINEEDIT(line, s, soap::MiddleSize);
    this->line_edits_ << line;

    G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
    this->buttons_ << button;
    connect(button, &QPushButton::clicked, this, &MultipleLineEditLayout::s_delete_line);

    G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line, button);
    this->row_layouts_ << row_layout;
};

void MultipleLineEditLayout::s_add_line(){

    this->add_item("");
};

QString MultipleLineEditLayout::current_value(){
    QString ret;
    for(const auto ptr : this->line_edits_){
        ret += SOAP_DELIMITER + ptr->text();
    }
    return ret;
};

