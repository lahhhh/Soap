#include "SimpleChoiceWithLineEditLayout.h"

#include "SoapGUI.h"

SimpleChoiceWithLineEditLayout::SimpleChoiceWithLineEditLayout(
    const QString& name, 
    const QString& hint, 
    const QStringList& item_names,
    QWidget* parent
) : 
    QWidget(parent), 
    name_(name), 
    hint_(hint), 
    item_names_(item_names)
{
    this->main_layout_ = new QVBoxLayout;
    setLayout(this->main_layout_);

    this->item_layout_ = new QVBoxLayout;

    G_SET_NEW_COMBOBOX(box, this->item_names_, 30)
    this->combo_boxes_ << box;

    QLineEdit* line = new QLineEdit(this);
    line->setFixedSize(200, 30);
    line->setPlaceholderText(this->hint_);
    line->setStyleSheet("::placeholder{color:grey; font-style:italic}");
    this->line_edits_ << line;

    G_SET_NEW_BUTTON(button, "-", QSize(30, 30))
    this->buttons_ << button;

    connect(button, &QPushButton::clicked, this, &SimpleChoiceWithLineEditLayout::s_delete_choice);

    G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, line, button)
    this->row_layouts_ << row_layout;

    this->main_layout_->addLayout(this->item_layout_);

    G_SET_BUTTON(this->add_button_, "+ " + this->name_, soap::MiddleSize)
    connect(this->add_button_, &QPushButton::clicked, this, &SimpleChoiceWithLineEditLayout::s_add_choice);

    G_ADD_SINGLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_)
}

void SimpleChoiceWithLineEditLayout::s_delete_choice() {
    QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
    int index = this->buttons_.indexOf(button);
    this->item_layout_->removeItem(this->row_layouts_[index]);

    delete this->combo_boxes_[index];
    delete this->line_edits_[index];
    delete this->buttons_[index];
    delete this->row_layouts_[index];

    this->combo_boxes_.remove(index);
    this->line_edits_.remove(index);
    this->buttons_.remove(index);
    this->row_layouts_.remove(index);
};

void SimpleChoiceWithLineEditLayout::s_add_choice() {
    G_SET_NEW_COMBOBOX(box, this->item_names_, 30)
    this->combo_boxes_ << box;

    QLineEdit* line = new QLineEdit(this);
    line->setFixedSize(200, 30);
    line->setPlaceholderText(this->hint_);
    line->setStyleSheet("::placeholder{color:grey; font-style:italic}");
    this->line_edits_ << line;

    G_SET_NEW_BUTTON(button, "-", QSize(30, 30))
    this->buttons_ << button;
    connect(button, &QPushButton::clicked, this, &SimpleChoiceWithLineEditLayout::s_delete_choice);

    G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, box, line, button)
    this->row_layouts_ << row_layout;
};

QString SimpleChoiceWithLineEditLayout::current_value() {
    QString ret;
    int size = this->buttons_.size();
    for (int i = 0; i < size; ++i) {
        ret += SOAP_DELIMITER + this->combo_boxes_[i]->currentText() + SOAP_DELIMITER + this->line_edits_[i]->text();
    }
    return ret;
};