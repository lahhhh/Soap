#include "MultipleDoubleLineEditWithCompleterLayout.h"
#include "Identifier.h"

MultipleDoubleLineEditWithCompleterLayout::MultipleDoubleLineEditWithCompleterLayout(
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

    this->completer_ = new QCompleter(this->item_names_, this);

    QHBoxLayout* row_layout = new QHBoxLayout;


    QLineEdit* line1 = new QLineEdit(this);
    line1->setCompleter(this->completer_);
    line1->adjustSize();
    line1->setFixedHeight(25);
    this->line_edits1_ << line1;

    QLineEdit* line2 = new QLineEdit(this);
    line2->setFixedSize(200, 25);
    line2->setPlaceholderText(this->hint_);
    this->line_edits2_ << line2;

    QPushButton* button = new QPushButton("-", this);
    button->setFixedSize(30, 25);
    this->buttons_ << button;

    connect(button, &QPushButton::clicked, this, &MultipleDoubleLineEditWithCompleterLayout::s_delete_line);

    row_layout->addWidget(line1);
    row_layout->addWidget(line2);
    row_layout->addWidget(button);


    this->row_layouts_ << row_layout;
    this->item_layout_->addLayout(row_layout);

    this->main_layout_->addLayout(this->item_layout_);

    this->add_button_ = new QPushButton("+" + this->name_, this);
    this->add_button_->setFixedSize(120, 25);
    connect(this->add_button_, &QPushButton::clicked, this, &MultipleDoubleLineEditWithCompleterLayout::s_add_line);

    row_layout = new QHBoxLayout;

    row_layout->addStretch();
    row_layout->addWidget(this->add_button_);
    row_layout->addStretch();

    this->main_layout_->addLayout(row_layout);

}

void MultipleDoubleLineEditWithCompleterLayout::s_delete_line() {
    QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
    int index = this->buttons_.indexOf(button);
    this->item_layout_->removeItem(this->row_layouts_[index]);

    delete this->line_edits1_[index];
    delete this->line_edits2_[index];
    delete this->buttons_[index];
    delete this->row_layouts_[index];

    this->line_edits1_.remove(index);
    this->line_edits2_.remove(index);
    this->buttons_.remove(index);
    this->row_layouts_.remove(index);
};

void MultipleDoubleLineEditWithCompleterLayout::s_add_line() {

    QHBoxLayout* row_layout = new QHBoxLayout;

    QLineEdit* line1 = new QLineEdit(this);
    line1->setCompleter(this->completer_);
    line1->adjustSize();
    line1->setFixedHeight(25);
    this->line_edits1_ << line1;

    QLineEdit* line2 = new QLineEdit(this);
    line2->setFixedSize(200, 25);
    line2->setPlaceholderText(this->hint_);
    this->line_edits2_ << line2;

    QPushButton* button = new QPushButton("-", this);
    button->setFixedSize(30, 25);
    this->buttons_ << button;

    connect(button, &QPushButton::clicked, this, &MultipleDoubleLineEditWithCompleterLayout::s_delete_line);

    row_layout->addWidget(line1);
    row_layout->addWidget(line2);
    row_layout->addWidget(button);


    this->row_layouts_ << row_layout;
    this->item_layout_->addLayout(row_layout);
};

QString MultipleDoubleLineEditWithCompleterLayout::current_value() {
    QString ret;
    int size = this->buttons_.size();
    for (int i = 0; i < size; ++i) {
        ret += SOAP_DELIMITER + this->line_edits1_[i]->text() + SOAP_DELIMITER + this->line_edits2_[i]->text();
    }
    return ret;
};