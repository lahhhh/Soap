#include "FactorDoubleLineEditWithCompleterLayout.h"

#include "Identifier.h"

#include "SoapGUI.h"

FactorDoubleLineEditWithCompleterLayout::FactorDoubleLineEditWithCompleterLayout(const QString& name, const QString& hint, const QMap<QString, QStringList>& factor_map, QWidget* parent) : 
    QWidget(parent), 
    name_(name), 
    hint_(hint), 
    factors_(factor_map.keys()), 
    factor_map_(factor_map)
{
    for (const auto& factor : this->factor_map_.keys()) {
        this->completers_[factor] = new QCompleter(this->factor_map_[factor], this);
    }

    this->main_layout_ = new QVBoxLayout;
    setLayout(this->main_layout_);

    G_SET_COMBOBOX(this->factor_box_, this->factor_map_.keys(), 30)
    connect(this->factor_box_, &QComboBox::currentIndexChanged, this, &FactorDoubleLineEditWithCompleterLayout::s_factor_changed);

    this->factor_line_edit_ = new QLineEdit(this);
    this->factor_line_edit_->setPlaceholderText(this->hint_);
    this->factor_line_edit_->setFixedSize(200, 30);

    G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, row_layout, this->factor_box_, this->factor_line_edit_)

    this->item_layout_ = new QVBoxLayout;
    this->main_layout_->addLayout(this->item_layout_);

    G_SET_BUTTON(this->add_button_, "+ " + this->name_, soap::MiddleSize)
    connect(this->add_button_, &QPushButton::clicked, this, &FactorDoubleLineEditWithCompleterLayout::s_add_choice);

    G_ADD_SINGLE_ITEM_ROWLAYOUT(this->main_layout_, row_layout, this->add_button_)
}

void FactorDoubleLineEditWithCompleterLayout::s_factor_changed() {
    int size = this->buttons_.size();
    if (size < 1)return;

    for (int i = 0; i < size; ++i) {
        this->item_layout_->removeItem(this->row_layouts_[i]);
        delete this->line_edits_1_[i];
        delete this->line_edits_2_[i];
        delete this->buttons_[i];
        delete this->row_layouts_[i];
    }
    
    this->line_edits_1_.clear();
    this->line_edits_2_.clear();
    this->buttons_.clear();
    this->row_layouts_.clear();
};

void FactorDoubleLineEditWithCompleterLayout::s_delete_choice() {
    QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
    int index = this->buttons_.indexOf(button);
    
    this->item_layout_->removeItem(this->row_layouts_[index]);
    
    delete this->line_edits_1_[index];
    delete this->line_edits_2_[index];
    delete this->buttons_[index];
    delete this->row_layouts_[index];
    
    this->line_edits_1_.remove(index);
    this->line_edits_2_.remove(index);
    this->buttons_.remove(index);
    this->row_layouts_.remove(index);
};

void FactorDoubleLineEditWithCompleterLayout::s_add_choice() {

    QLineEdit* line1 = new QLineEdit(this);
    line1->setCompleter(this->completers_[this->factor_box_->currentText()]);
    line1->adjustSize();
    line1->setFixedHeight(30);
    this->line_edits_1_ << line1;

    QLineEdit* line2 = new QLineEdit(this);
    line2->setFixedSize(200, 30);
    line2->setPlaceholderText(this->hint_);
    this->line_edits_2_ << line2;

    G_SET_NEW_BUTTON(button, "-", QSize(30, 30))
    this->buttons_ << button;
    connect(button, &QPushButton::clicked, this, &FactorDoubleLineEditWithCompleterLayout::s_delete_choice);

    G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, line1, line2, button)
    this->row_layouts_ << row_layout;
};

QString FactorDoubleLineEditWithCompleterLayout::current_value() {
    QString ret = this->factor_box_->currentText() + SOAP_DELIMITER;
    ret += this->factor_line_edit_->text();

    int size = this->buttons_.size();

    for (int i = 0; i < size; ++i) {
        ret += SOAP_DELIMITER + this->line_edits_1_[i]->text();
        ret += SOAP_DELIMITER + this->line_edits_2_[i]->text();
    }

    return ret;
};