#include "MultiCheckBoxLayout.h"

#include "Identifier.h"

MultiCheckBoxLayout::MultiCheckBoxLayout(
    const QStringList& item_names,
    QWidget* parent
):
    QWidget(parent),
    item_names_(item_names)
{
    this->main_layout_ = new QVBoxLayout;
    this->setLayout(this->main_layout_);

    this->all_box_ = new QCheckBox("Choose All");
    this->main_layout_->addWidget(this->all_box_);
    this->all_box_->setChecked(false);

    connect(this->all_box_, &QCheckBox::stateChanged, this, &MultiCheckBoxLayout::s_set_all);

    for (auto&& item_name : this->item_names_) {

        QStringList item = item_name.split(':');

        QCheckBox* check_box = new QCheckBox(item[0]);

        if (item.size() > 1 && item[1] == "true") {
            check_box->setChecked(true);
        }
        else {
            check_box->setChecked(false);
        }

        this->checkboxes_ << check_box;

        this->main_layout_->addWidget(check_box);
    }
};

void MultiCheckBoxLayout::s_set_all(int state) {
    
    bool checked = state == Qt::Checked;

    for (auto cb : this->checkboxes_) {
        cb->setChecked(checked);
    }
};

QString MultiCheckBoxLayout::current_value() {
    QString ret;

    for (auto check_box : this->checkboxes_) {

        if (check_box->isChecked()) {
            ret += SOAP_DELIMITER + check_box->text();
        }

    }

    return ret;
};