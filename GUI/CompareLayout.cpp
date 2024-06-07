#include "CompareLayout.h"

#include "SoapGUI.h"

CompareLayout::CompareLayout(
    const QString& name,
    const QMap<QString, QStringList>& factors,
    int type,
    QWidget* parent
):
    QWidget(parent),
    name_(name),
    factors_(factors),
    type_(type)
{
    this->main_layout_ = new QVBoxLayout;
    this->setLayout(this->main_layout_);

    G_SET_COMBOBOX(this->factor_name_box_, factors.keys(), 30);
    G_SET_COMBOBOX(this->group_box_1_, this->get_1_(), 30);
    G_SET_COMBOBOX(this->group_box_2_, this->get_2_(), 30);

    connect(this->factor_name_box_, &QComboBox::currentIndexChanged, this, &CompareLayout::s_check_factor);

    this->main_layout_->addWidget(this->factor_name_box_);
    this->main_layout_->addWidget(this->group_box_1_);
    this->main_layout_->addWidget(this->group_box_2_);
};

QString CompareLayout::current_value() {
    return this->factor_name_box_->currentText()
        + SOAP_DELIMITER
        + this->group_box_1_->currentText()
        + SOAP_DELIMITER
        + this->group_box_2_->currentText();
};

void CompareLayout::s_check_factor() {
    this->group_box_1_->clear();
    this->group_box_1_->addItems(this->get_1_());
    this->group_box_2_->clear();
    this->group_box_2_->addItems(this->get_2_());
}

QStringList CompareLayout::get_1_() {
    QString factor_name = this->factor_name_box_->currentText();

    if (this->type_ == 0) {
        return this->factors_[factor_name];
    }
    else if (this->type_ == 1) {
        return QStringList() << "ALL" << this->factors_[factor_name];
    }
    else if (this->type_ == 2) {
        return this->factors_[factor_name];
    }
    else {
        return {};
    }
};

QStringList CompareLayout::get_2_() {
    QString factor_name = this->factor_name_box_->currentText();

    if (this->type_ == 0) {
        return this->factors_[factor_name];
    }
    else if (this->type_ == 1) {
        return QStringList() << "REST" << this->factors_[factor_name];
    }
    else if (this->type_ == 2) {
        return QStringList() << "REST" << this->factors_[factor_name];
    }
    else {
        return {};
    }
};