#include "ChooseMultiFileLayout.h"

#include <QFileDialog>
#include "MessageDialog.h"

#include "SoapGUI.h"

ChooseMultiFileLayout::ChooseMultiFileLayout(
    const QString& file_type, 
    const QStringList& file_names, 
    QWidget* parent
) :
	QWidget(parent),
	file_type_(file_type)
{
    this->main_layout_ = new QVBoxLayout;
    this->setLayout(this->main_layout_);

    this->item_layout_ = new QVBoxLayout;

    this->main_layout_->addLayout(this->item_layout_);

    G_SET_BUTTON(this->and_button_, "Add File", soap::MiddleSize);
    connect(this->and_button_, &QPushButton::clicked, this, &ChooseMultiFileLayout::s_add_file);

    G_NEW_SINGLE_ITEM_ROWLAYOUT(row_layout, this->and_button_);

    this->main_layout_->addLayout(row_layout);

    if (!file_names.isEmpty()) {
        for (const auto& file_name : file_names) {
            this->add_file(file_name);
        }
    }
};

void ChooseMultiFileLayout::add_file(const QString& file_name) {

    QFileInfo file_info(file_name);
    if (!file_info.exists() || !file_info.isFile() || !file_info.isReadable()) {
        return;
    }

    G_SET_NEW_LABEL_ADJUST_SIZE(label, file_name);
    this->labels_ << label;

    G_SET_NEW_BUTTON(button, "-", QSize(30, 25));
    this->buttons_ << button;
    connect(button, &QPushButton::clicked, this, &ChooseMultiFileLayout::s_delete_file);

    G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->item_layout_, row_layout, label, button);
    this->row_layouts_ << row_layout;
};

void ChooseMultiFileLayout::s_add_file() {

    auto file_names = QFileDialog::getOpenFileNames(this, "Choose File", "", this->file_type_);
    if (file_names.isEmpty()) {
        return;
    }

    for (auto&& f : file_names) {

        this->add_file(f);
    }
};

void ChooseMultiFileLayout::s_delete_file() {
    QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
    int index = this->buttons_.indexOf(button);
    this->item_layout_->removeItem(this->row_layouts_[index]);

    delete this->labels_[index];
    delete this->buttons_[index];
    delete this->row_layouts_[index];

    this->labels_.remove(index);
    this->buttons_.remove(index);
    this->row_layouts_.remove(index);
};

QString ChooseMultiFileLayout::current_value() {

    QString ret;
    
    for (const auto ptr : this->labels_) {
        ret += SOAP_DELIMITER + ptr->text();
    }
    
    return ret;
};

