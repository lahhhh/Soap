#include "ChooseOpenFileLayout.h"

#include <QFileDialog>
#include "MessageDialog.h"

#include "SoapGUI.h"

ChooseOpenFileLayout::ChooseOpenFileLayout(
    const QString& file_type,
    QWidget* parent
) :
    QWidget(parent),
    file_type_(file_type) {
    
    this->main_layout_ = new QHBoxLayout;
    this->setLayout(this->main_layout_);

    G_SET_LABEL(this->file_name_label_, "", soap::MiddleSize);

    G_SET_BUTTON(this->choose_file_button_, "Choose File", soap::MiddleSize);
    connect(this->choose_file_button_, &QPushButton::clicked, this, &ChooseOpenFileLayout::s_choose_file);

    this->main_layout_->addWidget(this->file_name_label_);
    this->main_layout_->addWidget(this->choose_file_button_);
};

void ChooseOpenFileLayout::s_choose_file() {

    QString file_path = QFileDialog::getOpenFileName(this, "Choose File", "", this->file_type_);
    if (file_path.isEmpty())return;

    QFileInfo fi(file_path);
    if (!fi.exists() || !fi.isFile() || !fi.isReadable()) {
        MessageDialog::get_response("File Not Available", file_path + " is not available.");
        return;
    }

    this->file_name_label_->setText(file_path);
};