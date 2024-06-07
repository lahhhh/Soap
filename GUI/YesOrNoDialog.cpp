#include "YesOrNoDialog.h"
#include "Identifier.h"

#include "SoapGUI.h"

YesOrNoDialog::YesOrNoDialog(const QString &title, const QString &message) : QDialog()
{
    this->main_layout_ = new QVBoxLayout;

    this->message_label_ = new QLabel(message, this);
    this->message_label_->adjustSize();

    this->main_layout_->addWidget(this->message_label_);

    G_SET_BUTTON_FINISH_STYLE(this->yes_button_,"Yes")

    G_SET_BUTTON_CANCEL_STYLE(this->no_button_,"No")
    
    G_DOUBLE_ITEM_ROWLAYOUT(this->bottom_layout_,this->yes_button_,this->no_button_)

    this->main_layout_->addLayout(this->bottom_layout_);

    setLayout(this->main_layout_);

    connect(this->yes_button_, &QPushButton::clicked, this, &YesOrNoDialog::accept);
    connect(this->no_button_, &QPushButton::clicked, this, &YesOrNoDialog::reject);
    G_SET_ICON;
    this->setWindowTitle(title);
    this->exec();
}

bool YesOrNoDialog::get_response(const QString & title, const QString & message){
    YesOrNoDialog dlg(title, message);
    if(!dlg.is_accepted_)
        return false;
    else
        return true;
};

void YesOrNoDialog::accept(){
    QDialog::accept();
    this->is_accepted_ = true;
}

void YesOrNoDialog::reject(){
    QDialog::reject();
}
