#include "MessageDialog.h"
#include "Identifier.h"

#include "Custom.h"

#include "SoapGUI.h"

static const QFontMetrics fm = QFontMetrics(QFont("Arial"));
static constexpr int LabelWidth = 150;

MessageDialog::MessageDialog(const QString& title, const QString& message)
{
    QVBoxLayout* main_layout = new QVBoxLayout;
    QString processed_message = custom::string_next_line(fm, message, LabelWidth);
    QLabel* message_label = new QLabel(processed_message, this);
    message_label->adjustSize();

    main_layout->addWidget(message_label);

    QPushButton* finish_button;
    QHBoxLayout* bottom_layout;

    G_SET_BUTTON_FINISH_STYLE(finish_button, "OK");
    
    G_ADD_SINGLE_ITEM_ROWLAYOUT(main_layout, bottom_layout, finish_button);

    setLayout(main_layout);

    connect(finish_button, &QPushButton::clicked, this, &QDialog::accept);

    G_SET_ICON;
    this->setWindowTitle(title);
    this->exec();
}

void MessageDialog::get_response(const QString& title, const QString& message) {

    MessageDialog dlg(title, message);
};
