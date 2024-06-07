#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>

class YesOrNoDialog : public QDialog
{
    Q_OBJECT
public:
    YesOrNoDialog(const QString & title, const QString & message);

    QVBoxLayout* main_layout_;
    QHBoxLayout* bottom_layout_;

    QLabel* message_label_;

    QPushButton* yes_button_;
    QPushButton* no_button_;

    bool is_accepted_ = false;

    bool static get_response(const QString & title, const QString & message);

private slots:
    void accept();
    void reject();
};
