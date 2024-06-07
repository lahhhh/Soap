#pragma once

#include "Identifier.h"

#include <QWidget>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLabel>

class ChooseSaveFileLayout :
    public QWidget
{
public:
    ChooseSaveFileLayout(
        const QString& file_type,
        QWidget* parent
    );

    QString current_value() {

        return this->file_name_label_->text();
    };

private:

    QString file_type_;

    QHBoxLayout* main_layout_{ nullptr };

    QLabel* file_name_label_{ nullptr };

    QPushButton* choose_file_button_{ nullptr };

private slots:

    void s_choose_file();

};

