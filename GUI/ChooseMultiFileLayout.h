#pragma once

#include "Identifier.h"

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>

class ChooseMultiFileLayout :
    public QWidget
{
public:
    ChooseMultiFileLayout(
        const QString& file_type, 
        const QStringList& file_names, 
        QWidget* parent);

    QString current_value();

    void add_file(const QString& file_name);

private:

    QString file_type_;

    QVBoxLayout* main_layout_;

    QVBoxLayout* item_layout_;

    QPushButton* and_button_;

    QList<QHBoxLayout*> row_layouts_;

    QList<QLabel*> labels_;

    QList<QPushButton*> buttons_;    

private slots:

    void s_add_file();

    void s_delete_file();

};

