#pragma once

#include "Identifier.h"

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>

class CompareLayout :
    public QWidget
{
public:
    CompareLayout(
        const QString& name,
        const QMap<QString, QStringList>& factors,
        int type,
        QWidget* parent = nullptr
    );

    QString current_value();

private:

    QString name_;

    QMap<QString, QStringList> factors_;

    int type_{ 0 };

    QVBoxLayout* main_layout_{ nullptr };

    QComboBox* factor_name_box_{ nullptr };
    QComboBox* group_box_1_{ nullptr };
    QComboBox* group_box_2_{ nullptr };

    QStringList get_1_();
    QStringList get_2_();

private slots:

	void s_check_factor();
};

