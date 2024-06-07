#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>

#include "logicgroup.h"
#include "custommatrix.h"

#include "LogicHandler.h"

using OrLogic = QList<SingleLogic>;

class SubLogicLayout : public QWidget
{
    Q_OBJECT
public:

    SubLogicLayout(const CustomMatrix * metadata, const QStringList * gene_names, QWidget * parent);
    SubLogicLayout(LogicHandler* lh, QWidget* parent);

    const CustomMatrix * metadata_;
    const QStringList * gene_names_;

    LogicHandler* logic_handler_{ nullptr };

    QVBoxLayout * full_layout_;
    QVBoxLayout * main_layout_;
    QPushButton * or_button_;

    QList<QLabel *> or_labels_;
    QList<LogicGroup *> logic_groups_;

    OrLogic get_filters();

    QString current_value();

signals:
    void x_delete_this(SubLogicLayout *);

private slots:
    void s_or_clicked();

    void s_delete_logic_group(LogicGroup * lg);
};

