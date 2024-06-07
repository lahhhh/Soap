#pragma once

#include "Identifier.h"

#include <QDialog>

#include "EnrichmentGroupLayout.h"
#include "Switch.h"

class EnrichmentMultiGroupPlotDialog :
    public QDialog
{
public:
    EnrichmentMultiGroupPlotDialog(const QMap<QString, QStringList>& group_map_, QWidget* parent = nullptr);

    static QPair<bool, QList<QStringList>> get_plot_settings(const QMap<QString, QStringList>& group_map, QWidget* parent = nullptr);

    bool is_accepted_ = false;

    QMap<QString, QStringList> group_map_;

    QVBoxLayout* all_layout_;
    QVBoxLayout* main_layout_;

    QVBoxLayout* group_layout_;

    QList<QHBoxLayout*> row_layouts_;

    QList<EnrichmentGroupLayout*> layout_groups_;

    QList<QPushButton*> delete_buttons_;

    QPushButton* add_button_;

    QLabel* capitalize_label_;
    Switch* capitalize_switch_;

    QPushButton* finish_button_;
    QPushButton* cancel_button_;

private slots:

    void s_add_group();

    void s_delete_group();

    void accept();
    void reject();
};

