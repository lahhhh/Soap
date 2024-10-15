#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>

#include "Switch.h"

struct COMPARE_PLOT_ELEMENT {
    QString feature;
    bool normalize{ false };
    bool use_gene_activity{ false };
    bool show_p_value{ false };
    QString method;
    QString group;
    QString type;
    QList<QPair<QString, QString>> comparisons;
};

class ComparePlotDialog
    : public QDialog
{
    Q_OBJECT
public:

    ComparePlotDialog(
        const QMap<QString, QStringList>& factors,
        const QStringList& valid_features
    );

    QMap<QString, QStringList> factors_;

    QVBoxLayout* main_layout_{ nullptr };

    QLabel* feature_label_{ nullptr };
    QLineEdit* choose_feature_{ nullptr };

    QLabel* normalize_label_{ nullptr };
    Switch* normalize_{ nullptr };

    QLabel* gene_activity_label_{ nullptr };
    Switch* gene_activity_{ nullptr };

    QLabel* p_value_label_{ nullptr };
    Switch* p_value_{ nullptr };

    QLabel* method_label_{ nullptr };
    QComboBox* choose_method_{ nullptr };

    QLabel* group_label_{ nullptr };
    QComboBox* choose_group_{ nullptr };

    QLabel* type_label_{ nullptr };
    QComboBox* draw_type_{ nullptr };

    QWidget* compare_widget_{ nullptr };
    QLabel* comparison_label_{ nullptr };
    QVBoxLayout* choose_comparison_{ nullptr };
    QHBoxLayout* comparison_layout_{ nullptr };
    QVBoxLayout* comparison_inner_layout_{ nullptr };
    QHBoxLayout* comparison_outer_layout_{ nullptr };
    QList<QComboBox*> comp1_;
    QList<QComboBox*> comp2_;
    QList<QPushButton*> delete_buttons_;
    QList<QHBoxLayout*> row_layouts_;
    QPushButton* add_button_{ nullptr };

    QHBoxLayout* bottom_layout_{ nullptr };
    QPushButton* finish_button_{ nullptr };
    QPushButton* cancel_button_{ nullptr };

    bool is_accepted_{ false };
    bool is_expanding_{ false };

    COMPARE_PLOT_ELEMENT current_value();

    static COMPARE_PLOT_ELEMENT
        get_response(
        const QMap<QString, QStringList>& factors,
        const QStringList& valid_features
    );

private slots:

    void accept();

    void reject();

    void s_change_interface();

    void s_factor_changed();

    void s_add_comparison();
    void s_delete_comparison();
 };

