#include "ComparePlotDialog.h"

#include "SoapGUI.h"

#include <QCompleter>

ComparePlotDialog::ComparePlotDialog(
    const QMap<QString, QStringList>& factors,
    const QStringList& valid_features
):
    factors_(factors)
{
    this->main_layout_ = new QVBoxLayout;

    G_SET_LABEL(this->feature_label_, "Feature", soap::MiddleSize);
    G_SET_LINEEDIT_WITH_COMPLETER(this->choose_feature_, "", valid_features, soap::MiddleSize);
    G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, layout, this->feature_label_, this->choose_feature_);

    G_SET_LABEL(this->normalize_label_, "Normalize", soap::MiddleSize);
    G_SET_SWITCH(this->normalize_, true, this->normalize_label_, soap::MiddleSize);
    G_ADD_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, layout, this->normalize_label_, this->normalize_);

    G_SET_LABEL(this->gene_activity_label_, "Use Gene Activity", soap::MiddleSize);
    G_SET_SWITCH(this->gene_activity_, false, this->gene_activity_label_, soap::MiddleSize);
    G_ADD_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, layout, this->gene_activity_label_, this->gene_activity_);

    G_SET_LABEL(this->group_label_, "Factor", soap::MiddleSize);
    G_SET_COMBOBOX(this->choose_group_, factors.keys(), 30);
    G_ADD_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, layout, this->group_label_, this->choose_group_);
    connect(this->choose_group_, &QComboBox::currentTextChanged, this, &ComparePlotDialog::s_factor_changed);

    G_SET_LABEL(this->type_label_, "Show Significance", soap::MiddleSize);
    G_SET_COMBOBOX(this->draw_type_, QStringList() << "All" << "Significant Only" << "Custom", 30);
    G_ADD_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, layout, this->type_label_, this->draw_type_);
    connect(this->draw_type_, &QComboBox::currentTextChanged, this, &ComparePlotDialog::s_change_interface);

    G_SET_LABEL(this->p_value_label_, "Show P Value", soap::MiddleSize);
    G_SET_SWITCH(this->p_value_, false, this->p_value_label_, soap::MiddleSize);
    G_ADD_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(this->main_layout_, layout, this->p_value_label_, this->p_value_);

    G_SET_LABEL(this->comparison_label_, "Comparison", soap::MiddleSize);

    this->choose_comparison_ = new QVBoxLayout;
    this->add_button_ = new QPushButton("+", this);
    this->add_button_->setFixedSize(50, 30);
    connect(this->add_button_, &QPushButton::clicked, this, &ComparePlotDialog::s_add_comparison);
    this->comparison_layout_ = new QHBoxLayout;
    this->comparison_outer_layout_ = new QHBoxLayout;
    this->comparison_inner_layout_ = new QVBoxLayout;
    this->comparison_layout_->addWidget(this->comparison_label_);
    this->choose_comparison_->addLayout(this->comparison_inner_layout_);
    this->choose_comparison_->addWidget(this->add_button_);
    this->comparison_layout_->addLayout(this->choose_comparison_);
    this->compare_widget_ = new QWidget(this);
    this->compare_widget_->setLayout(this->comparison_layout_);
    this->compare_widget_->setVisible(false);
    this->main_layout_->addLayout(this->comparison_outer_layout_);

    G_SET_BUTTON_FINISH_STYLE(this->finish_button_, "Yes");
    G_SET_BUTTON_CANCEL_STYLE(this->cancel_button_, "No");
    G_DOUBLE_ITEM_ROWLAYOUT(this->bottom_layout_, this->finish_button_, this->cancel_button_);
    this->main_layout_->addLayout(this->bottom_layout_);

    this->setLayout(this->main_layout_);

    connect(this->finish_button_, &QPushButton::clicked, this, &ComparePlotDialog::accept);
    connect(this->cancel_button_, &QPushButton::clicked, this, &ComparePlotDialog::reject);
    G_SET_ICON;
    this->setWindowTitle("Compare Plot Settings");
    this->exec();
};

COMPARE_PLOT_ELEMENT
ComparePlotDialog::get_response(
    const QMap<QString, QStringList>& factors,
    const QStringList& valid_features
) {
    ComparePlotDialog dlg(factors, valid_features);

    if (dlg.is_accepted_) {
        return dlg.current_value();
    }
    else {
        return {};
    }
};

COMPARE_PLOT_ELEMENT
ComparePlotDialog::current_value() {

    QString feature = this->choose_feature_->text();
    bool normalize = this->normalize_->current_value() == SWITCH_ACCEPT;
    bool use_gene_activity = this->gene_activity_->current_value() == SWITCH_ACCEPT;
    bool show_p_value = this->p_value_->current_value() == SWITCH_ACCEPT;

    QString group = this->choose_group_->currentText();
    QString type = this->draw_type_->currentText();

    QList<QPair<QString, QString>> comparisons;
    int size = this->row_layouts_.size();
    for (int i = 0; i < size; ++i) {
        comparisons << qMakePair(this->comp1_[i]->currentText(), this->comp2_[i]->currentText());
    }

    return { feature, normalize, use_gene_activity, show_p_value, group, type, comparisons };
};

void ComparePlotDialog::s_factor_changed() {

    int size = this->row_layouts_.size();

    if (size < 1)return;

    for (int i = 0; i < size; ++i) {
        this->comparison_inner_layout_->removeItem(this->row_layouts_[i]);
        delete this->comp1_[i];
        delete this->comp2_[i];
        delete this->delete_buttons_[i];
        delete this->row_layouts_[i];
    }

    this->comp1_.clear();
    this->comp2_.clear();
    this->delete_buttons_.clear();
    this->row_layouts_.clear();
};

void ComparePlotDialog::s_delete_comparison() {

    QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());

    int index = this->delete_buttons_.indexOf(button);

    this->comparison_inner_layout_->removeItem(this->row_layouts_[index]);

    delete this->comp1_[index];
    delete this->comp2_[index];
    delete this->delete_buttons_[index];
    delete this->row_layouts_[index];

    this->comp1_.remove(index);
    this->comp2_.remove(index);
    this->delete_buttons_.remove(index);
    this->row_layouts_.remove(index);
};

void ComparePlotDialog::s_add_comparison() {

    G_SET_NEW_COMBOBOX(box1, this->factors_[this->choose_group_->currentText()], 30);
    G_SET_NEW_COMBOBOX(box2, this->factors_[this->choose_group_->currentText()], 30);
    this->comp1_ << box1;
    this->comp2_ << box2;

    G_SET_NEW_BUTTON(button, "-", QSize(30, 30));
    this->delete_buttons_ << button;
    connect(button, &QPushButton::clicked, this, &ComparePlotDialog::s_delete_comparison);

    G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(this->comparison_inner_layout_, row_layout, box1, box2, button);
    this->row_layouts_ << row_layout;
};

void ComparePlotDialog::s_change_interface() {

    QString ind = this->draw_type_->currentText();

    if (ind == "Custom") {
        if (this->is_expanding_) {
            return;
        }

        this->comparison_outer_layout_->addWidget(this->compare_widget_);
        this->compare_widget_->setVisible(true);
        this->is_expanding_ = true;
    }
    else {
        if (!this->is_expanding_) {
            return;
        }

        this->comparison_outer_layout_->removeWidget(this->compare_widget_);
        this->compare_widget_->setVisible(false);
        this->is_expanding_ = false;
    }
};

void ComparePlotDialog::accept() {
    QDialog::accept();
    this->is_accepted_ = true;
}

void ComparePlotDialog::reject() {
    QDialog::reject();
}