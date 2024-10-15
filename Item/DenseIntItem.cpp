#include "DenseIntItem.h"

#include "Custom.h"
#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CustomPlot.h"

void DenseIntItem::__s_update_interface() {

    this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " | " + QString::number(this->data()->mat_.cols()) + " ]");
};

void DenseIntItem::__show_this() {

    MatrixWindow::show_matrix(this->data(), this->title_, this->signal_emitter_, false, this->data_);
}

void DenseIntItem::__set_menu() {

    CREATE_ROOT_MENU;

    ADD_MAIN_ACTION("Rename", __s_rename);

    ADD_MAIN_MENU("Search");

    ADD_ACTION("by Row", "Search", s_search_row);
    ADD_ACTION("by Column", "Search", s_search_column);
    ADD_ACTION("by Block", "Search", s_search_block);

    ADD_MAIN_ACTION("Correlation Plot", s_correlation_plot);

    ADD_MAIN_MENU("Export");
    ADD_ACTION("as Item", "Export", __s_export_as_item);

    ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void DenseIntItem::update_shape() {

    this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " | " + QString::number(this->data()->mat_.cols()) + " ]");
};

void DenseIntItem::__s_delete_this() {

	if (this->attached_to(soap::VariableType::BulkRna)) {
		G_WARN("Count Item cannot be deleted from BulkRna.");
		return;
	}

	G_GETLOCK;

	G_UNLOCK;

	this->__remove_this();
};

void DenseIntItem::s_search_row() {
    QStringList quest = CommonDialog::get_response(
        this->signal_emitter_,
        "Enter Row Names Searched",
        { "Rownames" },
        { soap::InputStyle::MultipleLineEdit }
    );

    if (quest.isEmpty())return;

    quest = multiple_line_edit_to_list(quest[0]);

    QVector<int> row_index = custom::valid_index_of(quest, this->data()->rownames_);

    if (row_index.isEmpty()) {
        G_LOG("Feature Not Found.");
        return;
    }
    auto tmp = this->data()->row_reordered(row_index);
    MatrixWindow::show_matrix(
        &tmp,
        "Search Row from " + this->title_,
        this->signal_emitter_
    );
};

void DenseIntItem::s_search_column() {

    QStringList quest = CommonDialog::get_response(
        this->signal_emitter_,
        "Enter Column Names Searched",
        { "Colnames" },
        { soap::InputStyle::MultipleLineEdit }
    );

    if (quest.isEmpty())return;

    quest = multiple_line_edit_to_list(quest[0]);

    QVector<int> column_index = custom::valid_index_of(quest, this->data()->colnames_);

    if (column_index.isEmpty()) {
        G_LOG("Column Not Found.");
        return;
    }
    auto tmp = this->data()->col_reordered(column_index);
    MatrixWindow::show_matrix(
        &tmp,
        "Search Column from " + this->title_,
        this->signal_emitter_
    );
};

void DenseIntItem::s_search_block() {

    QStringList quest = CommonDialog::get_response(
        this->signal_emitter_,
        "Please enter names searched",
        { "Rownames", "Colnames" },
        { soap::InputStyle::MultipleLineEdit, soap::InputStyle::MultipleLineEdit }
    );
    if (quest.isEmpty())return;

    QStringList row_quests = multiple_line_edit_to_list(quest[0]);
    QStringList col_quests = multiple_line_edit_to_list(quest[1]);

    QVector<int> row_index = custom::valid_index_of(row_quests, this->data()->rownames_);
    QVector<int> column_index = custom::valid_index_of(col_quests, this->data()->colnames_);
    if (column_index.isEmpty() || row_index.isEmpty()) {
        G_LOG("Data Not Found.");
        return;
    }
    auto tmp = this->data()->reordered(row_index, column_index);
    MatrixWindow::show_matrix(
        &tmp,
        "Search Block from " + this->title_,
        this->signal_emitter_
    );
};

void DenseIntItem::s_correlation_plot() {

    auto&& valid_features = this->data()->rownames_;

    if (valid_features.size() < 2) {
        G_WARN("No Enough Data.");
        return;
    }

    auto settings = CommonDialog::get_response(
        this->signal_emitter_,
        "Plot Settings",
        { "Feature 1", "Feature 2" },
        { soap::InputStyle::LineEditWithCompleter, soap::InputStyle::LineEditWithCompleter },
        { valid_features, valid_features }
    );

    if (settings.isEmpty()) {
        return;
    }

    QString feature_1 = settings[0];
    QString feature_2 = settings[1];
    int index_1 = valid_features.indexOf(feature_1);
    int index_2 = valid_features.indexOf(feature_2);

    if (index_1 == -1 || index_2 == -1) {
        G_WARN("Invalid data.");
        return;
    }

    Eigen::ArrayXd f1 = this->data()->mat_.row(index_1).cast<double>();
    Eigen::ArrayXd f2 = this->data()->mat_.row(index_2).cast<double>();

    auto&& gs = this->draw_suite_->graph_settings_;
    auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);

    custom_plot::set_scatter_plot_axis_style(
        draw_area,
        axis_rect,
        feature_1,
        feature_2,
        f1,
        f2,
        gs
    );

    custom_plot::patch::scatter(
        draw_area,
        axis_rect,
        f1,
        f2,
        Qt::black,
        gs.get_scatter_point_size()
    );

    double cor = custom::correlation_pearson(f1, f2);

    custom_plot::add_title(draw_area, "Correlation : " + QString::number(cor), gs);

    this->draw_suite_->update(draw_area);

};