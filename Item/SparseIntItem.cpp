#include "SparseIntItem.h"

#include "Custom.h"
#include "MatrixWindow.h"
#include "ItemIOWorker.h"
#include "CommonDialog.h"

#include "SingleCellRnaCreateWorker.h"

void SparseIntItem::__show_this() {

	MatrixWindow::show_matrix(this->data(), this->title_, this->signal_emitter_, false, this->data_);
}

void SparseIntItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Search");

	ADD_ACTION("by Row", "Search", s_search_row);
	ADD_ACTION("by Column", "Search", s_search_column);
	ADD_ACTION("by Block", "Search", s_search_block);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Extract as SingleCellRna Data", s_extract_as_single_cell_rna);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

}

void SparseIntItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " | " + QString::number(this->data()->mat_.cols()) + " ]");
};


void SparseIntItem::s_search_row() {
	QStringList quest = CommonDialog::get_response(
		this->signal_emitter_,
		"Enter Row Names Searched",
		{ "Rownames" },
		{ soap::InputStyle::MultipleLineEdit}
	);

	if (quest.isEmpty())return;

	quest = multiple_line_edit_to_list(quest[0]);

	QVector<int> row_index = custom::valid_index_of(quest, this->data()->rownames_);

	if (row_index.isEmpty()) {
		G_NOTICE("Feature Not Found.");
		return;
	}
	auto tmp = this->data()->row_reordered(row_index);
	MatrixWindow::show_matrix(
		&tmp, 
		"Search Row from " + this->title_,
		this->signal_emitter_);
};


void SparseIntItem::s_search_column() {
	QStringList quest = CommonDialog::get_response(
		this->signal_emitter_,
		"Enter Column Names Searched",
		{ "Colnames" },
		{ soap::InputStyle::MultipleLineEdit}
	);

	if (quest.isEmpty())return;

	quest = multiple_line_edit_to_list(quest[0]);

	QVector<int> column_index = custom::valid_index_of(quest, this->data()->colnames_);

	if (column_index.isEmpty()) {
		G_NOTICE("Column Not Found.");
		return;
	}
	auto tmp = this->data()->col_reordered(column_index);
	MatrixWindow::show_matrix(
		&tmp, 
		"Search Column from " + this->title_, 
		this->signal_emitter_);
};

void SparseIntItem::s_search_block() {
	QStringList quest = CommonDialog::get_response(
		this->signal_emitter_,
		"Please enter names searched",
		{ "Rownames", "Colnames" },
		{ soap::InputStyle::MultipleLineEdit, soap::InputStyle::MultipleLineEdit}
	);
	if (quest.isEmpty())return;

	QStringList row_quests = multiple_line_edit_to_list(quest[0]);
	QStringList col_quests = multiple_line_edit_to_list(quest[1]);

	QVector<int> row_index = custom::valid_index_of(row_quests, this->data()->rownames_);
	QVector<int> column_index = custom::valid_index_of(col_quests, this->data()->colnames_);
	if (column_index.isEmpty() || row_index.isEmpty()) {
		G_NOTICE("Data Not Found.");
		return;
	}
	auto tmp = this->data()->reordered(row_index, column_index);
	MatrixWindow::show_matrix(
		&tmp, 
		"Search Block from " + this->title_, 
		this->signal_emitter_);
};

void SparseIntItem::__s_delete_this() {
	if (this->attached_to(soap::VariableType::SingleCellRna)) {
		G_WARN("Count Item cannot be deleted from SingleCellRna.");
		return;
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
		G_WARN("Count Item cannot be deleted from SingleCellAtac.");
		return;
	}
	else if (this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("Count Item cannot be deleted from SingleCellMultiome.");
		return;
	}
	else if (this->attached_to(soap::VariableType::VelocytoBase)) {
		G_WARN("Count Item cannot be deleted from VelocytoBase. Please delete the whole VelocytoBase directly.");
		return;
	}

	G_GETLOCK;

	G_UNLOCK;

	this->__remove_this();
};

void SparseIntItem::s_receive_single_cell_rna(SingleCellRna* data) {

	this->signal_emitter_->x_data_create_soon(data, soap::VariableType::SingleCellRna, "Created SingleCellRna");
}

void SparseIntItem::s_extract_as_single_cell_rna() {

	G_GETLOCK;

	auto worker = new SingleCellRnaCreateWorker(this->data());

	G_LINK_WORKER_THREAD(SingleCellRnaCreateWorker, x_single_cell_rna_created, SparseIntItem, s_receive_single_cell_rna);

};