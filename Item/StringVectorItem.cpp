#include "StringVectorItem.h"

#include "CustomMatrix.h"
#include "qCustomPlot.h"

#include "ItemIOWorker.h"
#include "FileWritingWorker.h"
#include "CommonDialog.h"
#include "MatrixWindow.h"
#include "TextEditWindow.h"


void StringVectorItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->data_.size()) + " ]");
};

void StringVectorItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->data_,
		this->title_,
		this->signal_emitter_
	);
}

void StringVectorItem::s_edit_manually() {
	G_GETLOCK;

	TextEditWindow::view(
		&this->data()->data_,
		this->data_,
		this->signal_emitter_
	);

	G_UNLOCK;
};

void StringVectorItem::s_edit_sort() {
	G_GETLOCK;

	this->data()->data_ = _Cs sorted(this->data()->data_);	

	G_UNLOCK;
};

void StringVectorItem::s_edit_remove_empty_elements() {

	G_GETLOCK;

	this->data()->data_.removeAll("");

	this->__s_update_interface();

	G_UNLOCK;
};

void StringVectorItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Edit");
	ADD_ACTION("Manually", "Edit", s_edit_manually);
	ADD_ACTION("Sort", "Edit", s_edit_sort);
	ADD_ACTION("Remove Empty Elements", "Edit", s_edit_remove_empty_elements);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};