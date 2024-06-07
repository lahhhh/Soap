#include "NumericMatrixItem.h"

#include "MatrixWindow.h"
#include "CustomMatrix.h"
#include "qCustomPlot.h"

void NumericMatrixItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void NumericMatrixItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->data_.rows()) + " | " + QString::number(this->data()->data_.cols()) + " ]");
};

void NumericMatrixItem::__show_this() {

	MatrixWindow::show_matrix(
		&this->data()->data_,
		this->title_,
		this->signal_emitter_
	);
}
