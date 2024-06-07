#include "FragmentsItem.h"

#include "ItemIOWorker.h"
#include "CommonDialog.h"

#include "customplot.h"

void FragmentsItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->cell_names_.size()) + " ]");
}

void FragmentsItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}
