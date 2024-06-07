#include "CoverageTrackItem.h"

#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CoverageTrackWindow.h"
#include "GenomeUtility.h"
#include "CustomPlot.h"

void CoverageTrackItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->levels_.size()) + " ]");
};

void CoverageTrackItem::__set_menu() {
	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Export");

	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void CoverageTrackItem::__show_this() {

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

		CoverageTrackWindow::view(single_cell_multiome, this->data(), this->signal_emitter_);
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);

		CoverageTrackWindow::view(single_cell_atac, this->data(), this->signal_emitter_);
	}
	else {

		CoverageTrackWindow::view(this->data(), this->signal_emitter_);
	}
}
