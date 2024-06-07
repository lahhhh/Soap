#include "VelocytoBaseItem.h"

#include "MatrixWindow.h"
#include "EnrichWorker.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CustomPlot.h"
#include "FileWritingWorker.h"

#include "VelocytoDownstreamWorker.h"
#include "ScveloWorker.h"

void VelocytoBaseItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(SparseInt));
	this->check_variable(DATA_SUBMODULES(VelocityEstimate));
	this->check_variable(DATA_SUBMODULES(ScveloEstimate));
}

void VelocytoBaseItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Estimate Velocity", s_estimate_velocity);
	ADD_MAIN_ACTION("Estimate Scvelo", s_estimate_scvelo);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
};

void VelocytoBaseItem::s_receive_estimate(VelocityEstimate* velocity_estimate) {
	QString title = this->signal_emitter_->get_unique_name(VARIABLE_VELOCITY_ESTIMATE);

	DATA_SUBMODULES(VelocityEstimate)[title] = std::move(*velocity_estimate);
	delete velocity_estimate;

	VelocityEstimateItem* item = new VelocityEstimateItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(VelocityEstimate)[title], 
		this->draw_suite_, 
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);
}

void VelocytoBaseItem::s_estimate_velocity() {

	G_GETLOCK;

	VelocytoDownstreamWorker* worker = new VelocytoDownstreamWorker(this->data());

	G_LINK_WORKER_THREAD(VelocytoDownstreamWorker, x_estimate_ready, VelocytoBaseItem, s_receive_estimate);
};

void VelocytoBaseItem::s_receive_scvelo(ScveloEstimate* scvelo_estimate) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_SCVELO_ESTIMATE);

	DATA_SUBMODULES(ScveloEstimate)[title] = std::move(*scvelo_estimate);
	delete scvelo_estimate;

	ScveloEstimateItem* item = new ScveloEstimateItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(ScveloEstimate)[title], 
		this->draw_suite_, 
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);
}

void VelocytoBaseItem::s_estimate_scvelo() {

	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		SingleCellRna* single_cell_rna = this->get_root<SingleCellRna>();

		auto pca = single_cell_rna->pca();
		if (pca == nullptr) {
			G_WARN("Please compute PCA before scvelo.");
			return;
		}

		G_GETLOCK;
		ScveloWorker* worker = new ScveloWorker(this->data(), pca->data_.mat_);
		G_LINK_WORKER_THREAD(ScveloWorker, x_scvelo_ready, VelocytoBaseItem, s_receive_scvelo);
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

		auto pca = single_cell_multiome->rna_pca();
		if (pca == nullptr) {
			G_WARN("Please compute RNA PCA before scvelo.");
			return;
		}

		G_GETLOCK;
		ScveloWorker* worker = new ScveloWorker(this->data(), pca->data_.mat_);
		G_LINK_WORKER_THREAD(ScveloWorker, x_scvelo_ready, VelocytoBaseItem, s_receive_scvelo);
	}
	else {
		G_WARN("You can only compute scvelo in a single-cell rna data or a single-cell multiome data.");
	}
};