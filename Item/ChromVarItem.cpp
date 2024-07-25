#include "ChromVARItem.h"

#include "MatrixWindow.h"

#include "CommonDialog.h"
#include "DifferentialAnalysisWorker.h"

#include "SimpleUmapWorker.h"

void ChromVARItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Differential Analysis", s_differential_analysis);

	ADD_MAIN_MENU("Export");

	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

	ADD_MAIN_ACTION("Run UMAP", s_umap);

	ADD_MAIN_ACTION("Embedding Plot", s_embedding_plot);
}

void ChromVARItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(DifferentialAnalysis));
	this->check_variable(DATA_SUBMODULES(Embedding));
};

void ChromVARItem::s_differential_analysis() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Data is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();

	G_GETLOCK;

	auto factor_map = single_cell_multiome->metadata()->mat_.get_factor_information(false);
	if (factor_map.isEmpty()) {
		G_LOG("No suitable feature for Differential Analysis");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Differential Expression Analysis settings",
		{ "Choose Group:1", "P Adjust Method" },
		{soap::InputStyle::CompareLayout, soap::InputStyle::ComboBox},
		{{ "Bonferroni", "FDR" }},
		{factor_map}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto comparison = compare_layouts_to_list(settings[0]);
	QString factor_name = comparison[0], comparison1 = comparison[1], comparison2 = comparison[2];

	if (comparison1 == comparison2) {
		G_WARN("Illegal Comparison!");
		G_UNLOCK;
		return;
	}

	QStringList metadata = single_cell_multiome->metadata()->mat_.get_qstring(factor_name);

	QString p_adjust_method = settings[1];

	G_LOG("Differential Analysis in " + factor_name + "...");
	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		p_adjust_method
	);

	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, ChromVARItem, s_receive_differential_analysis);
};

void ChromVARItem::s_embedding_plot() {

	auto emb_names = this->index_tree_->search(soap::VariableType::Embedding);

	if (emb_names.isEmpty()) {
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Settings",
		{ "Feature", "Embedding" },
		{soap::InputStyle::LineEditWithCompleter, soap::InputStyle::ComboBox},
		{this->data()->motif_names_, emb_names}
	);
	if (settings.isEmpty()) {
		return;
	}

	auto index = this->data()->motif_names_.indexOf(settings[0]);
	if (index == -1) {
		G_WARN("Invalid Feature.");
		return;
	}

	auto item = (EmbeddingItem*)this->signal_emitter_->get_item(settings[1]);

	item->show(this->data()->z_.row(index), settings[0], false, "Z Score");
};

void ChromVARItem::s_receive_differential_analysis(DifferentialAnalysis da, QString name) {

	QString title = this->signal_emitter_->get_unique_name(name);
	DATA_SUBMODULES(DifferentialAnalysis)[title] = da;

	DifferentialAnalysisItem* item = new DifferentialAnalysisItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(DifferentialAnalysis)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Differential Analysis finished");
};

void ChromVARItem::s_umap() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This data is not part of Single Cell Multiome Data.");
		return;
	}

	G_GETLOCK;

	SimpleUmapWorker* worker = new SimpleUmapWorker(this->data()->z_.transpose());
	G_LINK_WORKER_THREAD(SimpleUmapWorker, x_umap_ready, ChromVARItem, s_receive_umap);
};

void ChromVARItem::s_receive_umap(Eigen::MatrixXd emb) {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		return;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto rna_counts = single_cell_multiome->rna_counts();

	if (rna_counts == nullptr) {
		return;
	}

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_CHROMVAR_UMAP);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Umap,
		emb,
		rna_counts->colnames_,
		_Cs paste("ChromVAR UMAP-", _Cs cast<QString>(_Cs seq_n(1, emb.cols()))));

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(Embedding)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("UMAP finished.");
};