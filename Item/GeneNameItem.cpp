#include "GeneNameItem.h"

#include "EnrichmentItem.h"

#include "CustomMatrix.h"
#include "CustomPlot.h"

#include "ItemIOWorker.h"
#include "FileWritingWorker.h"
#include "CommonDialog.h"
#include "MatrixWindow.h"
#include "TextEditWindow.h"

#include "EnrichWorker.h"
#include "GenomeUtility.h"
#include "StringVector.h"

void GeneNameItem::__s_update_interface() {
	
	this->setText(2, "[ " + QString::number(this->data()->data_.size()) + " ]");
};

void GeneNameItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->data_,
		this->title_,
		this->signal_emitter_
	);
}

void GeneNameItem::s_edit_manually() {
	G_GETLOCK;

	TextEditWindow::view(
		&this->data()->data_,
		this->data_,
		this->signal_emitter_
	);

	G_UNLOCK;
};

void GeneNameItem::s_edit_sort() {
	G_GETLOCK;

	this->data()->data_ = custom::sorted(this->data()->data_);

	G_UNLOCK;
};

void GeneNameItem::s_edit_unique() {

	G_GETLOCK;

	this->data()->data_ = custom::unique(this->data()->data_);

	this->__s_update_interface();

	G_UNLOCK;
};

void GeneNameItem::s_remove_elements_by_regular_expression() {
	G_GETLOCK;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Enter Regular Expression",
		{ "Re" },
		{soap::InputStyle::StringLineEdit}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QRegularExpression re(settings[0]);
	auto filter = custom::sapply(this->data()->data_,
		[&re](auto&& s) {return s.contains(re); }
	);
	this->data()->data_ = custom::sliced(this->data()->data_, custom::flip(filter));
	this->__s_update_interface();

	G_UNLOCK;
};

void GeneNameItem::s_remove_empty_elements() {
	G_GETLOCK;

	this->data()->data_.removeAll("");

	this->__s_update_interface();

	G_UNLOCK;
};

void GeneNameItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Edit");
	ADD_ACTION("Manually", "Edit", s_edit_manually);
	ADD_ACTION("Sort", "Edit", s_edit_sort);
	ADD_ACTION("Unique", "Edit", s_edit_unique);

	ADD_MENU("Edit | Remove Elements", "Remove Elements", "Edit");
	ADD_ACTION("Empty", "Edit | Remove Elements", s_remove_empty_elements);
	ADD_ACTION("by Regular Expression", "Edit | Remove Elements", s_remove_elements_by_regular_expression);
	ADD_ACTION("by Chromosome Location", "Edit | Remove Elements", s_remove_elements_by_chromosome_location);

	ADD_MAIN_MENU("Enrich");

	ADD_ACTION("GO", "Enrich", s_enrich_go);
	ADD_ACTION("KEGG", "Enrich", s_enrich_kegg);

	ADD_MAIN_MENU("Convert");

	ADD_ACTION("String Vector", "Convert", s_convert_to_string_vector);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void GeneNameItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Enrichment));
}


void GeneNameItem::s_enrich_kegg() {

	if (this->data()->data_.isEmpty()) {
		G_WARN("Empty Gene Name");
		return;
	}

	G_GETLOCK;

	soap::Species species;
	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		species = (this->get_root<SingleCellRna>())->species_;
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		species = (this->get_root<SingleCellMultiome>())->species_;
	}
	else {
		QStringList res = CommonDialog::get_response(
			this->signal_emitter_,
			"Species Setting",
			{ "Species" },
			{ soap::InputStyle::ComboBox},
			{ { "Human", "Mouse" }}
		);
		if (res.isEmpty()) {
			G_UNLOCK;
			return;
		}
		species = res[0] == "Human" ? soap::Species::Human : soap::Species::Mouse;
	}
	if (species == soap::Species::Undefined) {
		G_NOTICE("Only human and mouse gene enrichment are supported");
		G_UNLOCK;
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Enrichment Setting",
		{ "P adjust method", "P threshold (filter pathways):0.05"},
		{ soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit},
		{{ "FDR", "Bonferroni" }}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString p_adjust_method = settings[0];
	double pathway_p_threshold = settings[1].toDouble();
	if (pathway_p_threshold <= 0 || pathway_p_threshold > 1) {
		G_NOTICE("Invalid p value to filter pathways! Reset to 0.05");
		pathway_p_threshold = 0.05;
	}

	EnrichWorker* worker = new EnrichWorker(
		"KEGG", 
		"",
		this->data()->data_,
		"ALL",
		species, 
		p_adjust_method, 
		pathway_p_threshold
	);
	G_LINK_WORKER_THREAD(EnrichWorker, x_enrichment_ready, GeneNameItem, s_receive_enrichment)
};

void GeneNameItem::s_enrich_go() {

	if (this->data()->data_.isEmpty()) {
		G_WARN("Empty Gene Name");
		return;
	}

	G_GETLOCK;

	soap::Species species;
	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		species = (this->get_root<SingleCellRna>())->species_;
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		species = (this->get_root<SingleCellMultiome>())->species_;
	}
	else {
		QStringList res = CommonDialog::get_response(
			this->signal_emitter_,
			"Species Setting",
			{ "Species" },
			{ soap::InputStyle::ComboBox},
			{ { "Human", "Mouse" }}
		);
		if (res.isEmpty()) {
			G_UNLOCK;
			return;
		}
		species = res[0] == "Human" ? soap::Species::Human : soap::Species::Mouse;
	}
	if (species == soap::Species::Undefined) {
		G_NOTICE("Only human and mouse gene enrichment are supported");
		G_UNLOCK;
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Enrichment Setting",
		{"Ontology", "P adjust method", "P threshold (pathway):0.05"},
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit},
		{ { "BP", "MF", "CC", "ALL"}, {"FDR", "Bonferroni"} }
	);
	
	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString ontology = settings[0], p_adjust_method = settings[1];

	double pathway_p_threshold = settings[2].toDouble();
	if (pathway_p_threshold <= 0 || pathway_p_threshold > 1) {
		G_NOTICE("Invalid p value to filter pathways! Reset to 0.05");
		pathway_p_threshold = 0.05;
	}

	EnrichWorker* worker = new EnrichWorker(
		"GO", 
		"",
		this->data()->data_, 
		ontology, 
		species, 
		p_adjust_method, 
		pathway_p_threshold
	);
	G_LINK_WORKER_THREAD(EnrichWorker, x_enrichment_ready, GeneNameItem, s_receive_enrichment);
};

void GeneNameItem::s_receive_enrichment(const CustomMatrix& matrix, QString name) {

	QString new_title = this->signal_emitter_->get_unique_name(name);
	DATA_SUBMODULES(Enrichment)[new_title] = Enrichment(matrix);

	EnrichmentItem* item = new EnrichmentItem(
		new_title, 
		this->index_tree_, 
		&DATA_SUBMODULES(Enrichment)[new_title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);
};

void GeneNameItem::s_remove_elements_by_chromosome_location() {

	G_GETLOCK;
	G_UNLOCK;

	G_NOTICE("Now only human genome is supported (default).");

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Location Settings",
		{ "Location" },
		{ soap::InputStyle::MultipleLineEdit }
	);

	if (settings.isEmpty()) {
		return;
	}

	auto locs = multiple_line_edit_to_list(settings[0]);

	if (locs.isEmpty()) {
		return;
	}

	auto hg38 = custom::get_hg38_gene_location();

	int n_loc = locs.size();

	QStringList to_remove;

	for (int i = 0; i < n_loc; ++i) {

		to_remove << hg38.find_gene(locs[i]);
		to_remove << hg38.find_chromosome_gene(locs[i]);
	}

	to_remove = custom::unique(to_remove);

	QStringList d = this->data()->data_;

	d = custom::sliced(d, custom::flip(custom::in(d, to_remove)));

	this->data()->data_ = d;

	this->__s_update_interface();
};

void GeneNameItem::s_convert_to_string_vector() {

	StringVector* sv = new StringVector(this->data()->data_);

	this->signal_emitter_->x_data_create_soon(sv, soap::VariableType::StringVector, "Converted From " + this->title_);
};