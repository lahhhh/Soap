#include "GenomicRangeItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "FileWritingWorker.h"
#include "ItemIOWorker.h"
#include "CalculateCountsByGenomicRangeWorker.h"

#include "SingleCellMultiomeItem.h"
#include "SingleCellAtacItem.h"

void GenomicRangeItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->size()) + " ]");
}

void GenomicRangeItem::__show_this() {
	MatrixWindow::show_matrix(
		this->data(), 
		custom::cast<QString>(custom::seq_n(1, this->data()->size())),
		this->data()->colnames(), 
		this->title_, 
		this->signal_emitter_,
		false, 
		this->data());
}

void GenomicRangeItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Change Sequence Name Style", s_change_sequence_name_style);

	ADD_MAIN_ACTION("Recalculate Counts by this Genomic Range", s_recalculate_counts);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void GenomicRangeItem::s_recalculate_counts() {

	G_GETLOCK;
	
	Fragments* fragments{ nullptr };

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

		fragments = single_cell_multiome->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments not loaded.");
			G_UNLOCK;
			return;
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);

		fragments = single_cell_atac->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments not loaded.");
			G_UNLOCK;
			return;
		}
	}
	else {
		G_WARN("Fragments Needed.");
		G_UNLOCK;
		return;
	}

	CalculateCountsByGenomicRangeWorker* worker = new CalculateCountsByGenomicRangeWorker(
		fragments, 
		*this->data()
	);

	G_LINK_WORKER_THREAD(CalculateCountsByGenomicRangeWorker, x_peak_counts_ready, GenomicRangeItem, s_recalculated_counts);
};

void GenomicRangeItem::s_recalculated_counts(SparseInt* data) {
	
	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

		auto single_cell_multiome_item = this->get_item<SingleCellMultiomeItem>(single_cell_multiome);

		auto atac_field = single_cell_multiome_item->get_atac_field();
		if (atac_field != nullptr) {
			atac_field->__remove_this();
		}

		atac_field = single_cell_multiome_item->create_field_item(DataField::DataType::Atac, "ATAC");

		QString title = this->signal_emitter_->get_unique_name(VARIABLE_ATAC_COUNTS);
		SUBMODULES(*atac_field->data(), SparseInt)[title] = std::move(*data);
		delete data;

		atac_field->__check_data();
		single_cell_multiome_item->atac_update_quality_control_information();

		G_LOG("Counts Updated.");
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);

		auto single_cell_atac_item = this->get_item<SingleCellAtacItem>(single_cell_atac);

		auto counts = single_cell_atac_item->counts();

		if (counts != nullptr) {
			counts->__remove_this();
		}

		auto normalized = single_cell_atac_item->normalized();

		if (normalized != nullptr) {
			normalized->__remove_this();
		}
		auto name = this->signal_emitter_->get_unique_name(VARIABLE_COUNTS);

		SUBMODULES(*single_cell_atac, SparseInt)[name] = std::move(*data);
		delete data;

		auto item = new SparseIntItem(
			name,
			single_cell_atac_item->index_tree_,
			single_cell_atac->counts(),
			single_cell_atac_item->draw_suite_,
			single_cell_atac_item->information_area_,
			single_cell_atac_item->signal_emitter_
		);

		single_cell_atac_item->set_item(item);

		single_cell_atac_item->update_quality_control_information();

		G_LOG("Counts Updated.");
	}
	else {
		delete data;
		return;
	}
};

QStringList GenomicRangeItem::ucsc_to_ensembl(const QStringList& ucsc_name) {
	return custom::sapply(ucsc_name, [](QString name)->QString {
		if (name == "chrM") {
			return "MT";
		}
		else {
			return name.remove("chr");
		}
	});
};

QStringList GenomicRangeItem::ensembl_to_ucsc(const QStringList& ensembl_name) {
	return custom::sapply(ensembl_name, [](QString name)->QString {
		if (name == "MT" || name == "Mt") {
			return "chrM";
		}
		else {
			return "chr" + name;
		}
	});
};

void GenomicRangeItem::s_change_sequence_name_style() {
	G_GETLOCK;
	static const QStringList style_list{ "UCSC", "ENSEMBL" };
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Settings for change name style",
		{ "From", "To" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox},
		{ style_list, style_list}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	if (settings[0] == settings[1]) {
		G_UNLOCK;
		return;
	}
	QString from = settings[0];
	QString to = settings[1];

	QStringList& sequence_names = this->data()->sequence_names_.data_;


	if (from == "UCSC" && to == "ENSEMBL") {
		sequence_names = ucsc_to_ensembl(sequence_names);
	}
	else if (from == "ENSEMBL" && to == "UCSC") {
		sequence_names = ensembl_to_ucsc(sequence_names);
	}

	G_UNLOCK;	
};