#include "SingleCellMultiomeItem.h"

#include "Custom.h"
#include "CustomPlot.h"
#include "FileIO.h"

#include "YesOrNoDialog.h"
#include "CommonDialog.h"

#include "ItemIOWorker.h"
#include "PcaWorker.h"
#include "UmapWorker.h"
#include "ScrubletWorker.h"
#include "DifferentialAnalysisWorker.h"
#include "IntegrateWorker.h"
#include "GseaWorker.h"
#include "HarmonyWorker.h"
#include "ScicnvWorker.h"
#include "InferCnvWorker.h"
#include "MacsCallPeakWorker.h"
#include "MotifLocateWorker.h"
#include "CreateCoverageTrackWorker.h"
#include "LoadFragmentsWorker.h"
#include "FragmentsQualityViewWorker.h"
#include "TssPlotWorker.h"
#include "StatisticsDialog.h"
#include "CalculateGeneActivityWorker.h"
#include "VelocytoWorker.h"
#include "PandoWorker.h"
#include "CellTypeAnnotationWorker.h"
#include "Monocle3Worker.h"
#include "CiceroWorker.h"
#include "ScentWorker.h"

DataFieldItem* SingleCellMultiomeItem::create_field_item(DataField::DataType type, const QString& name) {

	auto title = this->signal_emitter_->get_unique_name(name);
	auto data = &DATA_SUBMODULES(DataField)[title];
	data->data_type_ = type;

	auto item = new DataFieldItem(
		title,
		this->index_tree_,
		data,
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	return item;
};

DataFieldItem* SingleCellMultiomeItem::get_field_item(DataField::DataType type) {

	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::DataField) {
			DataField* data_field = static_cast<DataField*>(item->data_);
			if (data_field->data_type_ == type) {
				return static_cast<DataFieldItem*>(item);
			}
		}
	}

	return nullptr;
};

DataFieldItem* SingleCellMultiomeItem::get_rna_field() {

	return this->get_field_item(DataField::DataType::Rna);
};

DataFieldItem* SingleCellMultiomeItem::get_atac_field() {

	return this->get_field_item(DataField::DataType::Atac);
};

DataFieldItem* SingleCellMultiomeItem::get_trans_field() {
	
	return this->get_field_item(DataField::DataType::Trans);
};

void SingleCellMultiomeItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Metadata));
	this->check_variable(DATA_SUBMODULES(DataField));
	this->check_variable(DATA_SUBMODULES(GenomicRange));
	this->check_variable(DATA_SUBMODULES(GSEA));
	this->check_variable(DATA_SUBMODULES(CellChat));
	this->check_variable(DATA_SUBMODULES(CNV));
	this->check_variable(DATA_SUBMODULES(MotifPosition));	
	this->check_variable(DATA_SUBMODULES(CoverageTrack));
	this->check_variable(DATA_SUBMODULES(Fragments));
	this->check_variable(DATA_SUBMODULES(VelocytoBase));
	this->check_variable(DATA_SUBMODULES(Pando));
	this->check_variable(DATA_SUBMODULES(Monocle3));
	this->check_variable(DATA_SUBMODULES(Cicero));

	auto item = new NoteItem(&this->data()->string_information_["Note"], this->data(), this->signal_emitter_);
	this->addChild(item);
};

void SingleCellMultiomeItem::rna_update_quality_control_information() {

	auto rna_counts = this->data()->rna_counts();
	if (rna_counts == nullptr) {
		return;
	}

	Eigen::ArrayXi rna_count = _Cs col_sum_mt(rna_counts->mat_);
	const int ncol_rna = rna_counts->mat_.cols();
	Eigen::ArrayXi rna_gene(ncol_rna);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(ncol_rna);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(ncol_rna);

	for (std::size_t i = 0; i < ncol_rna; ++i) {
		rna_gene[i] = rna_counts->mat_.outerIndexPtr()[i + 1] - rna_counts->mat_.outerIndexPtr()[i];
	}

	QVector<int> mitochondrial_location, ribosomal_location;

	auto& gene_symbols = rna_counts->rownames_;
	if (this->data()->species_ == soap::Species::Human) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (this->data()->species_ == soap::Species::Mouse) {
		for (int i = 0; i < gene_symbols.length(); ++i) {
			if (gene_symbols.at(i).startsWith("mt-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("Rps") || gene_symbols.at(i).startsWith("Rpl")) {
				ribosomal_location << i;
			}
		}
	}

	if (mitochondrial_location.length() > 0) {
		for (int& i : mitochondrial_location) {
			mitochondrial_content += rna_counts->mat_.row(i).cast<double>();
		}
		mitochondrial_content /= rna_count.cast<double>();
	}
	if (ribosomal_location.length() > 0) {
		for (int& i : ribosomal_location) {
			ribosomal_content += rna_counts->mat_.row(i).cast<double>();
		}
		ribosomal_content /= rna_count.cast<double>();
	}
	_Cs remove_na(mitochondrial_content);
	_Cs remove_na(ribosomal_content);

	auto&& metadata = *this->data()->metadata();
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, _Cs cast<QVector>(rna_count));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, _Cs cast<QVector>(rna_gene));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, _Cs cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, _Cs cast<QVector>(ribosomal_content));
}

void SingleCellMultiomeItem::atac_update_quality_control_information() {

	auto atac_counts = this->data()->atac_counts();
	if (atac_counts == nullptr) {
		return;
	}

	Eigen::ArrayXi peak_count = _Cs col_sum_mt(atac_counts->mat_);
	const int ncol_atac = atac_counts->mat_.cols();
	Eigen::ArrayXi peak_gene(ncol_atac);

	for (std::size_t i = 0; i < ncol_atac; ++i) {
		peak_gene[i] = atac_counts->mat_.outerIndexPtr()[i + 1] - atac_counts->mat_.outerIndexPtr()[i];
	}

	Metadata& metadata = *this->data()->metadata();

	metadata.mat_.update(METADATA_ATAC_UMI_NUMBER, _Cs cast<QVector>(peak_count));
	metadata.mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, _Cs cast<QVector>(peak_gene));
}

void SingleCellMultiomeItem::update_quality_control_information(DataField::DataType type) {

	switch (type)
	{
	case DataField::DataType::Rna:
		this->rna_update_quality_control_information();
		break;
	case DataField::DataType::Atac:
		this->atac_update_quality_control_information();
		break;
	default:
		break;
	}
};

void SingleCellMultiomeItem::__set_menu() {
	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Statistics", s_statistics);

	//field
	ADD_MAIN_MENU("RNA");

	//fragments
	ADD_MAIN_MENU("Fragments");

	ADD_ACTION("Load Fragments", "Fragments", s_load_fragments);

	//coverage

	ADD_MAIN_MENU("Coverage");

	ADD_ACTION( "Coverage Plot", "Coverage", s_coverage_plot);
	ADD_ACTION( "Create Coverage Track", "Coverage", s_create_coverage_track);

	//sample
	ADD_MAIN_ACTION("Sample", s_sample);

	//setting
	ADD_MAIN_MENU("Settings");

	ADD_ACTION("Random Seed", "Settings", s_set_random_state);
	ADD_ACTION("Species", "Settings", s_set_species);

	// annotate
	ADD_MAIN_MENU("Annotate");

	ADD_MENU("Annotate | Add Metadata", "Add Metadata", "Annotate");

	ADD_ACTION("Add New Metadata", 
		"Annotate | Add Metadata", 
		s_add_new_metadata);

	ADD_ACTION("Combine Existed Metadata",
		"Annotate | Add Metadata",
		s_combine_existed_metadata);

	ADD_ACTION("Edit Metadata", "Annotate", s_edit_metadata);

	ADD_ACTION("Annotate Cell Type", "Annotate", s_fast_annotation);


	// qc
	ADD_MAIN_MENU("Quality Control");

	ADD_ACTION("Calculate Quality Metrics", "Quality Control", s_calculate_quality_metrics);

	ADD_ACTION("Show Quality Metrics", "Quality Control", s_show_quality_matrics);

	ADD_ACTION("Show Fragments Length Distribution", "Quality Control", s_show_fragments_length_distribution);

	ADD_ACTION( "TSS Plot", "Quality Control", s_tss_plot);

	// calculate gene activity
	ADD_MAIN_ACTION("Calculate Gene Activity", s_calculate_gene_activity);

	// call peak
	ADD_MAIN_MENU("Call Peaks");

	ADD_ACTION("MACS", "Call Peaks", s_call_peaks_by_macs);

	// find motif

	ADD_MAIN_ACTION("Find Motifs", s_find_motifs);

	// cicero

	ADD_MAIN_ACTION("Cicero", s_cicero);

	//Doublet Detection
	ADD_MAIN_MENU("Detect Doublets");

	ADD_MENU("Detect Doublets | Scrublet", "Scrublet", "Detect Doublets");

	ADD_ACTION("Run Scrublet", "Detect Doublets | Scrublet", s_rna_scrublet);
	ADD_ACTION("Set Threshold", "Detect Doublets | Scrublet", s_set_scrublet_threshold);

	//filter
	ADD_MAIN_MENU("Filter");

	ADD_ACTION("By Features", "Filter", s_filter_by_features);
	ADD_ACTION("By Parameters", "Filter", s_filter_by_parameters);
	ADD_ACTION("By Quality Metrics", "Filter", s_filter_by_quality_metrics);

	// Lineage Tracing
	ADD_MAIN_MENU("Lineage Tracing");

	ADD_ACTION("Monocle3", "Lineage Tracing", s_monocle3);
	ADD_ACTION("Velocyto", "Lineage Tracing", s_velocyto);


	// GSEA
	ADD_MENU("RNA | GSEA", "GSEA", "RNA");
	ADD_ACTION("in database", "RNA | GSEA", s_gsea_in_database);
	//ADD_ACTION("from input", "RNA | GSEA", s_gsea_from_input);

	// CNV
	ADD_MENU("RNA | CNV Detection", "CNV Detection", "RNA");

	ADD_ACTION("SciCnv", "RNA | CNV Detection", s_scicnv);

	ADD_ACTION("InferCNV [No HMM]", "RNA | CNV Detection", s_infercnv);

	ADD_MAIN_ACTION("Infer Gene Regulatory NetWork", s_pando);

	// landscape plot
	ADD_MAIN_ACTION("ATAC Landscape Plot", s_show_atac_landscape);

	// integrate
	ADD_MAIN_ACTION("Integrate With...", s_integrate);

	// duplicate
	ADD_MAIN_ACTION("Duplicate", __s_duplicate);

	// save
	ADD_MAIN_ACTION("Save", __s_export_as_item);

	// delete
	ADD_MAIN_ACTION("Delete", __s_delete_this);

	ADD_MAIN_MENU("More...");

	ADD_ACTION("Calculate Signalling Entropy[SCENT]", "More...", s_scent);
}

void SingleCellMultiomeItem::s_statistics() {
	StatisticsDialog::get_response(this->data());
};

void SingleCellMultiomeItem::s_load_fragments() {
	G_GETLOCK;

	QStringList fragments_file_path = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Fragments File",
		{ "Fragments (tsv.gz):Zipped TSV (*.tsv.gz)" },
		{ soap::InputStyle::MultiFile},
		{ {}}
	);

	if (fragments_file_path.isEmpty()) {
		G_UNLOCK;
		return;
	}
	fragments_file_path = multiple_file_to_list(fragments_file_path[0]);

	if (fragments_file_path.isEmpty()) {
		G_UNLOCK;
		return;
	}

	if (auto item = this->fragments()) {
		item->__remove_this();
	}

	QStringList barcodes;
	auto& metadata = this->data()->metadata()->mat_;
	if (metadata.contains(METADATA_BARCODES) && metadata.data_type_[METADATA_BARCODES] == CustomMatrix::DataType::QString) {
		barcodes = metadata.get_const_qstring_reference(METADATA_BARCODES);
	}
	else {
		G_NOTICE("Barcodes data is not found. Using cell names instead.");
		barcodes = this->data()->rna_counts()->colnames_;
	}
	LoadFragmentsWorker* worker = new LoadFragmentsWorker(barcodes, fragments_file_path);
	G_LINK_WORKER_THREAD(LoadFragmentsWorker, x_fragments_ready, SingleCellMultiomeItem, s_receive_loaded_fragments);
};

void SingleCellMultiomeItem::s_receive_macs_peaks(GenomicRange genomic_range) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_MACS_PEAKS);
	DATA_SUBMODULES(GenomicRange)[title] = std::move(genomic_range);

	GenomicRangeItem* item = new GenomicRangeItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(GenomicRange)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("MACS finished");
};

void SingleCellMultiomeItem::s_receive_loaded_fragments(Fragments* fragments) {

	fragments->cell_names_ = this->data()->rna_counts()->colnames_;

	auto title = this->signal_emitter_->get_unique_name(VARIABLE_FRAGMENTS);

	DATA_SUBMODULES(Fragments)[title] = std::move(*fragments);
	delete fragments;

	FragmentsItem* item = new FragmentsItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(Fragments)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Fragments : " + title + " is created.");
};

void SingleCellMultiomeItem::s_receive_coverage_track(CoverageTrack* d) {

	auto title = this->signal_emitter_->get_unique_name(VARIABLE_COVERAGE_TRACK);

	DATA_SUBMODULES(CoverageTrack)[title] = std::move(*d);
	delete d;

	CoverageTrackItem* item = new CoverageTrackItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(CoverageTrack)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("CoverageTrack : " + title + " is created.");
};

void SingleCellMultiomeItem::s_create_coverage_track() {

	auto fragments = this->data()->fragments();

	if (fragments == nullptr) {
		G_WARN("No Fragments loaded.");
		return;
	}

	G_GETLOCK;

	const auto& metadata = this->data()->metadata()->mat_;
	QStringList group_names = metadata.get_factor_name();
	if (group_names.isEmpty()) {
		G_WARN("No Suitable Metadata for Creating Coverage Track.");
		G_UNLOCK;
		return;
	}
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Coverage Track Settings",
		{ "Group" },
		{ soap::InputStyle::ComboBox},
		{ group_names}
	);
	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	CreateCoverageTrackWorker* worker = new CreateCoverageTrackWorker(
		this->data()->metadata(),
		fragments, settings[0]
	);
	G_LINK_WORKER_THREAD(CreateCoverageTrackWorker, x_coverage_track_ready, SingleCellMultiomeItem, s_receive_coverage_track);
};

void SingleCellMultiomeItem::s_tss_plot() {

	if(this->data()->atac_counts() == nullptr){
		G_WARN("No ATAC Data Found.");
		return;
	}

	auto fragments = this->data()->fragments();
	if (fragments == nullptr) {
		G_WARN("No Fragments Loaded.");
		return;
	}

	G_GETLOCK;

	G_LOG("Start TSS Estimation from Fragments.");

	TssPlotWorker* worker = new TssPlotWorker(this->data()->species_, fragments);
	G_LINK_WORKER_THREAD(TssPlotWorker, x_tss_ready, SingleCellMultiomeItem, s_receive_tss_plot_data);
};

void SingleCellMultiomeItem::s_receive_tss_plot_data(Eigen::ArrayXXd tss_matrix, QVector<double> tss_vector) {
	const auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, main_layout] = _Cp prepare_lg(gs);

	QCPMarginGroup* vertical_margin = new QCPMarginGroup(draw_area);
	QCPMarginGroup* horizontal_margin = new QCPMarginGroup(draw_area);

	QCPAxisRect* tss_line_rect = new QCPAxisRect(draw_area);
	main_layout->addElement(0, 0, tss_line_rect);

	tss_line_rect->setMarginGroup(QCP::msLeft | QCP::msRight, vertical_margin);

	QPen border_pen(Qt::black);
	border_pen.setWidth(2);

	tss_line_rect->setupFullAxesBox();

	tss_line_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
	tss_line_rect->axis(QCPAxis::atBottom)->setTicks(true);
	tss_line_rect->axis(QCPAxis::atBottom)->setSubTickPen(Qt::NoPen);
	tss_line_rect->axis(QCPAxis::atBottom)->setBasePen(border_pen);

	tss_line_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	tss_line_rect->axis(QCPAxis::atLeft)->setTicks(true);
	tss_line_rect->axis(QCPAxis::atLeft)->setSubTickPen(Qt::NoPen);
	tss_line_rect->axis(QCPAxis::atLeft)->setBasePen(border_pen);

	tss_line_rect->axis(QCPAxis::atTop)->grid()->setVisible(false);
	tss_line_rect->axis(QCPAxis::atTop)->setTicks(false);
	tss_line_rect->axis(QCPAxis::atTop)->setBasePen(border_pen);

	tss_line_rect->axis(QCPAxis::atRight)->grid()->setVisible(false);
	tss_line_rect->axis(QCPAxis::atRight)->setTicks(false);
	tss_line_rect->axis(QCPAxis::atRight)->setBasePen(border_pen);

	int n_tss = tss_vector.size();

	_CpPatch line(draw_area, tss_line_rect, _Cs linspaced(n_tss, 0, n_tss - 1), tss_vector, Qt::red, 3);
	_CpPatch set_range(tss_line_rect, QCPRange(0, n_tss - 1), QCPRange(0, 1.1 * std::ranges::max(tss_vector)));

	Eigen::ArrayXd bottom_axis_location(3);
	bottom_axis_location[0] = 0;
	bottom_axis_location[1] = (n_tss - 1) / 2.0;
	bottom_axis_location[2] = n_tss - 1;

	_Cp set_bottom_axis_label(
		tss_line_rect, 
		bottom_axis_location, 
		{ "-3kb" , "Start Site" , "3kb" }, 
		6,
		gs
	);

	QCPAxisRect* color_rect = new QCPAxisRect(draw_area);

	QCPColorMap* heatmap = new QCPColorMap(color_rect->axis(QCPAxis::atBottom), color_rect->axis(QCPAxis::atLeft));
	color_rect->setMarginGroup(QCP::msLeft | QCP::msRight, vertical_margin);
	color_rect->setMarginGroup(QCP::msTop | QCP::msBottom, horizontal_margin);
	const int n_row = tss_matrix.rows(), n_col = tss_matrix.cols();

	heatmap->data()->setSize(n_col, n_row);
	heatmap->data()->setRange(QCPRange(0, n_col - 1), QCPRange(0, n_row - 1));
	_CpPatch remove_left_bottom_axis(color_rect);
	_CpPatch set_range(color_rect, QCPRange(-0.5, n_col - 0.5), QCPRange(-0.5, n_row - 0.5));

	for (int i = 0; i < n_row; ++i) {
		for (int j = 0; j < n_col; ++j) {
			heatmap->data()->setCell(j, i, tss_matrix(i, j));
		}
	}
	heatmap->setInterpolate(false);
	heatmap->setTightBoundary(false);
	main_layout->addElement(1, 0, color_rect);

	QCPColorScale* color_scale = new QCPColorScale(draw_area);
	color_scale->setType(QCPAxis::atRight);
	main_layout->addElement(1, 1, color_scale);
	heatmap->setColorScale(color_scale);

	color_scale->axis()->setSubTickPen(Qt::NoPen);

	QCPColorGradient gradient;
	gradient.setColorStopAt(0.0, QColor(_CpColor firebrick3));
	gradient.setColorStopAt(0.5, QColor(_CpColor gold));
	gradient.setColorStopAt(1.0, QColor(_CpColor navy));
	heatmap->setGradient(gradient);
	heatmap->rescaleDataRange();

	color_scale->setMarginGroup(QCP::msTop | QCP::msBottom, horizontal_margin);

	_Cp add_title(draw_area, "TSS Score", gs);

	this->draw_suite_->update(draw_area);
};

void SingleCellMultiomeItem::s_coverage_plot() {

	auto fragments = this->data()->fragments();
	if (fragments == nullptr) {
		G_WARN("No Fragments Loaded.");
		return;
	}

	const auto& metadata = this->data()->metadata()->mat_;
	auto group_info = metadata.get_factor_information();
	if (group_info.isEmpty()) {
		G_WARN("No Suitable Metadata for Visualization.");
		return;
	}
	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Coverage Plot Settings",
		{ "Group", "Region", "Draw Gene:yes", "Draw Link:no", "Link Cutoff:0.5", "Draw Legend:no"},
		{ soap::InputStyle::FactorChoice, soap::InputStyle::StringLineEdit, soap::InputStyle::SwitchButton,
		soap::InputStyle::SwitchButton, soap::InputStyle::NumericLineEdit, soap::InputStyle::SwitchButton},
		{}, {group_info}
	);
	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto [name, levels] = factor_choice_to_pair(settings[0]);

	if (levels.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto factors = metadata.get_qstring(name);

	QString region = settings[1];

	bool draw_gene = switch_to_bool(settings[2]);
	bool draw_link = switch_to_bool(settings[3]);
	double link_cutoff = settings[4].toDouble();
	bool draw_legend = switch_to_bool(settings[5]);

	CoveragePlotWorker* worker = new CoveragePlotWorker(
		this->data(), 
		fragments, 
		factors,
		levels,
		region, 
		draw_gene, 
		draw_link, 
		link_cutoff, 
		draw_legend
	);
	G_LINK_WORKER_THREAD(CoveragePlotWorker, x_plot_ready, SingleCellMultiomeItem, s_receive_coverage_plot_data);
};

void SingleCellMultiomeItem::s_receive_coverage_plot_data(COVERAGE_PLOT_ELEMENTS res) {

	auto& gs = this->draw_suite_->graph_settings_;

	auto draw_area = _Cp coverage_plot(res, gs);

	this->draw_suite_->update(draw_area);
};

void SingleCellMultiomeItem::s_receive_motif_location(MotifPosition motif_position) {

	QString title = this->signal_emitter_->get_unique_name("Motif Position");
	DATA_SUBMODULES(MotifPosition)[title] = std::move(motif_position);

	MotifPositionItem* item = new MotifPositionItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(MotifPosition)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Motif Finding finished");
};

void SingleCellMultiomeItem::s_receive_fast_annotation(
	QStringList main_type,
	QStringList sub_type,
	QString main_type_name,
	QString sub_type_name
) {

	this->data()->metadata()->mat_.update(main_type_name, main_type, CustomMatrix::DataType::QStringFactor);
	this->data()->metadata()->mat_.update(sub_type_name, sub_type, CustomMatrix::DataType::QStringFactor);

	this->signal_emitter_->x_update_interface();

	G_NOTICE("Cell Type Annotation Finished.");
};

void SingleCellMultiomeItem::s_fast_annotation() {
	G_GETLOCK;

	auto rna_counts = this->data()->rna_counts();
	if (rna_counts == nullptr) {
		G_WARN("No RNA Counts Data Found.");
		G_UNLOCK;
		return;
	}

	if (this->data()->species_ != soap::Species::Human) {
		G_WARN("Fast Annotation Now Only Support Human");
		G_UNLOCK;
		return;
	}

	auto cluster_info = this->data()->metadata()->mat_.get_factor_name(false);

	CellTypeAnnotationWorker* worker;
	if (cluster_info.isEmpty()) {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Set Name",
			{ "Cell Type Name:Celltype-Soap", "Sub Type Name:Subtype-Soap" },
			{soap::InputStyle::StringLineEdit, soap::InputStyle::StringLineEdit}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		QString main_type_name = settings[0];
		QString sub_type_name = settings[1];

		worker = new CellTypeAnnotationWorker(
			rna_counts,
			main_type_name,
			sub_type_name,
			false
		);

	}
	else {

		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Set Name",
			{ "Cell Type Name:Celltype-Soap", "Sub Type Name:Subtype-Soap", "Annotation By Cluster:yes", "Annotate By: "},
			{soap::InputStyle::StringLineEdit, soap::InputStyle::StringLineEdit,
			soap::InputStyle::SwitchButton, soap::InputStyle::ComboBox},
			{cluster_info}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		QString main_type_name = settings[0];
		QString sub_type_name = settings[1];

		bool annotation_by_cluster = switch_to_bool(settings[2]);

		if (annotation_by_cluster) {

			worker = new CellTypeAnnotationWorker(
				rna_counts,
				main_type_name,
				sub_type_name,
				true,
				this->data()->metadata()->mat_.get_qstring(settings[3])
			);
		}
		else {

			worker = new CellTypeAnnotationWorker(
				rna_counts,
				main_type_name,
				sub_type_name,
				false
			);
		}
		
	}

	G_LINK_WORKER_THREAD(CellTypeAnnotationWorker, x_annotation_ready, SingleCellMultiomeItem, s_receive_fast_annotation);
};

void SingleCellMultiomeItem::s_cicero() {

	G_GETLOCK;

	auto atac_counts = this->data()->atac_counts();
	if (atac_counts == nullptr) {
		G_WARN("No Peak Data Found.");
		G_UNLOCK;
		return;
	}

	auto embedding_names = this->index_tree_->search(soap::VariableType::Embedding);
	if (embedding_names.isEmpty()) {
		G_WARN("No Embedding.");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Cicero Settings",
		{ "Embedding" },
		{soap::InputStyle::ComboBox},
		{embedding_names}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString embedding_name = settings[0];
	Embedding* embedding = static_cast<Embedding*>(this->index_tree_->search(embedding_name)->data_);

	CiceroWorker* worker = new CiceroWorker(atac_counts, embedding);
	G_LINK_WORKER_THREAD(CiceroWorker, x_cicero_ready, SingleCellMultiomeItem, s_receive_cicero);
};

void SingleCellMultiomeItem::s_receive_cicero(Cicero* cicero) {

	QString title = this->signal_emitter_->get_unique_name("Cicero");
	DATA_SUBMODULES(Cicero)[title] = std::move(*cicero);
	delete cicero;

	if (auto item = this->cicero()) {
		item->__remove_this();
	}

	CiceroItem* item = new CiceroItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(Cicero)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Cicero finished");
};

void SingleCellMultiomeItem::s_find_motifs() {
	
	G_GETLOCK;

	auto atac_counts = this->data()->atac_counts();
	if (atac_counts == nullptr) {
		G_WARN("No Peak Data Found.");
		G_UNLOCK;
		return;
	}

	if (this->data()->species_ != soap::Species::Human) {
		G_WARN("Motif Location Now Only Support Human");
		G_UNLOCK;
		return;
	}
	
	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Motif Location Settings",
		{ "Database" },
		{soap::InputStyle::ComboBox},
		{ {"JASPAR2024", "Pando"}}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString database_name = (settings[0] == "JASPAR2024") ? FILE_JASPAR2024_HUMAN : FILE_PANDO_MOTIF;

	if (auto item = this->motif_position()) {
		item->__remove_this();
	}
	
	MotifLocateWorker* worker = new MotifLocateWorker(atac_counts->rownames_, database_name, this->data()->species_);
	G_LINK_WORKER_THREAD(MotifLocateWorker, x_motif_location_ready, SingleCellMultiomeItem, s_receive_motif_location);
};

void SingleCellMultiomeItem::s_calculate_gene_activity() {

	auto fragments = this->data()->fragments();
	if (fragments == nullptr) {
		G_WARN("Fragments have not been loaded yet.");
	}

	G_GETLOCK;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Gene Activity Settings",
		{ "Use TSS:yes" },
		{soap::InputStyle::SwitchButton}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	CalculateGeneActivityWorker* worker = new CalculateGeneActivityWorker(
		this->data()->species_, 
		fragments, 
		switch_to_bool(settings[0])
	);
	G_LINK_WORKER_THREAD(CalculateGeneActivityWorker, x_gene_activity_ready, SingleCellMultiomeItem, s_receive_gene_activity);

};

void SingleCellMultiomeItem::s_receive_gene_activity(SparseInt* counts) {

	auto trans = this->get_field_item(DataField::DataType::Trans);

	if (trans != nullptr) {
		trans->__remove_this();
	}

	trans = this->create_field_item(DataField::DataType::Trans, "Gene Activity");

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_TRANS_COUNTS);
	counts->data_type_ = SparseInt::DataType::Counts;

	SUBMODULES(*trans->data(), SparseInt)[title] = std::move(*counts);
	delete counts;

	trans->__check_data();
};

void SingleCellMultiomeItem::s_call_peaks_by_macs() {
	G_GETLOCK;

	auto fragments = this->data()->fragments();
	if (fragments == nullptr) {
		G_WARN("No Fragments Loaded.");
		G_UNLOCK;
		return;
	}

	MacsCallPeakWorker* worker = new MacsCallPeakWorker({ fragments });
	G_LINK_WORKER_THREAD(MacsCallPeakWorker, x_genomic_range_ready, SingleCellMultiomeItem, s_receive_macs_peaks);
}

void SingleCellMultiomeItem::s_receive_qc(QMap<QString, QList<double>> qc) {

	auto&& metadata = this->data()->metadata()->mat_;

	metadata.update(METADATA_TSS_ENRICHMENT, qc[METADATA_TSS_ENRICHMENT]);
	metadata.update(METADATA_TSS_PERCENTILE, qc[METADATA_TSS_PERCENTILE]);
	metadata.update(METADATA_NUCLEOSOME_SIGNAL, qc[METADATA_NUCLEOSOME_SIGNAL]);
	metadata.update(METADATA_BLACKLIST_RATIO, qc[METADATA_BLACKLIST_RATIO]);
	metadata.update(METADATA_FRIP_SCORE, qc[METADATA_FRIP_SCORE]);

	this->data()->double_vectors_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION] = qc[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION];
};

void SingleCellMultiomeItem::s_calculate_quality_metrics() {
	G_GETLOCK;
		
	if (this->data()->atac_counts() == nullptr) {
		G_WARN("No ATAC Data Found.");
		G_UNLOCK;
		return;
	}

	auto fragments = this->data()->fragments();
	if (fragments == nullptr) {
		G_WARN("No Fragments Found.");
		G_UNLOCK;
		return;
	}

	G_LOG("Start Quality Control from Fragments.");

	FragmentsQualityViewWorker* worker = new FragmentsQualityViewWorker(this->data(), fragments);
	G_LINK_WORKER_THREAD(FragmentsQualityViewWorker, x_qc_ready, SingleCellMultiomeItem, s_receive_qc);
};

void SingleCellMultiomeItem::s_filter_by_quality_metrics() {

	const auto& metadata = this->data()->metadata()->mat_;

	if (!metadata.contains(METADATA_TSS_ENRICHMENT, CustomMatrix::DataType::DoubleNumeric)
		|| !metadata.contains(METADATA_NUCLEOSOME_SIGNAL, CustomMatrix::DataType::DoubleNumeric)
		|| !metadata.contains(METADATA_BLACKLIST_RATIO, CustomMatrix::DataType::DoubleNumeric)
		|| !metadata.contains(METADATA_FRIP_SCORE, CustomMatrix::DataType::DoubleNumeric)
		)
	{

		G_WARN("Quality information has been lost.");
		return;
	}

	G_GETLOCK;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Filter Settings",
		{ "Minimum TSS Enrichment Score:3", "Maximum Nucleosome Signal:4",
		"Minimum Frip Score:0.15", "Maximum Black List reads ratio:0.05", "In Place:yes" },
		{soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit,
		soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit, soap::InputStyle::SwitchButton}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	double min_tss = settings[0].toDouble();
	double max_ns = settings[1].toDouble();
	double min_frip = settings[2].toDouble();
	double max_black = settings[3].toDouble();
	bool in_place = switch_to_bool(settings[4]);

	Eigen::ArrayX<bool> filter = _Cs greater_equal(metadata.get_const_double_reference(METADATA_TSS_ENRICHMENT), min_tss);
	filter *= _Cs less_equal(metadata.get_const_double_reference(METADATA_NUCLEOSOME_SIGNAL), max_ns);
	filter *= _Cs less_equal(metadata.get_const_double_reference(METADATA_BLACKLIST_RATIO), max_black);
	filter *= _Cs greater_equal(metadata.get_const_double_reference(METADATA_FRIP_SCORE), min_frip);

	if (filter.count() == 0) {
		G_WARN("No cell remain after filter.");
		G_UNLOCK;
		return;
	}
	else {
		this->col_slice(filter, in_place);
	}
};

void SingleCellMultiomeItem::s_show_quality_matrics() {

	const auto& metadata = this->data()->metadata()->mat_;

	if (!metadata.contains(METADATA_TSS_ENRICHMENT, CustomMatrix::DataType::DoubleNumeric)
		|| !metadata.contains(METADATA_NUCLEOSOME_SIGNAL, CustomMatrix::DataType::DoubleNumeric)
		|| !metadata.contains(METADATA_BLACKLIST_RATIO, CustomMatrix::DataType::DoubleNumeric)
		|| !metadata.contains(METADATA_FRIP_SCORE, CustomMatrix::DataType::DoubleNumeric)
		) 
	{

		G_WARN("Quality information has been lost.");
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, layout] = _Cp prepare_lg(gs);

	auto colors = gs.palette({ METADATA_TSS_ENRICHMENT , METADATA_NUCLEOSOME_SIGNAL ,
		METADATA_BLACKLIST_RATIO , METADATA_FRIP_SCORE });
	const int control_point_number = 24;

	auto trans = metadata.get_double(METADATA_TSS_ENRICHMENT);
	auto axis_rect = _CpPatch new_axis_rect(draw_area);
	layout->addElement(0, 0, axis_rect);
	if (gs.use_boxplot()) {
		_CpPatch single_box_plot(draw_area, axis_rect, trans, colors[0], "", "TSS Score");
	}
	else {
		_CpPatch single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[0], "", "TSS Score");
	}
	_CpPatch clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_NUCLEOSOME_SIGNAL);
	axis_rect = _CpPatch new_axis_rect(draw_area);
	layout->addElement(0, 1, axis_rect);
	if (gs.use_boxplot()) {
		_CpPatch single_box_plot(draw_area, axis_rect, trans, colors[1], "", "Nucleosome Signal");
	}
	else {
		_CpPatch single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[1], "", "Nucleosome Signal");
	}
	_CpPatch clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_BLACKLIST_RATIO);
	axis_rect = _CpPatch new_axis_rect(draw_area);
	layout->addElement(0, 2, axis_rect);
	if (gs.use_boxplot()) {
		_CpPatch single_box_plot(draw_area, axis_rect, trans, colors[2], "", "Blacklist Percentage");
	}
	else {
		_CpPatch single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[2], "", "Blacklist Percentage");
	}
	_CpPatch clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_FRIP_SCORE);
	axis_rect = _CpPatch new_axis_rect(draw_area);
	layout->addElement(0, 3, axis_rect);
	if (gs.use_boxplot()) {
		_CpPatch single_box_plot(draw_area, axis_rect, trans, colors[3], "", "In-Peak Percentage");
	}
	else {
		_CpPatch single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[3], "", "In-Peak Percentage");
	}
	_CpPatch clear_bottom_axis(axis_rect);

	_Cp add_title(draw_area, "Quality Metrics", gs);

	this->draw_suite_->update(draw_area);

};

void SingleCellMultiomeItem::s_show_fragments_length_distribution() {

	if (!this->data()->double_vectors_.contains(VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION) || 
		this->data()->double_vectors_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION].size() != 1000) {
		G_WARN("Fragments Length Distribution Data has lost. Please run Quality Control again.");
		return;
	}

	const auto& distribution = this->data()->double_vectors_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION];

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect] = _Cp prepare_ar(gs);

	double max_length = distribution.size() - 1;
	double minimum_x = 0, maximum_x = max_length + 1, maximum_y = std::ranges::max(distribution) * 1.1;

	draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
	QPen pen;
	pen.setWidth(1);
	pen.setColor(_Cs random_color());

	draw_area->graph()->setPen(pen);
	draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
	draw_area->graph()->setData(_Cs linspaced(distribution.size(), 0, max_length), distribution);

	axis_rect->axis(QCPAxis::atLeft)->setRange(0, maximum_y);
	axis_rect->axis(QCPAxis::atBottom)->setRange(minimum_x, maximum_x);

	_Cp set_simple_axis(axis_rect, "Fragment Length(bp)", "Proportion", gs);
	
	_Cp add_title(draw_area, "Fragments Length Distribution", gs);

	this->draw_suite_->update(draw_area);
};

void SingleCellMultiomeItem::s_sample() {

	G_GETLOCK;

	QStringList factors = this->data()->metadata()->mat_.get_factor_name();

	if (factors.isEmpty()) {
		G_UNLOCK;
		G_NOTICE("No Metadata for sampling.");
		return;
	}
	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Sample Setting",
		{ "Sample by", "Downsample Number" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox},
		{ factors, { "100", "200", "500", "1000", "2000", "5000", "10000", "20000", "ALL" }}
	);
	if (setting.isEmpty()) {
		G_UNLOCK;
		return;
	}
	QString feature = setting[0], downsample = setting[1];
	QStringList metadata = this->data()->metadata()->mat_.get_qstring(feature);
	QStringList levels;
	if (this->data()->metadata()->mat_.data_type_[feature] == CustomMatrix::DataType::QStringFactor) {
		levels = this->data()->metadata()->mat_.string_factors_[feature];
	}
	else {
		levels = _Cs cast<QString>(this->data()->metadata()->mat_.integer_factors_[feature]);
	}
	QVector<int> selected;
	if (downsample != "ALL") {
		int downsample_number = downsample.toInt();
		for (const auto& factor : levels) {
			auto index = _Cs match(metadata, factor);
			if (index.size() > downsample_number) {
				index = _Cs sample(index, downsample_number, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		selected = _Cs seq_n(0, metadata.size());
	}

	if (selected.count() == 0) {
		G_UNLOCK;
		G_NOTICE("No cell is selected.")
			return;
	}

	this->signal_emitter_->x_data_create_soon(this->data()->col_reordered(selected), soap::VariableType::SingleCellMultiome, "Sampled SingleCellMultiome");

	G_UNLOCK;
};

void SingleCellMultiomeItem::__s_rename() {
	G_GETLOCK;

	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Rename SingleCellMultiome",
		{ "New Name:" + this->title_, "Overwrite Metadata:Source?:yes" },
		{ soap::InputStyle::StringLineEdit, soap::InputStyle::SwitchButton}
	);

	if (setting.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString new_name = setting[0];

	if (new_name == this->title_) {
		G_UNLOCK;
		return;
	}

	new_name = this->signal_emitter_->get_unique_name(new_name);

	G_UNLOCK;

	this->__clear_reserve_data();

	this->__change_data_name(new_name);

	this->title_ = new_name;
	this->index_tree_->name_ = new_name;
	this->index_tree_->children_.clear();
	this->signal_emitter_->update_information(new_name, this->data_type_, this->data_, this->is_atomic());
	this->setText(0, new_name);

	this->__check_data();

	if (switch_to_bool(setting[1])) {
		this->data()->metadata()->mat_.update("Source", QStringList(this->data()->rna_counts()->mat_.cols(), this->title_), CustomMatrix::DataType::QStringFactor);
	}

};

void SingleCellMultiomeItem::s_filter_by_features() {

	G_GETLOCK;

	LogicHandler lh(this->data());

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Filter Settings",
		{ "Set Standard", "Filter in place" },
		{soap::InputStyle::LogicLayout, soap::InputStyle::SwitchButton},
		{},
		{},
		{&lh}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto filter = lh.resolve(settings[0]);
	bool in_place = switch_to_bool(settings[1]);

	if (filter.count() == 0) {
		G_WARN("No Cell Meets Requirements.");
		G_UNLOCK;
		return;
	}
	if (filter.count() == filter.size()) {
		G_WARN("No Cell is Excluded.");
		G_UNLOCK;
		return;
	}

	this->col_slice(filter, in_place);

	G_UNLOCK;
}

void SingleCellMultiomeItem::s_column_slice(Eigen::ArrayX<bool> filter, bool in_place) {
	col_slice(filter, in_place);
};

void SingleCellMultiomeItem::s_filter_by_parameters() {
	G_GETLOCK;

	QStringList parameters = CommonDialog::get_response(
		this->signal_emitter_,
		"Set filter parameters",
		{ 
			"Max RNA UMI Count:30000" , 
			"Max Unique Gene:10000" , 
			"Min RNA UMI Count:1000" ,
			"Min Unique Gene:500" , 
			"Max ATAC UMI Count:60000" , 
			"Max ATAC Peak Count:30000" ,
			"Min ATAC UMI Count:3000" , 
			"Min ATAC Peak Count:2000" , 
			"Max Mitochondrial Content(%):10" ,
			"Min Cell Per Gene:5" , 
			"In Place:no" },
		{ 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::SwitchButton
		}
	);
	if (parameters.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto filter = this->get_parameter_filter(parameters.sliced(0, 10));

	if (filter.second.count() == filter.second.size()) {
		G_NOTICE("No cell is excluded.");
	}

	if (filter.second.count() == 0) {
		G_NOTICE("No cell left after filtering.");
		G_UNLOCK;
		return;
	}

	if (filter.first.count() == 0) {
		G_NOTICE("No Gene left after filtering.");
		G_UNLOCK;
		return;
	}

	bool filter_row = (filter.first.count() != filter.first.size());

	bool in_place = switch_to_bool(parameters[10]);
	if (in_place) {
		G_LOG("Filter Object by Parameters...");
	}
	else {
		G_LOG("Filter Object and Create New Object by Parameters...");
	}
	this->slice(filter.first, filter.second, in_place);

}

void SingleCellMultiomeItem::slice(const Eigen::ArrayX<bool>& row_slice, const Eigen::ArrayX<bool>& col_slice, bool in_place) {

	if (in_place) {

		this->__clear_reserve_data();

		this->signal_emitter_->update_information(this->title_, soap::VariableType::SingleCellMultiome, this->data(), true);

		this->data()->col_slice(col_slice);

		this->data()->rna_counts()->row_slice(row_slice);

		this->__check_data();

		G_UNLOCK;
	}
	else {
		SingleCellMultiome* data = this->data()->col_sliced(col_slice);

		data->rna_counts()->row_slice(row_slice);

		G_UNLOCK;

		this->signal_emitter_->x_data_create_soon(data, soap::VariableType::SingleCellMultiome, "Sliced scMultiome");
	}
	G_LOG("Filter finished.");
};

void SingleCellMultiomeItem::col_slice(const Eigen::ArrayX<bool>& filter, bool in_place) {

	if (in_place) {

		this->__clear_reserve_data();

		this->signal_emitter_->update_information(this->title_, soap::VariableType::SingleCellMultiome, this->data(), true);

		this->data()->col_slice(filter);

		this->__check_data();

		G_UNLOCK;
	}
	else {
		SingleCellMultiome* data = this->data()->col_sliced(filter);

		G_UNLOCK;

		this->signal_emitter_->x_data_create_soon(data, soap::VariableType::SingleCellMultiome, "Sliced scMultiome");
	}

	G_LOG("Filter finished.");
};

std::pair<Eigen::ArrayX<bool>, Eigen::ArrayX<bool> > SingleCellMultiomeItem::get_parameter_filter(const QStringList& parameters) {

	const int nrow = this->data()->rna_counts()->mat_.rows();
	const int ncol = this->data()->rna_counts()->mat_.cols();
	Eigen::ArrayX<bool> passed_column = Eigen::ArrayX<bool>::Constant(ncol, true);
	Eigen::ArrayX<bool> passed_row = Eigen::ArrayX<bool>::Constant(nrow, true);

	auto& metadata = this->data()->metadata()->mat_;

	if (!parameters[0].isEmpty() || !parameters[2].isEmpty()) {
		if (metadata.contains(METADATA_RNA_UMI_NUMBER) &&
			metadata.data_type_[METADATA_RNA_UMI_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {
			auto& n_umi = metadata.get_const_integer_reference(METADATA_RNA_UMI_NUMBER);
			if (!parameters[0].isEmpty()) {
				int maximum_count_number = parameters[0].toInt();
				passed_column *= _Cs less_equal(n_umi, maximum_count_number);
			}
			if (!parameters[2].isEmpty()) {
				int minimum_count_number = parameters[2].toInt();
				passed_column *= _Cs greater_equal(n_umi, minimum_count_number);
			}

		}
	}

	if (!parameters[1].isEmpty() || !parameters[3].isEmpty()) {
		if (metadata.contains(METADATA_RNA_UNIQUE_GENE_NUMBER) &&
			metadata.data_type_[METADATA_RNA_UNIQUE_GENE_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {
			auto& n_gene = metadata.get_const_integer_reference(METADATA_RNA_UNIQUE_GENE_NUMBER);

			if (!parameters[1].isEmpty()) {
				int maximum_gene_number = parameters[1].toInt();
				passed_column *= _Cs less_equal(n_gene, maximum_gene_number);
			}
			if (!parameters[3].isEmpty()) {
				int minimum_gene_number = parameters[3].toInt();
				passed_column *= _Cs greater_equal(n_gene, minimum_gene_number);
			}

		}
	}

	if (!parameters[4].isEmpty() || !parameters[6].isEmpty()) {
		if (metadata.contains(METADATA_ATAC_UMI_NUMBER) &&
			metadata.data_type_[METADATA_ATAC_UMI_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {
			auto& n_umi = metadata.get_const_integer_reference(METADATA_ATAC_UMI_NUMBER);
			if (!parameters[4].isEmpty()) {
				int maximum_count_number = parameters[4].toInt();
				passed_column *= _Cs less_equal(n_umi, maximum_count_number);
			}
			if (!parameters[6].isEmpty()) {
				int minimum_count_number = parameters[6].toInt();
				passed_column *= _Cs greater_equal(n_umi, minimum_count_number);
			}

		}
	}

	if (!parameters[5].isEmpty() || !parameters[7].isEmpty()) {
		if (metadata.contains(METADATA_ATAC_UNIQUE_PEAK_NUMBER) &&
			metadata.data_type_[METADATA_ATAC_UNIQUE_PEAK_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {
			auto& n_gene = metadata.get_const_integer_reference(METADATA_ATAC_UNIQUE_PEAK_NUMBER);

			if (!parameters[5].isEmpty()) {
				int maximum_gene_number = parameters[5].toInt();
				passed_column *= _Cs less_equal(n_gene, maximum_gene_number);
			}
			if (!parameters[7].isEmpty()) {
				int minimum_gene_number = parameters[7].toInt();
				passed_column *= _Cs greater_equal(n_gene, minimum_gene_number);
			}

		}
	}

	if (!parameters[8].isEmpty()) {
		if (this->data()->metadata()->mat_.contains(METADATA_RNA_MITOCHONDRIAL_CONTENT) &&
			this->data()->metadata()->mat_.data_type_[METADATA_RNA_MITOCHONDRIAL_CONTENT] == CustomMatrix::DataType::DoubleNumeric) {
			auto& mitochondrial_content = metadata.get_const_double_reference(METADATA_RNA_MITOCHONDRIAL_CONTENT);
			double maximum_mitochondrial_content = parameters[8].toDouble() / 100;
			passed_column *= _Cs less_equal(mitochondrial_content, maximum_mitochondrial_content);
		}
	}

	if (!parameters[9].isEmpty()) {
		passed_row *= _Cs row_count<int, true, false>(this->data()->rna_counts()->mat_, 0.) >= parameters[9].toInt();
	}

	return std::make_pair(passed_row, passed_column);
}

void SingleCellMultiomeItem::s_set_scrublet_threshold() {

	if (!this->data()->double_vectors_.contains(VECTOR_SCRUBLET_SCORES) ||
		!this->data()->double_vectors_.contains(VECTOR_SCRUBLET_SCORES_SIMULATED)) {

		G_WARN("Scrublet is not performed");
	}

	QStringList ret = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Threshold for Scrublet",
		{ "Threshold (0~1)" },
		{ soap::InputStyle::NumericLineEdit}
	);

	if (ret.isEmpty())return;

	double threshold = ret[0].toDouble();
	if (threshold <= 0 || threshold >= 1) {
		G_LOG("Please set a threshold between 0 and 1.");
		return;
	}
	QVector<double> scores = this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES];
	const int n_cell = scores.length();

	if (n_cell != this->data()->rna_counts()->cols()) {
		G_WARN("Invalid Scrublet Data.");
		return;
	}

	QStringList labels;

	for (int i = 0; i < n_cell; ++i) {
		if (scores[i] < threshold) {
			labels << "Singlet";
		}
		else {
			labels << "Doublet";
		}
	}

	Eigen::ArrayXd original_score = _Cs cast<Eigen::ArrayX>(scores);
	Eigen::ArrayXd simulate_score = _Cs cast<Eigen::ArrayX>(this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES_SIMULATED]);

	auto [counts, edges, locs] = _Cs histogram(original_score, 32);
	Eigen::ArrayXd normed_counts = log10(counts.cast<double>() + 1);

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect] = _Cp prepare_ar(gs);
	_Cp bar_plot(draw_area, axis_rect, Qt::blue, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);
	_CpPatch line(draw_area, axis_rect, QVector<double>{threshold, threshold}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
		
	draw_area->plotLayout()->insertRow(0);
	SoapTextElement* title = new SoapTextElement(draw_area, "Observed transcriptomes", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(0, 0, title);

	std::tie(counts, edges, locs) = _Cs histogram(simulate_score, 32);
	normed_counts = log10(counts.cast<double>() + 1);
	axis_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(2, 0, axis_rect);
	_Cp bar_plot(draw_area, axis_rect, Qt::red, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);
	_CpPatch line(draw_area, axis_rect, QVector<double>{threshold, threshold}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);

	draw_area->plotLayout()->insertRow(2);
	SoapTextElement* title2 = new SoapTextElement(draw_area, "Simulated doublets", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(2, 0, title2);

	this->draw_suite_->update(draw_area);

	this->data()->metadata()->mat_.update(METADATA_SCRUBLET_LABELS, labels, CustomMatrix::DataType::QStringFactor);
	G_LOG("Threshold : " + QString::number(threshold) + " has been set for Scrublet.");

	double original_rate = _Cs greater_than(original_score, threshold).count() / (double)original_score.size();
	double simulate_rate = _Cs greater_than(simulate_score, threshold).count() / (double)simulate_score.size();

	G_LOG("Detected " + QString::number(original_rate * 100) + " % doublets in original data and " + QString::number(simulate_rate * 100) + " % doublets in simulated doublets.");
};

void SingleCellMultiomeItem::s_receive_scrublet(Eigen::ArrayXd original_score, Eigen::ArrayXd simulate_score) {
	
	this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES] = _Cs cast<QVector>(original_score);
	this->data()->metadata()->mat_.update(METADATA_SCRUBLET_SCORES, this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES]);
	this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES_SIMULATED] = _Cs cast<QVector>(simulate_score);

	auto threshold = _Cs threshold_minimum(simulate_score);

	if (!threshold.first) {
		G_LOG("Failed to automatically find a threshold. Please determine the threshold yourself by [Set Threshold]");
	}
	else {
		G_LOG("The threshold has been automatically set to " + QString::number(threshold.second) + ".");
		double original_rate = (original_score > threshold.second).count() / (double)original_score.size();
		double simulate_rate = (simulate_score > threshold.second).count() / (double)simulate_score.size();

		G_LOG("Detected " + QString::number(original_rate * 100) + " % doublets in original data and " + QString::number(simulate_rate * 100) + " % doublets in simulated doublets.");

		double threshold_value = threshold.second;
		const int n_cell = original_score.size();
		QStringList labels;
		for (int i = 0; i < n_cell; ++i) {
			if (original_score[i] < threshold_value) {
				labels << "Singlet";
			}
			else {
				labels << "Doublet";
			}
		}
		this->data()->metadata()->mat_.update(METADATA_SCRUBLET_LABELS, labels, CustomMatrix::DataType::QStringFactor);
	}

	auto [counts, edges, locs] = _Cs histogram(original_score, 32);
	Eigen::ArrayXd normed_counts = log10(counts.cast<double>() + 1);

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect] = _Cp prepare_ar(gs);
	_Cp bar_plot(draw_area, axis_rect, Qt::blue, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);

	if (threshold.first) {
		double threshold_value = threshold.second;
		_CpPatch line(draw_area, axis_rect, QVector<double>{threshold_value, threshold_value}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
	}

	draw_area->plotLayout()->insertRow(0);
	SoapTextElement* title = new SoapTextElement(draw_area, "Observed transcriptomes", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(0, 0, title);

	std::tie(counts, edges, locs) = _Cs histogram(simulate_score, 32);
	normed_counts = log10(counts.cast<double>() + 1);
	axis_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(2, 0, axis_rect);
	_Cp bar_plot(draw_area, axis_rect, Qt::red, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);
	if (threshold.first) {
		double threshold_value = threshold.second;
		_CpPatch line(draw_area, axis_rect, QVector<double>{threshold_value, threshold_value}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
	}

	draw_area->plotLayout()->insertRow(2);
	SoapTextElement* title2 = new SoapTextElement(draw_area, "Simulated doublets", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(2, 0, title2);

	this->draw_suite_->update(draw_area);

	G_LOG("Detect doublets by Scrublet finished.");
}

void SingleCellMultiomeItem::s_rna_scrublet() {
	G_GETLOCK;

	auto rna_counts = this->data()->rna_counts();
	if (rna_counts == nullptr) {
		G_WARN("Counts data missed.");
		G_UNLOCK;
		return;
	}
	G_LOG("Detect doublets by Scrublet...");
	ScrubletWorker* worker = new ScrubletWorker(rna_counts->mat_);
	G_LINK_WORKER_THREAD(ScrubletWorker, x_scrublet_ready, SingleCellMultiomeItem, s_receive_scrublet)
};

void SingleCellMultiomeItem::__s_delete_this() {
	G_GETLOCK;

	G_UNLOCK;

	if (!YesOrNoDialog::get_response("Delete SingleCellMultiome", "This data will be deleted.")) {
		return;
	}

	this->__remove_this();
};

void SingleCellMultiomeItem::s_set_random_state() {
	G_GETLOCK;
	QStringList seed = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Random Seed",
		{ "New Random Seed:" + QString::number(this->data()->random_state_) },
		{ soap::InputStyle::IntegerLineEdit}
	);
	if (seed.isEmpty()) {
		G_UNLOCK;
		return;
	}
	this->data()->random_state_ = seed[0].toInt();
	G_UNLOCK;
};

void SingleCellMultiomeItem::s_set_species() {
	G_GETLOCK;
	QStringList species = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Species",
		{ "New Species" },
		{ soap::InputStyle::ComboBox},
		{ { "Human", "Mouse", "Undefined" }}
	);
	if (species.isEmpty()) {
		G_UNLOCK;
		return;
	}
	if (species[0] == "Human") {
		this->data()->species_ = soap::Species::Human;
	}
	else if (species[0] == "Mouse") {
		this->data()->species_ = soap::Species::Mouse;
	}
	else {
		this->data()->species_ = soap::Species::Undefined;
	}

	G_UNLOCK;
};

void SingleCellMultiomeItem::s_receive_gsea(GSEA gsea) {

	QString title = "GSEA " + gsea.database_ + " | " + gsea.comparison_[1] + " | " + gsea.comparison_[2];
	title = this->signal_emitter_->get_unique_name(title);

	DATA_SUBMODULES(GSEA)[title] = gsea;

	GSEAItem* item = new GSEAItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(GSEA)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);
	G_LOG("GSEA between " + gsea.comparison_[1] + " and " + gsea.comparison_[2] + " in " + gsea.comparison_[0] + " finished");
};

void SingleCellMultiomeItem::s_gsea_from_input() {};

void SingleCellMultiomeItem::s_gsea_in_database() {

	auto rna_normalized = this->data()->rna_normalized();
	if (rna_normalized == nullptr) {
		G_WARN("RNA data has not been normalized.");
		return;
	}
	G_GETLOCK;

	auto factor_map = this->data()->metadata()->mat_.get_factor_information(false);
	if (factor_map.isEmpty()) {
		G_LOG("No suitable feature for GSEA");
		G_UNLOCK;
		return;
	}

	auto species = this->data()->species_;

	if (species != soap::Species::Human && species != soap::Species::Mouse) {
		G_WARN("Only Human and Mouse data GSEA is supported.");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"GSEA settings",
		{ "Choose Comparison:2", "Database", "Min Overlap:10",
		"Min Pathway Size:10", "Max Pathway Size:500", "Random State:" + QString::number(this->data()->random_state_),
		"Permutation Type", "Downsample Size:1000" },
		{ soap::InputStyle::CompareLayout, soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit},
		{ { "Curated", "Ontology" }, { "Phenotype", "Geneset" }},
		{ factor_map }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto comparison = compare_layouts_to_list(settings[0]);
	QString factor_name = comparison[0], comparison1 = comparison[1], comparison2 = comparison[2];
	if (comparison1 == comparison2) {
		G_LOG("Invalid comparison!");
		G_UNLOCK;
		return;
	}
	QString database = settings[1];
	int minimum_overlap = settings[2].toInt();
	int minimum_size = settings[3].toInt();
	int maximum_size = settings[4].toInt();
	unsigned int random_state = settings[5].toUInt();
	QString permutation_type = settings[6];

	int downsample = settings[7].toInt();
	if (downsample <= 10) {
		G_WARN("Illegal downsample number!");
		G_UNLOCK;
		return;
	}

	QVector<int> selected;
	QStringList metadata = this->data()->metadata()->mat_.get_qstring(factor_name);

	auto index = _Cs match(metadata, comparison1);
	if (index.size() > downsample) {
		index = _Cs sample(index, downsample, random_state);
	}
	selected << index;

	if (comparison2 == "REST") {
		index = _Cs which(!_Cs equal(metadata, comparison1));
	}
	else {
		index = _Cs match(metadata, comparison2);
	}
	if (index.size() > downsample) {
		index = _Cs sample(index, downsample, random_state);
	}
	selected << index;

	metadata = _Cs reordered(metadata, selected);
	G_LOG("GSEAing in " + factor_name + "...");
	SparseDouble tmp = rna_normalized->col_reordered(selected);

	GseaWorker* worker = new GseaWorker(
		tmp, 
		metadata, 
		comparison, 
		species,
		database, 
		minimum_overlap, 
		minimum_size, 
		maximum_size, 
		random_state, 
		permutation_type,
		1000
	);
	G_LINK_WORKER_THREAD(GseaWorker, x_gsea_ready, SingleCellMultiomeItem, s_receive_gsea)
}

void SingleCellMultiomeItem::s_infercnv() {

	auto counts = this->data()->rna_counts();
	if (counts == nullptr) {
		G_WARN("RNA data not found.");
		return;
	}
	G_GETLOCK;
	auto factor_map = this->data()->metadata()->mat_.get_factor_information();
	if (factor_map.isEmpty()) {
		G_LOG("CNV analysis requires at least one factor for visualisation.");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"InferCNV Settings",
		{ "Feature:Reference", "Downsample Number",},
		{ soap::InputStyle::FactorChoice, soap::InputStyle::ComboBox},
		{ { "200", "300", "500", "1000", "ALL" } },
		{ factor_map }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	};

	auto [feature, reference] = factor_choice_to_pair(settings[0]);

	QStringList metadata = this->data()->metadata()->mat_.get_qstring(feature);
	QStringList levels = _Cs unique(metadata);

	QVector<int> selected;
	if (settings[1] != "ALL") {
		int downsample = settings[1].toInt();
		for (const auto& factor : levels) {
			auto index = _Cs match(metadata, factor);
			if (index.size() > downsample) {
				index = _Cs sample(index, downsample, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		for (const auto& factor : levels) {
			auto index = _Cs match(metadata, factor);
			selected << index;
		}
	}

	if (reference.isEmpty()) {
		reference = _Cs unique(metadata);
	}
	else {
		reference = _Cs unique(reference);
	}

	G_LOG("InferCnv start...");
	InferCnvWorker* worker = new InferCnvWorker(
		counts->col_reordered(selected),
		_Cs reordered(metadata, selected),
		reference,
		this->data()->species_
	);

	G_LINK_WORKER_THREAD(InferCnvWorker, x_cnv_ready, SingleCellMultiomeItem, s_receive_infercnv);
};

void SingleCellMultiomeItem::s_receive_infercnv(CNV* cnv) {

	QString title = this->signal_emitter_->get_unique_name("InferCNV");
	DATA_SUBMODULES(CNV)[title] = std::move(*cnv);
	delete cnv;

	CNVItem* item = new CNVItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(CNV)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("InferCNV finished");
};

void SingleCellMultiomeItem::s_scicnv() {

	auto rna_normalized = this->data()->rna_normalized();
	if (rna_normalized == nullptr) {
		G_WARN("RNA data has not been normalized.");
		return;
	}

	G_GETLOCK;
	auto factor_map = this->data()->metadata()->mat_.get_factor_information();
	if (factor_map.isEmpty()) {
		G_LOG("CNV analysis requires at least one factor for visualisation.");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Scicnv Settings",
		{ "Feature:Reference", "Downsample Number", "Filter Threshold:0.4", "Sharpness(0.6~1.4):1.0" },
		{ soap::InputStyle::FactorChoice, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit},
		{ { "200", "300", "500", "1000", "ALL" }},
		{ factor_map}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	};

	auto [feature, reference] = factor_choice_to_pair(settings[0]);

	double threshold = settings[2].toDouble();
	double sharpness = settings[3].toDouble();
	if (sharpness < 0.6 || sharpness > 1.4) {
		G_LOG("Sharpness must be between 0.6 and 1.4");
		G_UNLOCK;
		return;
	}

	QStringList metadata = this->data()->metadata()->mat_.get_qstring(feature);
	QStringList levels;
	if (this->data()->metadata()->mat_.data_type_[feature] == CustomMatrix::DataType::QStringFactor) {
		levels = this->data()->metadata()->mat_.string_factors_[feature];
	}
	else {
		levels = _Cs cast<QString>(this->data()->metadata()->mat_.integer_factors_[feature]);
	}
	QVector<int> selected;
	if (settings[1] != "ALL") {
		int downsample = settings[1].toInt();
		for (const auto& factor : levels) {
			auto index = _Cs match(metadata, factor);
			if (index.size() > downsample) {
				index = _Cs sample(index, downsample, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		for (const auto& factor : levels) {
			auto index = _Cs match(metadata, factor);
			selected << index;
		}
	}
	G_LOG("SciCnv start...");
	ScicnvWorker* worker = new ScicnvWorker(
		rna_normalized->col_reordered(selected),
		_Cs reordered(metadata, selected),
		reference,
		this->data()->species_,
		threshold,
		sharpness
	);

	G_LINK_WORKER_THREAD(ScicnvWorker, x_cnv_ready, SingleCellMultiomeItem, s_receive_scicnv);
};

void SingleCellMultiomeItem::s_receive_scicnv(CNV* cnv) {

	QString title = this->signal_emitter_->get_unique_name("SciCnv");
	DATA_SUBMODULES(CNV)[title] = std::move(*cnv);
	delete cnv;

	CNVItem* item = new CNVItem(
		title, 
		this->index_tree_,
		&DATA_SUBMODULES(CNV)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("SciCNV finished");
};

void SingleCellMultiomeItem::s_edit_metadata() {

	G_GETLOCK;

	LogicHandler lh(this->data());

	auto&& metadata = this->data()->metadata()->mat_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Filter Settings",
		{ "Select Cell", "Feature to edit", "New Value"},
		{soap::InputStyle::LogicLayout, soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit},
		{metadata.colnames_},
		{},
		{&lh}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto filter = lh.resolve(settings[0]);
	QString feature_name = settings[1];
	QString new_val = settings[2];

	if (filter.count() == 0) {
		G_WARN("No Cell Meets Requirements.");
		G_UNLOCK;
		return;
	}

	metadata.edit(feature_name, filter, new_val);
	G_UNLOCK;
};

void SingleCellMultiomeItem::s_combine_existed_metadata() {
	G_GETLOCK;

	QStringList metadata_names = this->data()->metadata()->mat_.colnames_;

	if (metadata_names.isEmpty()) {
		G_WARN("Metadata is Empty.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Metadata Combine Settings",
		{ "Choose Metadata", "Combine Type", "New Metadata Name", "New Metadata Type", "Join Character: "},
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::ComboBox,
		soap::InputStyle::StringLineEdit, soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ metadata_names,
		{ "Text Concatenate", "Numeric Add", "Numeric Minus", "Numeric Multiply", "Numeric Divide" },
		{ "String Factor", "String", "Numeric", "Integer", "Integer Factor" }}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QStringList components = simple_choice_to_list(settings[0]);

	if (components.size() < 2) {
		G_UNLOCK;
		return;
	}

	QString combine_type = settings[1];
	QString target_data_type = settings[3];
	QString new_metadata_name = settings[2];
	QString joiner = settings[4];

	auto& metadata = this->data()->metadata()->mat_;

	if (metadata.contains(new_metadata_name)) {
		if (!YesOrNoDialog::get_response("Warning: Name Duplicated!", "Overwrite existing metadata?")) {
			G_UNLOCK;
			return;
		}
	}

	if (combine_type == "Text Concatenate") {
		QList<QStringList> choosed;
		for (auto&& name : components) {
			choosed << metadata.get_qstring(name);
		}

		const int nrow = metadata.rows();
		const int n_choosed = choosed.size();

		QStringList res = choosed[0];
		for (int j = 1; j < n_choosed; ++j) {
			for (int i = 0; i < nrow; ++i) {
				res[i] += (joiner + choosed[j][i]);
			}
		}

		if (target_data_type == "Integer") {
			this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<int>(res), CustomMatrix::DataType::IntegerNumeric);
		}
		else if (target_data_type == "Numeric") {
			this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<double>(res), CustomMatrix::DataType::DoubleNumeric);
		}
		else if (target_data_type == "Integer Factor") {
			this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<int>(res), CustomMatrix::DataType::IntegerFactor);
		}
		else if (target_data_type == "String Factor") {
			this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::QStringFactor);
		}
		else if (target_data_type == "String") {
			this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::QString);
		}
	}
	else {
		auto data_types = _Cs sapply(components, [&metadata](auto&& name) {return metadata.data_type_[name]; });

		const int nrow = metadata.rows();
		const int n_choosed = data_types.size();

		if (data_types.contains(CustomMatrix::DataType::DoubleNumeric) || target_data_type == "Numeric") {
			QList<QVector<double>> choosed;
			for (auto&& name : components) {
				choosed << metadata.get_double(name);
			}

			QVector<double> res = choosed[0];

			if (combine_type == "Numeric Add") {
				for (int j = 1; j < n_choosed; ++j) {
					for (int i = 0; i < nrow; ++i) {
						res[i] += choosed[j][i];
					}
				}
			}
			else if(combine_type == "Numeric Minus") {
				if (n_choosed != 2) {
					G_WARN("Numeric Minus only support two metadata operation.");
					G_UNLOCK;
					return;
				}
				for (int i = 0; i < nrow; ++i) {
					res[i] -= choosed[1][i];
				}
			}
			else if (combine_type == "Numeric Multiply") {
				for (int j = 1; j < n_choosed; ++j) {
					for (int i = 0; i < nrow; ++i) {
						res[i] *= choosed[j][i];
					}
				}
			}
			else if (combine_type == "Numeric Divide") {
				if (n_choosed != 2) {
					G_WARN("Numeric Divide only support two metadata operation.");
					G_UNLOCK;
					return;
				}
				for (int i = 0; i < nrow; ++i) {
					if (choosed[1][i] == 0.0) {
						G_WARN("Meet 0.0 in divisor, process terminate");
						G_UNLOCK;
						return;
					}
					res[i] /= choosed[1][i];
				}
			}

			if (target_data_type == "Integer") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<int>(res), CustomMatrix::DataType::IntegerNumeric);
			}
			else if (target_data_type == "Numeric") {
				this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::DoubleNumeric);
			}
			else if (target_data_type == "Integer Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<int>(res), CustomMatrix::DataType::IntegerFactor);
			}
			else if (target_data_type == "String Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<QString>(res), CustomMatrix::DataType::QStringFactor);
			}
			else if (target_data_type == "String") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<QString>(res), CustomMatrix::DataType::QString);
			}
		}
		else {
			QList<QVector<int>> choosed;
			for (auto&& name : components) {
				choosed << metadata.get_integer(name);
			}

			QVector<int> res = choosed[0];

			if (combine_type == "Numeric Add") {
				for (int j = 1; j < n_choosed; ++j) {
					for (int i = 0; i < nrow; ++i) {
						res[i] += choosed[j][i];
					}
				}
			}
			else if (combine_type == "Numeric Minus") {
				if (n_choosed != 2) {
					G_WARN("Numeric Minus only support two metadata operation.");
					G_UNLOCK;
					return;
				}
				for (int i = 0; i < nrow; ++i) {
					res[i] -= choosed[1][i];
				}
			}
			else if (combine_type == "Numeric Multiply") {
				for (int j = 1; j < n_choosed; ++j) {
					for (int i = 0; i < nrow; ++i) {
						res[i] *= choosed[j][i];
					}
				}
			}
			else if (combine_type == "Numeric Divide") {
				if (n_choosed != 2) {
					G_WARN("Numeric Divide only support two metadata operation.");
					G_UNLOCK;
					return;
				}
				for (int i = 0; i < nrow; ++i) {
					if (choosed[1][i] == 0) {
						G_WARN("Meet 0 in divisor, process terminate");
						G_UNLOCK;
						return;
					}
					res[i] /= choosed[1][i];
				}
			}

			if (target_data_type == "Integer") {
				this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::IntegerNumeric);
			}
			else if (target_data_type == "Numeric") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<double>(res), CustomMatrix::DataType::DoubleNumeric);
			}
			else if (target_data_type == "Integer Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::IntegerFactor);
			}
			else if (target_data_type == "String Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<QString>(res), CustomMatrix::DataType::QStringFactor);
			}
			else if (target_data_type == "String") {
				this->data()->metadata()->mat_.update(new_metadata_name, _Cs cast<QString>(res), CustomMatrix::DataType::QString);
			}
		}
	}	

	this->signal_emitter_->x_update_interface();
	G_UNLOCK;
};

void SingleCellMultiomeItem::s_add_new_metadata() {
	G_GETLOCK;

	QStringList metadata_information = CommonDialog::get_response(
		this->signal_emitter_,
		"New Metadata Setting",
		{ "Name", "Type", "Initial Value", "Initiate from existing metadata?:no", "Metadata from" },
		{ soap::InputStyle::StringLineEdit, soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit
		, soap::InputStyle::SwitchButton, soap::InputStyle::ComboBox},
		{ { "String Factor", "String", "Numeric", "Integer", "Integer Factor" }, this->data()->metadata()->mat_.colnames_ }
	);

	if (metadata_information.isEmpty()) {
		G_UNLOCK;
		return;
	}
	if (this->data()->metadata()->mat_.contains(metadata_information[0])) {
		if (!YesOrNoDialog::get_response("Warning : Name Duplicated!", "Overwrite existing metadata?")) {
			G_UNLOCK;
			return;
		}
	}
	bool existed = switch_to_bool(metadata_information[3]);
	if (existed) {
		QString existed_name = metadata_information[4];
		CustomMatrix::DataType datatype = this->data()->metadata()->mat_.data_type_[existed_name];
		if (datatype == CustomMatrix::DataType::DoubleNumeric) {
			this->data()->metadata()->mat_.update(metadata_information[0], this->data()->metadata()->mat_.get_double(existed_name), CustomMatrix::DataType::DoubleNumeric);
		}
		else if (datatype == CustomMatrix::DataType::IntegerNumeric) {
			this->data()->metadata()->mat_.update(metadata_information[0], this->data()->metadata()->mat_.get_integer(existed_name), CustomMatrix::DataType::IntegerNumeric);
		}
		else if (datatype == CustomMatrix::DataType::QStringFactor) {
			this->data()->metadata()->mat_.update(metadata_information[0], this->data()->metadata()->mat_.get_qstring(existed_name), CustomMatrix::DataType::QStringFactor);
		}
		else if (datatype == CustomMatrix::DataType::QString) {
			this->data()->metadata()->mat_.update(metadata_information[0], this->data()->metadata()->mat_.get_qstring(existed_name), CustomMatrix::DataType::QString);
		}
		else if (datatype == CustomMatrix::DataType::IntegerFactor) {
			this->data()->metadata()->mat_.update(metadata_information[0], this->data()->metadata()->mat_.get_integer(existed_name), CustomMatrix::DataType::IntegerFactor);
		}
	}
	else {
		QString datatype = metadata_information[1];
		int size = this->data()->metadata()->mat_.rownames_.size();
		if (datatype == "Integer") {
			this->data()->metadata()->mat_.update(metadata_information[0], QVector<int>(size, metadata_information[2].toInt()), CustomMatrix::DataType::IntegerNumeric);
		}
		else if (datatype == "Numeric") {
			this->data()->metadata()->mat_.update(metadata_information[0], QVector<double>(size, metadata_information[2].toDouble()), CustomMatrix::DataType::DoubleNumeric);
		}
		else if (datatype == "Integer Factor") {
			this->data()->metadata()->mat_.update(metadata_information[0], QVector<int>(size, metadata_information[2].toInt()), CustomMatrix::DataType::IntegerFactor);
		}
		else if (datatype == "String Factor") {
			this->data()->metadata()->mat_.update(metadata_information[0], QStringList(size, metadata_information[2]), CustomMatrix::DataType::QStringFactor);
		}
		else if (datatype == "String") {
			this->data()->metadata()->mat_.update(metadata_information[0], QStringList(size, metadata_information[2]), CustomMatrix::DataType::QString);
		}
	}
	this->signal_emitter_->x_update_interface();
	G_UNLOCK;
};

void SingleCellMultiomeItem::s_receive_integrated_data(SingleCellMultiome* data, QList<const SingleCellMultiome* > items) {
	this->signal_emitter_->unlock(this->signal_emitter_->search(_Cs sapply(items, [](auto* data) {return (void*)data; })));

	this->signal_emitter_->x_data_create_soon(data, soap::VariableType::SingleCellMultiome, "Integrated scMultiome");
};

void SingleCellMultiomeItem::s_integrate() {

	QMap<QString, SingleCellMultiome*> available_data;
	for (const auto& [variable_name, data_info] : this->signal_emitter_->variable_information_) {
		if (data_info.first == soap::VariableType::SingleCellMultiome) {
			available_data[variable_name] = static_cast<SingleCellMultiome*>(data_info.second);
		}
	}
	if (available_data.size() < 2) {
		G_LOG("No enough single cell multiome data found.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Integrate SingleCellMultiome",
		{ "Data", "Style" },
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::ComboBox},
		{ available_data.keys(), { "Distinguish cells from different data", "Merge the same cell name" } }
	);

	if (settings.isEmpty())return;

	std::map<QString, SingleCellMultiome*> choosed_data;

	auto choosed_name = simple_choice_to_list(settings[0]);

	if (choosed_name.size() < 2)return;

	if (!_Cs is_unique(choosed_name)) {
		G_WARN("Can not integrate the same data!");
		return;
	}

	std::ranges::for_each(choosed_name, [&choosed_data, &available_data](const QString& name) {choosed_data[name] = available_data[name]; });

	QList<const SingleCellMultiome*> all_data;
	int velocyto_base_loaded = 0;

	QStringList fragments_missed_data;
	for (const auto& [data_name, ptr] : choosed_data) {

		auto fragments = ptr->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments of " + data_name + " has not been loaded.");
			return;
		}

		if (ptr->velocyto_base() != nullptr) {
			++velocyto_base_loaded;
		}
		all_data << ptr;
	}

	if (velocyto_base_loaded > 0 && velocyto_base_loaded < all_data.size()) {
		auto res = YesOrNoDialog::get_response(
			"Velocyto settings",
			"Some data did not load velocyto matrix. Integrated data cannot generate velocyto matrix anymore due to barcode collision. Still continue?"
		);

		if (!res) {
			return;
		}
	}

	auto species = _Cs unique(_Cs sapply(all_data, [](const SingleCellMultiome* data) {return data->species_; }));

	if (species.size() != 1) {
		G_WARN("Unmatched species!");
		return;
	}

	if (!this->signal_emitter_->try_lock(_Cs sapply(all_data, [this](const SingleCellMultiome* data) {return this->signal_emitter_->search((void*)data); }))) {
		G_WARN("Please waiting for the computation in progress.");
		return;
	}

	bool distinguish = (settings[1] == "Distinguish cells from different data");
	if (!distinguish) {
		G_NOTICE("Metadata will not be inherited to avoid collision.");
	}

	G_LOG("Start integrating data...");
	IntegrateWorker* worker = new IntegrateWorker(all_data, distinguish, species[0]);
	G_LINK_WORKER_THREAD(IntegrateWorker, x_scmultiome_ready, SingleCellMultiomeItem, s_receive_integrated_data);
};

void SingleCellMultiomeItem::s_scent() {

	if (this->data()->species_ != soap::Species::Human) {
		G_WARN("SCENT now only support human");
		return;
	}

	auto rna_counts = this->data()->rna_counts();
	if (rna_counts == nullptr) {
		G_WARN("No RNA Data Found.");
		return;
	}

	G_GETLOCK;

	ScentWorker* worker = new ScentWorker(rna_counts);
	G_LINK_WORKER_THREAD(ScentWorker, x_sr_ready, SingleCellMultiomeItem, s_receive_scent);
};

void SingleCellMultiomeItem::s_receive_scent(Eigen::ArrayXd data) {

	G_LOG("Signalling Entropy Calculation Finished.");

	this->data()->metadata()->mat_.update("Signalling Entropy", _Cs cast<QVector>(data));

};

void SingleCellMultiomeItem::s_pando() {

	if (this->data()->rna_counts() == nullptr) {
		G_WARN("Rna Counts Data is Missing.");
		return;
	}

	if (this->data()->rna_normalized() == nullptr) {
		G_WARN("Rna Normalized Data is Missing.");
		return;
	}

	if (this->data()->atac_counts() == nullptr) {
		G_WARN("Peak Counts Data is Missing.");
		return;
	}

	if (this->data()->atac_normalized() == nullptr) {
		G_WARN("Peak Normalized Data is Missing.");
		return;
	}

	if (this->data()->species_ != soap::Species::Human) {
		G_WARN("Now only human data is supported.");
		return;
	}

	G_GETLOCK;

	PandoWorker* worker{ nullptr };

	auto factors = this->data()->metadata()->mat_.get_factor_information(false);
	if (factors.isEmpty()) {

		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Pando Settings",
			{ "Peak Correlation Threshold:0.05", "Transcriptional Factor Threshold:0.05",
			"Feature selected:6000", "Use Trancriptional Start Site:yes",
			"P Adjust Method", "Database"},
			{soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit,
			soap::InputStyle::IntegerLineEdit, soap::InputStyle::SwitchButton, 
			soap::InputStyle::ComboBox, soap::InputStyle::ComboBox},
			{ { "fdr", "Bonferroni" }, {"PANDO", "JASPAR2024"}}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		double peak_cor = settings[0].toDouble();
		double tf_cor = settings[1].toDouble();

		int n_feature = settings[2].toInt();
		if (n_feature <= 0) {
			G_UNLOCK;
			return;
		}

		bool use_tss = switch_to_bool(settings[3]);
		QString p_adjust_method = settings[4];

		QString database_name = (settings[5] == "JASPAR2024") ? FILE_JASPAR2024_HUMAN : FILE_PANDO_MOTIF;

		worker = new PandoWorker(
			this->data(),
			database_name,
			false,
			use_tss,
			Eigen::ArrayX<bool>(),
			n_feature,
			peak_cor,
			tf_cor,
			p_adjust_method
		);
	}
	else {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Pando Settings",
			{  "Peak Correlation Threshold:0.05", "Transcriptional Factor Threshold:0.05",
			"Select cell:true", "Cell Group", "Feature selected:6000", "Use Trancriptional Start Site:yes",
			"P Adjust Method", "Database"},
			{
				soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit,
				soap::InputStyle::SwitchButton, soap::InputStyle::FactorChoice,
				soap::InputStyle::IntegerLineEdit, soap::InputStyle::SwitchButton, 
				soap::InputStyle::ComboBox, soap::InputStyle::ComboBox
				},
			{  { "fdr", "Bonferroni" }, {"Pando", "JASPAR2024"}},
			{factors}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		double peak_cor = settings[0].toDouble();
		double tf_cor = settings[1].toDouble();

		bool filter_cell = switch_to_bool(settings[2]);
		Eigen::ArrayX<bool> cell_filter;
		if (filter_cell) {

			auto [factor, levels] = factor_choice_to_pair(settings[3]);
			if (levels.isEmpty()) {
				G_UNLOCK;
				return;
			}
			cell_filter = _Cs in(this->data()->metadata()->mat_.get_qstring(factor), levels);
		}

		int n_feature = settings[4].toInt();
		if (n_feature <= 0) {
			G_UNLOCK;
			return;
		}

		bool use_tss = switch_to_bool(settings[5]);
		QString p_adjust_method = settings[6];
		QString database_name = (settings[7] == "JASPAR2024") ? FILE_JASPAR2024_HUMAN : FILE_PANDO_MOTIF;

		worker = new PandoWorker(
			this->data(),
			database_name,
			filter_cell,
			use_tss,
			cell_filter,
			n_feature,
			peak_cor,
			tf_cor,
			p_adjust_method
		);
	}

	G_LINK_WORKER_THREAD(PandoWorker, x_pando_ready, SingleCellMultiomeItem, s_receive_pando);
};

void SingleCellMultiomeItem::s_receive_pando(Pando pando) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_PANDO);

	DATA_SUBMODULES(Pando)[title] = std::move(pando);

	PandoItem* item = new PandoItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(Pando)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Pando finished.");
};


void SingleCellMultiomeItem::s_monocle3() {

	G_GETLOCK;

	auto embedding_names = this->data()->embedding_names();

	if (embedding_names.isEmpty()) {
		G_WARN("No Embedding.");
		G_UNLOCK;
		return;
	}

	auto factor_info = this->data()->metadata()->mat_.get_factor_information(false);

	if (factor_info.isEmpty()) {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Monocle3 Settings",
			{ "Embedding" },
			{soap::InputStyle::ComboBox},
			{embedding_names}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		QString embedding_name = settings[0];
		int n_cell = this->data()->rna_counts()->cols();

		Monocle3Worker* worker = new Monocle3Worker(
			*this->data()->embedding(embedding_name),
			Eigen::ArrayX<bool>::Constant(n_cell, true)
		);

		G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellMultiomeItem, s_receive_monocle3);
	}
	else {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Monocle3 Settings",
			{ "Embedding" , "use part", "choose cell group" },
			{soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton, soap::InputStyle::FactorChoice},
			{embedding_names},
			{factor_info}
		);

		if (settings.isEmpty()) {
			G_UNLOCK;
			return;
		}

		QString embedding_name = settings[0];
		bool use_part = switch_to_bool(settings[1]);
		if (use_part) {

			auto [factor_name, levels] = factor_choice_to_pair(settings[2]);

			if (levels.isEmpty()) {
				G_UNLOCK;
				return;
			}

			auto factor = this->data()->metadata()->mat_.get_qstring(factor_name);

			Monocle3Worker* worker = new Monocle3Worker(
				*this->data()->embedding(embedding_name),
				_Cs in(factor, levels)
			);

			G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellMultiomeItem, s_receive_monocle3);
		}
		else {
			int n_cell = this->data()->rna_counts()->cols();

			Monocle3Worker* worker = new Monocle3Worker(
				*this->data()->embedding(embedding_name),
				Eigen::ArrayX<bool>::Constant(n_cell, true)
			);

			G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellMultiomeItem, s_receive_monocle3);
		}
	}

};

void SingleCellMultiomeItem::s_receive_monocle3(Monocle3* monocle3) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_MONOCLE3);

	DATA_SUBMODULES(Monocle3)[title] = *monocle3;
	delete monocle3;

	Monocle3Item* item = new Monocle3Item(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(Monocle3)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Monocle3 computation finished.");
};

void SingleCellMultiomeItem::s_velocyto() {

	if (this->data()->data_type_ == SingleCellMultiome::DataType::Integrated) {
		G_WARN("Integrated object can not calculate velocyto due to barcode collision.");
		return;
	}

	G_GETLOCK;

	QStringList bam_file_path = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Bam File",
		{ "Bam (*.bam):Bam (*.bam)" },
		{ soap::InputStyle::MultiFile},
		{ QStringList()}
	);

	if (bam_file_path.isEmpty()) {
		G_UNLOCK;
		return;
	}

	bam_file_path = multiple_file_to_list(bam_file_path[0]);

	if (bam_file_path.isEmpty()) {
		G_UNLOCK;
		return;
	}
	
	if (auto item = this->velocyto_base()) {
		item->__remove_this();
	}

	VelocytoWorker* worker = new VelocytoWorker(this->data(), bam_file_path);
	G_LINK_WORKER_THREAD(VelocytoWorker, x_velocyto_ready, SingleCellMultiomeItem, s_receive_velocyto);
};

void SingleCellMultiomeItem::s_receive_velocyto(VelocytoBase* velocyto_base) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_VELOCYTO_BASE);

	DATA_SUBMODULES(VelocytoBase)[title] = std::move(*velocyto_base);
	delete velocyto_base;

	VelocytoBaseItem* item = new VelocytoBaseItem(
		title,
		this->index_tree_, 
		&DATA_SUBMODULES(VelocytoBase)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Velocyto Base computation finished.");
};

void SingleCellMultiomeItem::s_show_atac_landscape() {

	if (this->data()->species_ != soap::Species::Human) {
		G_WARN("Landscape now only support human.");
		return;
	}

	auto fragments = this->data()->fragments();

	if (fragments == nullptr) {
		G_WARN("Fragments has not been loaded.");
		return;
	}

	G_GETLOCK;

	auto&& metadata = this->data()->metadata()->mat_;

	auto factor_info = metadata.get_factor_information();

	if (factor_info.isEmpty()) {
		G_WARN("No Suitable Metadata for visualisation.");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"ATAC Landscape Plot Settings",
		{ "Choose Factor", "Choose Chromosome", "Normalize by", "Scale:yes", "Log Transform:yes"},
		{soap::InputStyle::FactorChoice, soap::InputStyle::SimpleChoice, 
		soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton},
		{soap::HumanChromosomeOrder, {"Sequencing Depth", "Choosed Chromosomes"}},
		{factor_info}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto [factor_name, levels] = factor_choice_to_pair(settings[0]);

	if (levels.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto chrs = simple_choice_to_list(settings[1]);

	if (chrs.isEmpty()) {
		G_UNLOCK;
		return;
	}
	QString normalize_method = settings[2];
	bool scale = switch_to_bool(settings[3]);
	bool log_transform = switch_to_bool(settings[4]);
		
	QStringList factors = metadata.get_qstring(factor_name);

	AtacLandscapePlotWorker* worker = new AtacLandscapePlotWorker(
		fragments,
		normalize_method,
		levels,
		factors,
		chrs,
		scale,
		log_transform
	);

	G_LINK_WORKER_THREAD(AtacLandscapePlotWorker, x_plot_ready, SingleCellMultiomeItem, s_receive_atac_landscape_plot);
};

void SingleCellMultiomeItem::s_receive_atac_landscape_plot(ATAC_LANDSCAPE_PLOT_ELEMENTS ele) {

	auto&& gs = this->draw_suite_->graph_settings_;

	int nrow = ele.mat.rows(), ncol = ele.mat.cols();

	QCustomPlot* draw_area = _Cp initialize_plot(gs);
	SoapTextElement* title = new SoapTextElement(draw_area, gs.get_title("ATAC Landscape"), gs.get_title_font());
	draw_area->plotLayout()->addElement(0, 0, title);

	QCPLayoutGrid* main_layout = new QCPLayoutGrid;
	draw_area->plotLayout()->addElement(1, 0, main_layout);

	QCPLayoutGrid* right_layout = new QCPLayoutGrid;

	draw_area->plotLayout()->addElement(1, 1, right_layout);

	auto legend_layout = _CpPatch set_legend_layout(draw_area, right_layout);

	if (ele.scale) {
		_CpPatch add_gradient_legend(
			draw_area,
			legend_layout,
			-1.0,
			1.0,
			"Accessibiliy",
			"Low",
			"High",
			gs.get_legend_title_font(),
			gs.get_legend_label_font(),
			_CpColor navy,
			Qt::white,
			_CpColor firebrick3
		);
	}
	else {

		_Cp add_gradient_legend(draw_area, legend_layout, ele.mat.minCoeff(), ele.mat.maxCoeff(), "Accessibility", gs);
	}

	QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	axis_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	axis_rect->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);
	QCPAxisRect* left_bottom_legend = new QCPAxisRect(draw_area, true);

	main_layout->addElement(0, 0, left_bottom_legend);
	main_layout->addElement(0, 1, axis_rect);

	auto cluster_names = _Cs sapply(ele.cell_loc, [](auto&& loc) {return std::get<0>(loc); });

	auto cluster_colors = gs.palette(cluster_names);
	int index = 0;
	for (const auto& [cluster, start, n] : ele.cell_loc) {

		double end = start + n;

		_CpPatch rectangle_borderless(
			draw_area, left_bottom_legend, 0, start, 1, end - start, cluster_colors[index++]
		);

		QCPItemText* label_x = new QCPItemText(draw_area);
		label_x->setClipToAxisRect(false);
		label_x->position->setAxisRect(left_bottom_legend);
		label_x->position->setAxes(left_bottom_legend->axis(QCPAxis::atBottom), left_bottom_legend->axis(QCPAxis::atLeft));
		label_x->position->setType(QCPItemPosition::ptPlotCoords);
		label_x->setPositionAlignment(Qt::AlignRight | Qt::AlignVCenter);
		label_x->position->setCoords(-1, (start + end) / 2);
		label_x->setText(cluster);
		label_x->setFont(gs.get_left_label_font());
	}

	_CpPatch remove_left_bottom_axis(left_bottom_legend);
	left_bottom_legend->setMinimumSize(_CpUtility get_max_text_width(_Cs sapply(ele.cell_loc, [](auto&& t) {return std::get<0>(t); }), gs.get_left_label_font()) * 1.2, 200);

	left_bottom_legend->axis(QCPAxis::atBottom)->setRange(-11, 1);
	left_bottom_legend->axis(QCPAxis::atLeft)->setRange(0, nrow);


	QCPColorMap* heatmap = new QCPColorMap(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
	heatmap->data()->setSize(ncol, nrow);
	heatmap->data()->setRange(QCPRange(0, ncol - 1), QCPRange(0, nrow - 1));
	_CpPatch remove_left_bottom_axis(axis_rect);
	_CpPatch set_range(axis_rect, QCPRange(-0.5, ncol - 0.5), QCPRange(-0.5, nrow - 0.5));

	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j) {
			heatmap->data()->setCell(j, i, ele.mat(i, j));
		}
	}
	heatmap->setInterpolate(false);
	heatmap->setTightBoundary(false);
	QCPColorGradient gradient;

	if (ele.scale) {
		gradient.setColorStopAt(0.5, Qt::white);
		gradient.setColorStopAt(1.0, _CpColor firebrick3);
		gradient.setColorStopAt(0.0, _CpColor navy);
	}
	else {
		gradient.setColorStopAt(0.5, gs.get_gradient_middle_color());
		gradient.setColorStopAt(1.0, gs.get_gradient_high_color());
		gradient.setColorStopAt(0.0, gs.get_gradient_low_color());
	}

	heatmap->setGradient(gradient);

	if (ele.scale) {
		heatmap->setDataRange({ -1, 1 });
	}
	else {
		heatmap->setDataRange({ ele.mat.minCoeff(), ele.mat.maxCoeff() });
	}

	QCPAxisRect* bottom_legend = new QCPAxisRect(draw_area, true);
	QCPAxisRect* bottom_left = new QCPAxisRect(draw_area, false);
	main_layout->addElement(1, 0, bottom_left);
	main_layout->addElement(1, 1, bottom_legend);

	auto chrColors = _CpUtility kmeans_palette(ele.chr_loc.size());
	index = 0;
	for (const auto& [chr, start, n] : ele.chr_loc) {

		double end = start + n;
		draw_area->addGraph(bottom_legend->axis(QCPAxis::atBottom), bottom_legend->axis(QCPAxis::atLeft));
		QCPCurve* legend = new QCPCurve(bottom_legend->axis(QCPAxis::atBottom), bottom_legend->axis(QCPAxis::atLeft));
		QVector<QCPCurveData> data(5);
		data[0] = QCPCurveData(0, start, 0);
		data[1] = QCPCurveData(1, start, 1);
		data[2] = QCPCurveData(2, end, 1);
		data[3] = QCPCurveData(3, end, 0);
		data[4] = QCPCurveData(0, start, 0);
		legend->setPen(Qt::NoPen);
		legend->setBrush(QBrush(chrColors[index++]));
		legend->data()->set(data, true);

		QCPItemText* label_x = new QCPItemText(draw_area);
		label_x->setClipToAxisRect(false);
		label_x->position->setAxisRect(bottom_legend);
		label_x->position->setAxes(bottom_legend->axis(QCPAxis::atBottom), bottom_legend->axis(QCPAxis::atLeft));
		label_x->position->setType(QCPItemPosition::ptPlotCoords);
		label_x->setPositionAlignment(Qt::AlignTop | Qt::AlignHCenter);
		label_x->position->setCoords((start + end) / 2, -(index % 2) * 1.5);
		label_x->setText(chr);
		label_x->setFont(gs.get_bottom_label_font());
	}

	bottom_legend->axis(QCPAxis::atBottom)->grid()->setVisible(false);
	bottom_legend->axis(QCPAxis::atLeft)->grid()->setVisible(false);
	bottom_legend->axis(QCPAxis::atBottom)->setBasePen(Qt::NoPen);
	bottom_legend->axis(QCPAxis::atBottom)->setTicks(false);
	bottom_legend->axis(QCPAxis::atLeft)->setBasePen(Qt::NoPen);
	bottom_legend->axis(QCPAxis::atLeft)->setTicks(false);

	bottom_legend->axis(QCPAxis::atBottom)->setRange(0, ncol);
	bottom_legend->axis(QCPAxis::atLeft)->setRange(-10, 1);

	bottom_legend->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);

	main_layout->setColumnStretchFactor(0, 1);
	main_layout->setColumnStretchFactor(1, 6);

	main_layout->setRowStretchFactor(0, 5);
	main_layout->setRowStretchFactor(1, 1);

	draw_area->plotLayout()->setRowStretchFactor(0, 1);
	draw_area->plotLayout()->setRowStretchFactor(1, 5);

	this->draw_suite_->update(draw_area);
};