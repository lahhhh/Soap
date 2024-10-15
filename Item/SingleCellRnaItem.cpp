#include "SingleCellRnaItem.h"

#include "Custom.h"

#include "ItemIOWorker.h"
#include "PcaWorker.h"
#include "UmapWorker.h"
#include "TsneWorker.h"
#include "LogNormalizeWorker.h"
#include "ScrubletWorker.h"
#include "DifferentialAnalysisWorker.h"
#include "IntegrateWorker.h"
#include "GseaWorker.h"
#include "HarmonyWorker.h"
#include "ScicnvWorker.h"
#include "InferCnvWorker.h"
#include "SlmWorker.h"
#include "LeidenPartitionWorker.h"
#include "VelocytoWorker.h"
#include "SvdWorker.h"
#include "CellTypeAnnotationWorker.h"
#include "Monocle3Worker.h"

#include "CustomPlot.h"
#include "FileIO.h"

#include "CommonDialog.h"
#include "StatisticsDialog.h"
#include "YesOrNoDialog.h"

void SingleCellRnaItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Metadata));
	this->check_variable(DATA_SUBMODULES(SparseInt));
	this->check_variable(DATA_SUBMODULES(SparseDouble));
	this->check_variable(DATA_SUBMODULES(Embedding));
	this->check_variable(DATA_SUBMODULES(DifferentialAnalysis));
	this->check_variable(DATA_SUBMODULES(GSEA));
	this->check_variable(DATA_SUBMODULES(CellChat));
	this->check_variable(DATA_SUBMODULES(CNV));
	this->check_variable(DATA_SUBMODULES(VelocytoBase));
	this->check_variable(DATA_SUBMODULES(Monocle3));

	auto item = new NoteItem(&this->data()->string_information_["Note"], this->data(), this->signal_emitter_);
	this->addChild(item);
};

void SingleCellRnaItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Statistics", s_statistics);

	// visualize
	ADD_MAIN_MENU("Visualize");

	ADD_ACTION("Quality View", "Visualize", s_view_quality);
	ADD_ACTION("Distribution Plot", "Visualize", s_distribution_plot);
	ADD_ACTION("Bubble Plot", "Visualize", s_bubble_plot);

	//sample
	ADD_MAIN_ACTION("Sample", s_sample);

	//setting
	ADD_MAIN_MENU("Settings");

	ADD_ACTION("Random Seed", "Settings", s_set_random_state);
	ADD_ACTION("Species", "Settings", s_set_species);

	// annotate
	ADD_MAIN_MENU("Annotate");

	ADD_ACTION("Add Metadata", "Annotate", s_add_metadata);
	ADD_ACTION("Edit Metadata", "Annotate", s_edit_metadata);
	ADD_ACTION("Annotate Cell Type", "Annotate", s_fast_annotation);

	//Doublet Detection
	ADD_MAIN_MENU("Detect Doublets");

	ADD_MENU("Detect Doublets | Scrublet", "Scrublet", "Detect Doublets");

	ADD_ACTION("Run Scrublet", "Detect Doublets | Scrublet", s_scrublet);
	ADD_ACTION("Set Threshold", "Detect Doublets | Scrublet", s_set_scrublet_threshold);

	//Filter
	ADD_MAIN_MENU("Filter");

	ADD_ACTION("By Features", "Filter", s_filter_by_features);
	ADD_ACTION("By Parameters", "Filter", s_filter_by_parameters);

	//Normalize
	ADD_MAIN_MENU("Normalize");

	ADD_MENU("Normalize | Log Normalize", "Log Normalize", "Normalize");

	ADD_ACTION("Default", "Normalize | Log Normalize", s_log_normalize_default);
	ADD_ACTION("Custom", "Normalize | Log Normalize", s_log_normalize_custom);

	// Dimensional Reduction
	ADD_MAIN_MENU("Dimension Reduction");

	ADD_MENU("Dimension Reduction | PCA", "PCA", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | PCA", s_pca_default);
	ADD_ACTION("Custom", "Dimension Reduction | PCA", s_pca_custom);

	ADD_MENU("Dimension Reduction | UMAP", "UMAP", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | UMAP", s_umap_default);
	ADD_ACTION("Custom", "Dimension Reduction | UMAP", s_umap_custom);

	ADD_MENU("Dimension Reduction | tSNE", "tSNE", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | tSNE", s_tsne_default);
	ADD_ACTION("Custom", "Dimension Reduction | tSNE", s_tsne_custom);


	//cluster
	ADD_MAIN_MENU("Cluster");

	ADD_MENU("Cluster | Louvain", "Louvain", "Cluster");

	ADD_MENU("Cluster | Louvain | Louvain", "Louvain", "Cluster | Louvain");

	ADD_ACTION("Default", "Cluster | Louvain | Louvain", s_louvain_default);
	ADD_ACTION("Custom", "Cluster | Louvain | Louvain", s_louvain_custom);

	ADD_MENU("Cluster | Louvain | Modified Louvain", "Modified Louvain", "Cluster | Louvain");

	ADD_ACTION("Default", "Cluster | Louvain | Modified Louvain", s_modified_louvain_default);
	ADD_ACTION("Custom", "Cluster | Louvain | Modified Louvain", s_modified_louvain_custom);

	ADD_MENU("Cluster | Louvain | Smart Local Moving", "Smart Local Moving", "Cluster | Louvain");

	ADD_ACTION("Default", "Cluster | Louvain | Smart Local Moving", s_smart_local_moving_default);
	ADD_ACTION("Custom", "Cluster | Louvain | Smart Local Moving", s_smart_local_moving_custom);

	ADD_MENU("Cluster | Leiden", "Leiden", "Cluster");

	ADD_ACTION("Default", "Cluster | Leiden", s_leiden_default);
	ADD_ACTION("Custom", "Cluster | Leiden", s_leiden_custom);

	// Lineage Tracing
	ADD_MAIN_MENU("Lineage Tracing");

	ADD_ACTION("Monocle3", "Lineage Tracing", s_monocle3);
	ADD_ACTION("Velocyto", "Lineage Tracing", s_velocyto);

	// DEG
	ADD_MAIN_MENU("Expression Analysis");

	ADD_ACTION("Find DEG", "Expression Analysis", s_find_deg);

	// Remove Batch Effect
	ADD_MAIN_MENU("Remove Batch Effect");

	ADD_ACTION("Harmony", "Remove Batch Effect", s_harmony);


	// GSEA
	ADD_MAIN_ACTION("GSEA", s_gsea);

	// CNV
	ADD_MAIN_MENU("CNV Detection");

	ADD_ACTION("SciCnv", "CNV Detection", s_scicnv);
	ADD_ACTION("InferCNV [No HMM]", "CNV Detection", s_infercnv);

	// integrate
	ADD_MAIN_ACTION("Integrate With...", s_integrate);

	// duplicate
	ADD_MAIN_ACTION("Duplicate", __s_duplicate);

	// save
	ADD_MAIN_ACTION("Save", __s_export_as_item);

	// delete
	ADD_MAIN_ACTION("Delete", __s_delete_this);

	// more
	ADD_MAIN_MENU("More...");

	ADD_ACTION("Recalculate quality parameters", "More...", s_recalculate_quality_parameters);

}

void SingleCellRnaItem::s_receive_fast_annotation(
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

void SingleCellRnaItem::s_fast_annotation() {
	G_GETLOCK;

	auto rna_counts = this->data()->counts();
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
			{ "Cell Type Name:Celltype-Soap", "Sub Type Name:Subtype-Soap", "Annotation By Cluster:yes", "Annotate By: " },
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

	G_LINK_WORKER_THREAD(CellTypeAnnotationWorker, x_annotation_ready, SingleCellRnaItem, s_receive_fast_annotation);
};

void SingleCellRnaItem::s_recalculate_quality_parameters() {
	G_GETLOCK;

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("No Counts data available.");
		return;
	}

	Eigen::ArrayXi col_count = custom::col_sum_mt(counts->mat_);
	const int ncol = counts->mat_.cols();
	Eigen::ArrayXi col_gene(ncol);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(ncol);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(ncol);

	for (int i = 0; i < ncol; ++i) {
		col_gene[i] = counts->mat_.outerIndexPtr()[i + 1] - counts->mat_.outerIndexPtr()[i];
	}

	QList<int> mitochondrial_location, ribosomal_location;

	auto& gene_symbols = counts->rownames_;
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
			mitochondrial_content += counts->mat_.row(i).cast<double>();
		}
		mitochondrial_content /= col_count.cast<double>();
	}

	if (ribosomal_location.length() > 0) {
		for (int& i : ribosomal_location) {
			ribosomal_content += counts->mat_.row(i).cast<double>();
		}
		ribosomal_content /= col_count.cast<double>();
	}

	custom::remove_na(mitochondrial_content);
	custom::remove_na(ribosomal_content);

	Metadata& metadata = *this->data()->metadata();
	metadata.mat_.update(METADATA_RNA_UMI_NUMBER, custom::cast<QVector>(col_count));
	metadata.mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, custom::cast<QVector>(col_gene));
	metadata.mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata.mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content));

	G_UNLOCK;
};

void SingleCellRnaItem::s_statistics() {

	StatisticsDialog::get_response(this->data());
};

void SingleCellRnaItem::s_sample() {
	G_GETLOCK;

	auto& metadata = this->data()->metadata()->mat_;
	QStringList factors = metadata.get_factor_name();

	if (factors.isEmpty()) {
		G_NOTICE("No Metadata for sampling.");
		G_UNLOCK;
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
	QStringList factor = metadata.get_qstring(feature);
	QStringList levels;
	if (metadata.data_type_[feature] == CustomMatrix::DataType::QStringFactor) {
		levels = metadata.string_factors_[feature];
	}
	else {
		levels = custom::cast<QString>(metadata.integer_factors_[feature]);
	}
	QVector<int> selected;
	if (downsample != "ALL") {
		int downsample_number = downsample.toInt();
		for (const auto& level : levels) {
			auto index = custom::match(factor, level);
			if (index.size() > downsample_number) {
				index = custom::sample(index, downsample_number, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		selected = custom::seq_n(0, factor.size());
	}

	this->signal_emitter_->x_data_create_soon(this->data()->col_reordered(selected), soap::VariableType::SingleCellRna, "Sampled SingleCellRna");

	G_UNLOCK;

};

void SingleCellRnaItem::__s_rename() {
	G_GETLOCK;

	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Rename SingleCellRna",
		{ "New Name:" + this->title_, "Overwrite Metadata:Source:yes" },
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
		this->data()->metadata()->mat_.update("Source", QStringList(this->data()->counts()->mat_.cols(), this->title_), CustomMatrix::DataType::QStringFactor);
	}
};

void SingleCellRnaItem::slice(const Eigen::ArrayX<bool>& row_slice, const Eigen::ArrayX<bool>& col_slice) {

	this->__clear_reserve_data();

	this->data()->slice(row_slice, col_slice);

	this->signal_emitter_->update_information(this->title_, soap::VariableType::SingleCellRna, this->data(), true);

	this->__check_data();

	G_LOG("Filter finished.");
};

void SingleCellRnaItem::col_slice(const Eigen::ArrayX<bool>& col_slice) {

	this->__clear_reserve_data();

	this->data()->col_slice(col_slice);

	this->signal_emitter_->update_information(this->title_, soap::VariableType::SingleCellRna, this->data(), true);

	this->__check_data();

	G_LOG("Filter finished.");
};

void SingleCellRnaItem::s_receive_normalize(SparseDouble* data) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_NORMALIZED);

	DATA_SUBMODULES(SparseDouble)[title] = std::move(*data);
	delete data;

	auto item = new SparseDoubleItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(SparseDouble)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Normalize finished.");
};

void SingleCellRnaItem::s_log_normalize_default() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Raw counts data missed.");
		return;
	}
	G_GETLOCK;

	if (auto item = this->normalized()) {
		item->__remove_this();
	}

	G_LOG("Log normalize start...");
	LogNormalizeWorker* worker = new LogNormalizeWorker(counts, 10000.0);
	G_LINK_WORKER_THREAD(LogNormalizeWorker, x_log_normalize_ready, SingleCellRnaItem, s_receive_normalize)
};

void SingleCellRnaItem::s_log_normalize_custom() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Raw counts data missed.");
		return;
	}
	G_GETLOCK;

	QStringList ret = CommonDialog::get_response(
		this->signal_emitter_,
		"Set a Scale Factor",
		{ "Scale Factor(>0):10000" },
		{ soap::InputStyle::NumericLineEdit}
	);

	if (ret.isEmpty()) {
		G_UNLOCK;
		return;
	}

	double factor = ret[0].toDouble();
	if (factor <= 0.0) {
		G_WARN("Scale factor cannot be negative.");
		G_UNLOCK;
		return;
	}

	if (auto item = this->normalized()) {
		item->__remove_this();
	}

	G_LOG("Log normalize start...");
	LogNormalizeWorker* worker = new LogNormalizeWorker(counts, factor);
	G_LINK_WORKER_THREAD(LogNormalizeWorker, x_log_normalize_ready, SingleCellRnaItem, s_receive_normalize)
};

void SingleCellRnaItem::s_filter_by_features() {
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
	if (in_place) {
		G_LOG("Filter Object by Features...");
		this->col_slice(filter);
	}
	else {
		G_LOG("Filter Object and Create New Object by Features...");
		this->signal_emitter_->x_data_create_soon(this->data()->col_sliced(filter), soap::VariableType::SingleCellRna, "Sliced SingleCellRna");
	}

	G_UNLOCK;
}

void SingleCellRnaItem::s_filter_by_parameters() {
	G_GETLOCK;

	if (this->data()->counts() == nullptr) {
		G_WARN("RNA Count data is missed.");
		return;
	}

	auto& metadata = this->data()->metadata()->mat_;
	if (!metadata.contains(METADATA_RNA_UMI_NUMBER) ||
		metadata.data_type_[METADATA_RNA_UMI_NUMBER] != CustomMatrix::DataType::IntegerNumeric) {

		G_NOTICE("Data of UMI number sums is not available, this setting of minimum/maximum count will be not applied.")

	}

	if (!metadata.contains(METADATA_RNA_UNIQUE_GENE_NUMBER) ||
		metadata.data_type_[METADATA_RNA_UNIQUE_GENE_NUMBER] != CustomMatrix::DataType::IntegerNumeric) {
		G_NOTICE("Data of Gene number sums is not available, this setting of minimum/maximum unique gene number will be not applied.")

	}

	if (!metadata.contains(METADATA_RNA_MITOCHONDRIAL_CONTENT) ||
		metadata.data_type_[METADATA_RNA_MITOCHONDRIAL_CONTENT] != CustomMatrix::DataType::DoubleNumeric) {
		G_NOTICE("Data of Mitochondrial Content is not available, this setting of Mitochondrial Content will be not applied.")

	}
	QStringList parameters = CommonDialog::get_response(
		this->signal_emitter_,
		"Set filter parameters",
		{ "Max Count:30000", "Max Unique Gene:10000", "Min Count:1000",
		"Min Unique Gene:500", "Max Mitochondrial Content(Percentage):10", "Min Cell Per Gene:5", "In place:no" },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::SwitchButton}
	);

	if (parameters.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto filter = this->get_parameter_filter(parameters.sliced(0, 6));
	if (filter.first.count() == 0 || filter.second.count() == 0) {
		G_WARN("No Data remained after filtering!");
		G_UNLOCK;
		return;
	}
	bool in_place = switch_to_bool(parameters[6]);
	if (in_place) {
		G_LOG("Filter Object by Parameters...");
		this->slice(filter.first, filter.second);
	}
	else {
		G_LOG("Filter Object and Create New Object by Parameters...");
		this->signal_emitter_->x_data_create_soon(this->data()->sliced(filter.first, filter.second), soap::VariableType::SingleCellRna, "Sliced SingleCellRna");
	}

	G_UNLOCK;
}

std::pair<Eigen::ArrayX<bool>, Eigen::ArrayX<bool> > SingleCellRnaItem::get_parameter_filter(const QStringList& parameters) {

	const int nrow = this->data()->counts()->mat_.rows();
	const int ncol = this->data()->counts()->mat_.cols();
	Eigen::ArrayX<bool> passed_column = Eigen::ArrayX<bool>::Constant(ncol, true);
	Eigen::ArrayX<bool> passed_row = Eigen::ArrayX<bool>::Constant(nrow, true);

	auto& metadata = this->data()->metadata()->mat_;

	if (!parameters[0].isEmpty() || !parameters[2].isEmpty()) {
		if (metadata.contains(METADATA_RNA_UMI_NUMBER) &&
			metadata.data_type_[METADATA_RNA_UMI_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {

			auto& n_umi = metadata.get_const_integer_reference(METADATA_RNA_UMI_NUMBER);
			if (!parameters[0].isEmpty()) {
				int maximum_count_number = parameters[0].toInt();
				passed_column *= custom::less_equal(n_umi, maximum_count_number);
			}
			if (!parameters[2].isEmpty()) {
				int minimum_count_number = parameters[2].toInt();
				passed_column *= custom::greater_equal(n_umi, minimum_count_number);
			}

		}
	}

	if (!parameters[1].isEmpty() || !parameters[3].isEmpty()) {
		if (metadata.contains(METADATA_RNA_UNIQUE_GENE_NUMBER) &&
			metadata.data_type_[METADATA_RNA_UNIQUE_GENE_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {

			auto& n_gene = metadata.get_const_integer_reference(METADATA_RNA_UNIQUE_GENE_NUMBER);

			if (!parameters[1].isEmpty()) {
				int maximum_gene_number = parameters[1].toInt();
				passed_column *= custom::less_equal(n_gene, maximum_gene_number);
			}
			if (!parameters[3].isEmpty()) {
				int minimum_gene_number = parameters[3].toInt();
				passed_column *= custom::greater_equal(n_gene, minimum_gene_number);
			}

		}
	}

	if (!parameters[4].isEmpty()) {
		if (metadata.contains(METADATA_RNA_MITOCHONDRIAL_CONTENT) &&
			metadata.data_type_[METADATA_RNA_MITOCHONDRIAL_CONTENT] == CustomMatrix::DataType::DoubleNumeric) {

			auto& mitochondrial_content = metadata.get_const_double_reference(METADATA_RNA_MITOCHONDRIAL_CONTENT);
			double maximum_mitochondrial_content = parameters[4].toDouble() / 100;
			passed_column *= custom::less_equal(mitochondrial_content, maximum_mitochondrial_content);
		}
	}

	if (!parameters[5].isEmpty()) {
		passed_row *= custom::row_count<int, true, false>(this->data()->counts()->mat_, 0.) >= parameters[5].toInt();
	}

	return std::make_pair(passed_row, passed_column);
}

void SingleCellRnaItem::s_pca_custom() {
	G_GETLOCK;
	if (this->data()->counts() == nullptr) {
		G_WARN("Raw RNA counts data missed!");
		G_UNLOCK;
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"PCA settings",
		{ "Number of dimension:50", "Random State:" + QString::number(this->data()->random_state_),
		"Variable Feature Proportion:0.1", "Fixed Variable Feature Number?:yes", "Fixed Variable Feature Number:2000" },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::NumericLineEdit, soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	int n_dimension = settings[0].toInt();
	if (n_dimension < 10 || n_dimension > 100) {
		G_WARN("Number of dimension should be between 10 and 100.")
			G_UNLOCK;
		return;
	}

	int random_state = settings[1].toInt();
	bool variable_feature_number_fixed = switch_to_bool(settings[3]);
	double variable_proportion = 0.;
	int variable_number = 0;
	if (variable_feature_number_fixed) {
		variable_number = settings[4].toInt();
		if (variable_number < 1000) {
			G_WARN("Number of Variable Features shold be more than 1000.")
				G_UNLOCK;
			return;
		}
	}
	else {
		variable_proportion = settings[2].toDouble();
		if (variable_proportion < 0.01 || variable_proportion > 1) {
			G_WARN("Number of Variable Features Proportion should be between 0.01 and 1.")
				G_UNLOCK;
			return;
		}
	}

	G_LOG("PCA start...");

	PcaWorker* worker = new PcaWorker(
		&this->data()->counts()->mat_, 
		variable_feature_number_fixed ? variable_number : variable_proportion, 
		n_dimension,
		random_state);
	G_LINK_WORKER_THREAD(PcaWorker, x_pca_ready, SingleCellRnaItem, s_receive_pca)
};

void SingleCellRnaItem::s_pca_default() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Raw counts data missed!");
		return;
	}

	G_GETLOCK;

	if (auto item = this->pca()) {
		item->__remove_this();
	}

	G_LOG("PCA by tSVD start...");

	PcaWorker* worker = new PcaWorker(&counts->mat_, 2000, 50, this->data()->random_state_);

	G_LINK_WORKER_THREAD(PcaWorker, x_pca_ready, SingleCellRnaItem, s_receive_pca);
}

void SingleCellRnaItem::s_receive_pca(Eigen::MatrixXd mat, QVector<double> sd) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_PCA);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Pca,
		mat,
		this->data()->counts()->colnames_,
		custom::paste("PCA-", custom::cast<QString>(custom::seq_n(1, mat.cols())))
	);

	this->data()->double_vectors_[VARIABLE_RNA_PCA_STANDARD_DEVIATION] = sd;

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		this->data()->pca(),
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("PCA finished.");
}

void SingleCellRnaItem::s_harmony() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}
	G_GETLOCK;

	if (auto item = this->harmony()) {
		item->__remove_this();
	}

	QStringList factor_names = this->data()->metadata()->mat_.get_factor_name(false);

	if (factor_names.isEmpty()) {
		G_LOG("No suitable feature for Harmony");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Harmony Settings",
		{ "Batch Factor" , "Dimension Start:1", "Dimension End:50" },
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit },
		{ factor_names }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	factor_names = simple_choice_to_list(settings[0]);
	if (factor_names.isEmpty()) {
		G_UNLOCK;
		return;
	}

	factor_names = custom::unique(factor_names);
	QList<QStringList> factor_list;
	for (const auto& factor_name : factor_names) {
		factor_list << this->data()->metadata()->mat_.get_qstring(factor_name);
	}

	int dim_start = settings[1].toInt();
	int dim_end = settings[2].toInt();

	int ndim = pca->data_.mat_.cols();

	if (dim_start < 1 || dim_start > ndim || dim_end < 1 || dim_end > ndim || dim_start >= dim_end) {
		G_WARN("Illegal Dimension Setting.");
		G_UNLOCK;
		return;
	}

	HarmonyWorker* worker = new HarmonyWorker(
		pca->data_.mat_.block(0, dim_start - 1, pca->data_.mat_.rows(), dim_end - dim_start + 1),
		factor_list
	);
	G_LINK_WORKER_THREAD(HarmonyWorker, x_harmony_ready, SingleCellRnaItem, s_receive_harmony)
};

void SingleCellRnaItem::s_receive_harmony(Eigen::MatrixXd mat) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_HARMONY);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Harmony,
		mat,
		this->data()->counts()->colnames_,
		custom::paste("Harmony-", custom::cast<QString>(custom::seq_n(1, mat.cols()))));

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		this->data()->harmony(),
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Harmony finished");
};

void SingleCellRnaItem::s_tsne_default() {
	G_GETLOCK;

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("Data requires PCA.");
		G_UNLOCK;
		return;
	}

	if (auto item = this->tsne()) {
		item->__remove_this();
	}

	int start_dimension = 0;
	int last_dimension = pca->data_.mat_.cols() - 1;
	if (last_dimension > 19) {
		last_dimension = 19;
	}
	G_LOG("Default tSNE start, using 1~" + QString::number(last_dimension + 1) + " PCA components...");

	Eigen::MatrixXd tsne_input = pca->data_.mat_.block(0, start_dimension, pca->data_.mat_.rows(), last_dimension - start_dimension + 1);
	
	TsneWorker* worker = new TsneWorker(tsne_input, this->data()->random_state_);

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, SingleCellRnaItem, s_receive_tsne)
};

void SingleCellRnaItem::s_tsne_custom() {

	G_GETLOCK;

	QStringList from_list;

	auto pca = this->data()->pca();
	if (pca != nullptr) {
		from_list << "PCA";
	}

	auto harmony = this->data()->harmony();
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No PCA or Harmony result for tSNE");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"tSNE settings",
		{
			"Dimension start:1",
			"Dimension end:20",
			"from",
			"Random State:" + QString::number(this->data()->random_state_)
		},
		{
		soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox,
			soap::InputStyle::IntegerLineEdit
		},
		{ from_list }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	Eigen::MatrixXd* mat;
	QString from = settings[2];
	if (from == "PCA") {
		mat = &pca->data_.mat_;
	}
	else {
		mat = &harmony->data_.mat_;
	}
	const int nrow = mat->rows(), ncol = mat->cols();

	int start_dimension = settings[0].toInt();
	if (start_dimension < 1 || start_dimension > ncol - 1) {
		G_WARN("Illegal start dimension value. Reset to 1.");
		start_dimension = 1;
	}

	int last_dimension = settings[1].toInt();
	if (last_dimension <= start_dimension || last_dimension > ncol) {
		G_WARN("Illegal end dimension value. Please reset.");
		G_UNLOCK;
		return;
	}

	unsigned int random_state = settings[3].toUInt();

	if (auto item = this->tsne()) {
		item->__remove_this();
	}

	Eigen::MatrixXd tsne_input = mat->block(0, start_dimension - 1, nrow, last_dimension - start_dimension + 1);

	G_LOG("Customized tSNE start using " + from +
		" (" + QString::number(start_dimension) + "~" + QString::number(last_dimension) + ") dimensions.");


	TsneWorker* worker = new TsneWorker(tsne_input, random_state);

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, SingleCellRnaItem, s_receive_tsne)

};

void SingleCellRnaItem::s_receive_tsne(Eigen::MatrixXd mat) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_TSNE);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Tsne,
		mat,
		this->data()->counts()->colnames_,
		custom::paste("tSNE-", custom::cast<QString>(custom::seq_n(1, mat.cols()))));

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		this->data()->tsne(),
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("t-SNE finished.");
};

void SingleCellRnaItem::s_umap_default() {
	G_GETLOCK;

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("Data requires PCA.");
		G_UNLOCK;
		return;
	}

	if (auto item = this->umap()) {
		item->__remove_this();
	}

	int start_dimension = 0;
	int last_dimension = pca->data_.mat_.cols() - 1;
	if (last_dimension > 19) {
		last_dimension = 19;
	}
	G_LOG("Default UMAP start, using 1~" + QString::number(last_dimension + 1) + " PCA components...");

	Eigen::MatrixXd umap_input = pca->data_.mat_.block(0, start_dimension, pca->data_.mat_.rows(), last_dimension - start_dimension + 1);

	UmapWorker* worker = new UmapWorker(
		umap_input,
		30,
		"Angular",
		1.0,
		"Random",
		0.3,
		1.0,
		1.0,
		1.0,
		5,
		this->data()->random_state_,
		50
	);

	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, SingleCellRnaItem, s_receive_umap)
}

void SingleCellRnaItem::s_umap_custom() {

	G_GETLOCK;

	QStringList from_list;

	auto pca = this->data()->pca();
	if (pca != nullptr) {
		from_list << "PCA";
	}

	auto harmony = this->data()->harmony();
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No PCA or Harmony result for UMAP");
		G_UNLOCK;
		return;
	}

	/*
		parameters setting
	*/
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"UMAP settings",
		{
			"Dimension start:1",
			"Dimension end:20",
			"Number of neighbors:30",
			"Distance metric",
			"Learning rate(>0):1.0",
			"Init type",
			"Minimum distance(>0):0.3",
			"Spread:1.0",
			"Set op mix ratio(0~1):1.0",
			"Repulsion strength(>0):1.0",
			"Negative sample rate(>0):5",
			"Random state:" + QString::number(this->data()->random_state_),
			"From"
		},
		{
		soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox,
			soap::InputStyle::NumericLineEdit,
			soap::InputStyle::ComboBox,
			soap::InputStyle::NumericLineEdit,
			soap::InputStyle::NumericLineEdit,
			soap::InputStyle::NumericLineEdit,
			soap::InputStyle::NumericLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox
	},
		{
			{ "Angular", "Euclidean", "Manhattan" },
				{ "Random" },
				from_list
		}
		);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	Eigen::MatrixXd* mat;
	QString from = settings[12];
	if (from == "PCA") {
		mat = &pca->data_.mat_;
	}
	else {
		mat = &harmony->data_.mat_;
	}
	const int nrow = mat->rows(), ncol = mat->cols();

	int start_dimension = settings[0].toInt();
	if (start_dimension < 1 || start_dimension > ncol - 1) {
		G_WARN("Illegal start dimension value. Reset to 1.");
		start_dimension = 1;
	}

	int last_dimension = settings[1].toInt();
	if (last_dimension <= start_dimension || last_dimension > ncol) {
		G_WARN("Illegal end dimension value. Please reset.");
		G_UNLOCK;
		return;
	}

	int n_neighbors = settings[2].toInt();
	if (n_neighbors < 2 || n_neighbors > ncol || n_neighbors > 100) {
		G_WARN("Illegal neighbor value : " + QString::number(n_neighbors) + ".Please reset.");
		G_UNLOCK;
		return;
	}

	QString metric = settings[3];

	double learning_rate = settings[4].toDouble();
	if (learning_rate <= 0) {
		G_WARN("Learning rate must be positive. Reset to 1.0");
		learning_rate = 1.0;
	}

	QString initialize_type = settings[5];

	double minimum_distance = settings[6].toDouble(), spread = settings[7].toDouble();
	if (minimum_distance < 0 || minimum_distance > spread) {
		G_LOG("Minimum distance must be less than or equal to spread and can not be negative. Please reset");
		G_UNLOCK;
		return;
	}

	double set_op_mix_ratio = settings[8].toDouble();
	if (set_op_mix_ratio < 0 || set_op_mix_ratio > 1.0) {
		G_LOG("Set op mix ratio must be between 0.0 and 1.0. Reset to 1.0");
		set_op_mix_ratio = 1.0;
	}

	double gamma = settings[9].toDouble();
	if (gamma < 0) {
		G_LOG("Repulsion strength cannot be negative. Reset to 1.0");
		gamma = 1.0;
	}

	int negative_sample_rate = settings[10].toInt();
	if (negative_sample_rate <= 0) {
		G_LOG("Negative sample rate must be positive. Reset to 5.");
		negative_sample_rate = 5;
	}
	else if (negative_sample_rate > 10) {
		G_LOG("Larger negative sample rate may cost more time.");
	}

	if (auto item = this->umap()) {
		item->__remove_this();
	}

	int random_state = settings[11].toInt();
	Eigen::MatrixXd umap_input = mat->block(0, start_dimension - 1, nrow, last_dimension - start_dimension + 1);

	G_LOG("Customized UMAP start using " + from +
		" (" + QString::number(start_dimension) + "~" + QString::number(last_dimension) + ") dimensions.");

	UmapWorker* worker = new UmapWorker(
		umap_input,
		n_neighbors,
		metric,
		learning_rate,
		initialize_type,
		minimum_distance,
		spread,
		set_op_mix_ratio,
		gamma,
		negative_sample_rate,
		random_state,
		50
	);

	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, SingleCellRnaItem, s_receive_umap)
}

void SingleCellRnaItem::s_receive_umap(Eigen::MatrixXd mat) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_UMAP);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Umap,
		mat, 
		this->data()->counts()->colnames_,
		custom::paste("UMAP-", custom::cast<QString>(custom::seq_n(1, mat.cols()))));

	auto item = new EmbeddingItem(
		title, 
		this->index_tree_, 
		this->data()->umap(), 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("UMAP finished.");
}

void SingleCellRnaItem::s_leiden_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	LeidenPartitionWorker* worker = new LeidenPartitionWorker(pca->data_.mat_, "Modularity", "Euclidean", 30, 50, 0.6);
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, SingleCellRnaItem, s_receive_leiden);
}

void SingleCellRnaItem::s_louvain_default() {
	
	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "Louvain", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellRnaItem, s_receive_louvain);
}

void SingleCellRnaItem::s_modified_louvain_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;
	
	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "Modified Louvain", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellRnaItem, s_receive_modified_louvain);
}

void SingleCellRnaItem::s_smart_local_moving_default() {
	
	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "SLM", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellRnaItem, s_receive_slm);
}

void SingleCellRnaItem::s_smart_local_moving_custom() {

	QStringList from_list;

	auto pca = this->data()->pca();
	if (pca != nullptr) {
		from_list << "PCA";
	}

	auto harmony = this->data()->harmony();
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No PCA or Harmony result for Cluster");
		G_UNLOCK;
		return;
	}

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Louvain Clustering Settings",

		{ 
			"Resolution:0.6",
			"Random State:" + QString::number(this->data()->random_state_),
			"Number of Neighbors:30",
			"Number of Trees:50", 
			"Metric", 
			"From",
			"Number of Start:10", 
			"Number of Iterations:10", 
			"Modularity Function" 
		},
		{ 
			soap::InputStyle::NumericLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ComboBox, 
			soap::InputStyle::ComboBox,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ComboBox
		},

		{ { "Euclidean", "Angular", "Manhattan" }, from_list, { "1", "2" }}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	double resolution = settings[0].toDouble();
	if (resolution <= 0 || resolution > 4) {
		G_WARN("Illegal resolution! reset to 0.6");
		resolution = 0.6;
	}

	int random_state = settings[1].toInt();
	int n_neighbors = settings[2].toInt();

	int n_trees = settings[3].toInt();
	if (n_trees < 10 || n_trees > 100) {
		G_WARN("Number of trees should be betweem 10 and 100, reset to 50");
		n_trees = 50;
	}

	QString metric = settings[4];

	Eigen::MatrixXd* from_matrix = nullptr;
	QString from = settings[5];
	if (from == "PCA") {
		from_matrix = &pca->data_.mat_;
	}
	else {
		from_matrix = &harmony->data_.mat_;
	}

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	int n_start = settings[6].toInt();
	if (n_start < 1 || n_start > 100) {
		G_WARN("Number of start should be between 1 and 100, reset to 10");
		n_start = 10;
	}

	int n_iteration = settings[7].toInt();
	if (n_iteration < 1 || n_iteration > 100) {
		G_WARN("Number of iteration should be between 1 and 100, reset to 10");
		n_iteration = 10;
	}

	int modularity_function = settings[8].toInt();

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(*from_matrix, "SLM", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellRnaItem, s_receive_slm);
}

void SingleCellRnaItem::s_modified_louvain_custom() {

	QStringList from_list;

	auto pca = this->data()->pca();
	if (pca != nullptr) {
		from_list << "PCA";
	}

	auto harmony = this->data()->harmony();
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No PCA or Harmony result for Cluster");
		G_UNLOCK;
		return;
	}

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Louvain Clustering Settings",
		{  
			"Resolution:0.6", 
			"Random State:" + QString::number(this->data()->random_state_),
			"Number of Neighbors:30",
			"Number of Trees:50", 
			"Metric", 
			"From",
			"Number of Start:10", 
			"Number of Iterations:10", 
			"Modularity Function" 
		},
		{ 
			soap::InputStyle::NumericLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox, 
			soap::InputStyle::ComboBox,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ComboBox
		},
		{ { "Euclidean", "Angular", "Manhattan" }, from_list, { "1", "2" }}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	double resolution = settings[0].toDouble();
	if (resolution <= 0 || resolution > 4) {
		G_WARN("Illegal resolution! reset to 0.6");
		resolution = 0.6;
	}

	int random_state = settings[1].toInt();
	int n_neighbors = settings[2].toInt();

	int n_trees = settings[3].toInt();
	if (n_trees < 10 || n_trees > 100) {
		G_WARN("Number of trees should be betweem 10 and 100, reset to 50");
		n_trees = 50;
	}

	QString metric = settings[4];

	Eigen::MatrixXd* from_matrix = nullptr;
	QString from = settings[5];
	if (from == "PCA") {
		from_matrix = &pca->data_.mat_;
	}
	else {
		from_matrix = &harmony->data_.mat_;
	}

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	int n_start = settings[6].toInt();
	if (n_start < 1 || n_start > 100) {
		G_WARN("Number of start should be between 1 and 100, reset to 10");
		n_start = 10;
	}

	int n_iteration = settings[7].toInt();
	if (n_iteration < 1 || n_iteration > 100) {
		G_WARN("Number of iteration should be between 1 and 100, reset to 10");
		n_iteration = 10;
	}

	int modularity_function = settings[8].toInt();

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(*from_matrix, "Modified Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellRnaItem, s_receive_modified_louvain);
}

void SingleCellRnaItem::s_leiden_custom() {

	QStringList from_list;

	auto pca = this->data()->pca();
	if (pca != nullptr) {
		from_list << "PCA";
	}

	auto harmony = this->data()->harmony();
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No PCA or Harmony result for Cluster");
		G_UNLOCK;
		return;
	}
	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Leiden Clustering Settings",
		{ "Method", "Metric", "Resolution:0.6", "Number of Neighbors:30", "Number of Trees:50", "From" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit
		, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit, soap::InputStyle::ComboBox},
		{ { "Modularity", "CPM", "RBConfiguration", "RBER", "Significance", "Surprise" }
		, { "Euclidean", "Angular", "Manhattan" }, from_list }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString method = settings[0];
	QString metric = settings[1];

	double resolution = settings[2].toDouble();
	if (resolution <= 0 || resolution > 4) {
		G_WARN("Illegal resolution! reset to 0.6");
		resolution = 0.6;
	}

	int n_neighbors = settings[3].toInt();

	int n_trees = settings[4].toInt();
	if (n_trees < 10 || n_trees > 100) {
		G_WARN("Number of trees should be betweem 10 and 100, reset to 50");
		n_trees = 50;
	}

	Eigen::MatrixXd* from_matrix = nullptr;
	QString from = settings[5];
	if (from == "PCA") {
		from_matrix = &pca->data_.mat_;
	}
	else /* if(from == "Harmony") */ {
		from_matrix = &harmony->data_.mat_;
	}

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	G_LOG("Start Clustering...");

	LeidenPartitionWorker* worker = new LeidenPartitionWorker(*from_matrix, method, metric, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, SingleCellRnaItem, s_receive_leiden);
}

void SingleCellRnaItem::s_louvain_custom() {

	QStringList from_list;

	auto pca = this->data()->pca();
	if (pca != nullptr) {
		from_list << "PCA";
	}

	auto harmony = this->data()->harmony();
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No PCA or Harmony result for Cluster");
		G_UNLOCK;
		return;
	}

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Louvain Clustering Settings",

		{ 
			"Resolution:0.6", 
			"Random State:" + QString::number(this->data()->random_state_),
			"Number of Neighbors:30",
			"Number of Trees:50", 
			"Metric", 
			"From",
			"Number of Start:10", 
			"Number of Iterations:10", 
			"Modularity Function" 
		},
		{ 
			soap::InputStyle::NumericLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ComboBox, 
			soap::InputStyle::ComboBox,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ComboBox
		},
		{ { "Euclidean", "Angular", "Manhattan" }, from_list, { "1", "2" }}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	double resolution = settings[0].toDouble();
	if (resolution <= 0 || resolution > 4) {
		G_WARN("Illegal resolution! reset to 0.6");
		resolution = 0.6;
	}

	int random_state = settings[1].toInt();
	int n_neighbors = settings[2].toInt();

	int n_trees = settings[3].toInt();
	if (n_trees < 10 || n_trees > 100) {
		G_WARN("Number of trees should be betweem 10 and 100, reset to 50");
		n_trees = 50;
	}

	QString metric = settings[4];

	Eigen::MatrixXd* from_matrix = nullptr;
	QString from = settings[5];
	if (from == "PCA") {
		from_matrix = &pca->data_.mat_;
	}
	else {
		from_matrix = &harmony->data_.mat_;
	}

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	int n_start = settings[6].toInt();
	if (n_start < 1 || n_start > 100) {
		G_WARN("Number of start should be between 1 and 100, reset to 10");
		n_start = 10;
	}

	int n_iteration = settings[7].toInt();
	if (n_iteration < 1 || n_iteration > 100) {
		G_WARN("Number of iteration should be between 1 and 100, reset to 10");
		n_iteration = 10;
	}

	int modularity_function = settings[8].toInt();

	G_LOG("Start Clustering...");
	SlmWorker* worker = new SlmWorker(*from_matrix, "Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellRnaItem, s_receive_louvain);
}

void SingleCellRnaItem::s_receive_leiden(QVector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_RNA_LEIDEN_CLUSTER, cluster, CustomMatrix::DataType::IntegerFactor);
};

void SingleCellRnaItem::s_receive_louvain(std::vector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_RNA_LOUVAIN_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void SingleCellRnaItem::s_receive_modified_louvain(std::vector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_RNA_MODIFIED_LOUVAIN_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void SingleCellRnaItem::s_receive_slm(std::vector<int> cluster) {


	this->data()->metadata()->mat_.update(METADATA_RNA_SLM_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void SingleCellRnaItem::s_set_scrublet_threshold() {

	if (!this->data()->double_vectors_.contains(VECTOR_SCRUBLET_SCORES) ||
		!this->data()->double_vectors_.contains(VECTOR_SCRUBLET_SCORES_SIMULATED)) {

		G_WARN("Scrublet is not performed");
	}

	int n_cell = this->data()->counts()->mat_.cols();

	QStringList ret = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Threshold for Scrublet",
		{ "Threshold (0~1)" },
		{ soap::InputStyle::NumericLineEdit}
	);
	if (ret.isEmpty())return;

	double threshold = ret[0].toDouble();
	if (threshold <= 0.0 || threshold >= 1.0) {
		G_LOG("Please set a threshold between 0 and 1.");
		return;
	}
	QVector<double> scores = this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES];

	if (n_cell != scores.size()) {
		G_WARN("Invalid Scrublet Data.");
		return;
	}
	QStringList labels;
	for (int i = 0; i < n_cell; ++i) {
		labels << (scores[i] < threshold ? "Singlet" : "Doublet");
	}

	Eigen::ArrayXd original_score = custom::cast<Eigen::ArrayX>(scores);
	Eigen::ArrayXd simulate_score = custom::cast<Eigen::ArrayX>(this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES_SIMULATED]);

	auto [counts, edges, locs] = custom::histogram(original_score, 32);
	Eigen::ArrayXd normed_counts = log10(counts.cast<double>() + 1);

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);
	custom_plot::bar_plot(draw_area, axis_rect, Qt::blue, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);
	custom_plot::patch::line(draw_area, axis_rect, QVector<double>{threshold, threshold}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
	
	draw_area->plotLayout()->insertRow(0);
	SoapTextElement* title = new SoapTextElement(draw_area, "Observed transcriptomes", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(0, 0, title);

	std::tie(counts, edges, locs) = custom::histogram(simulate_score, 32);
	normed_counts = log10(counts.cast<double>() + 1);
	axis_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(2, 0, axis_rect);
	custom_plot::bar_plot(draw_area, axis_rect, Qt::red, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);
	custom_plot::patch::line(draw_area, axis_rect, QVector<double>{threshold, threshold}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
	

	draw_area->plotLayout()->insertRow(2);
	SoapTextElement* title2 = new SoapTextElement(draw_area, "Simulated doublets", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(2, 0, title2);

	this->draw_suite_->update(draw_area);

	this->data()->metadata()->mat_.update(METADATA_SCRUBLET_LABELS, labels, CustomMatrix::DataType::QStringFactor);
	G_LOG("Threshold : " + QString::number(threshold) + " has been set for Scrublet.");

	double original_rate = custom::greater_than(original_score, threshold).count() / (double)original_score.size();
	double simulate_rate = custom::greater_than(simulate_score, threshold).count() / (double)simulate_score.size();

	G_LOG("Detected " + QString::number(original_rate * 100) + " % doublets in original data and " + QString::number(simulate_rate * 100) + " % doublets in simulated doublets.");
};

void SingleCellRnaItem::s_receive_scrublet(const Eigen::ArrayXd& original_score, const Eigen::ArrayXd& simulate_score) {
	this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES] = custom::cast<QVector>(original_score);
	this->data()->metadata()->mat_.update(METADATA_SCRUBLET_SCORES, this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES]);
	this->data()->double_vectors_[VECTOR_SCRUBLET_SCORES_SIMULATED] = custom::cast<QVector>(simulate_score);
	auto threshold = custom::threshold_minimum(simulate_score);

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

	auto [counts, edges, locs] = custom::histogram(original_score, 32);
	Eigen::ArrayXd normed_counts = log10(counts.cast<double>() + 1);


	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);
	custom_plot::bar_plot(draw_area, axis_rect, Qt::blue, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);

	if (threshold.first) {
		double threshold_value = threshold.second;
		custom_plot::patch::line(draw_area, axis_rect, QVector<double>{threshold_value, threshold_value}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
	}

	draw_area->plotLayout()->insertRow(0);
	SoapTextElement* title = new SoapTextElement(draw_area, "Observed transcriptomes", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(0, 0, title);

	std::tie(counts, edges, locs) = custom::histogram(simulate_score, 32);
	normed_counts = log10(counts.cast<double>() + 1);
	axis_rect = new QCPAxisRect(draw_area, true);
	draw_area->plotLayout()->addElement(2, 0, axis_rect);
	custom_plot::bar_plot(draw_area, axis_rect, Qt::blue, locs, normed_counts, 12, "Doublet score", "Prob. density (log10 + 1)", gs);
	if (threshold.first) {
		double threshold_value = threshold.second;
		custom_plot::patch::line(draw_area, axis_rect, QVector<double>{threshold_value, threshold_value}, QVector<double>{0, normed_counts.maxCoeff()}, Qt::black, 4, Qt::DashLine);
	}

	draw_area->plotLayout()->insertRow(2);
	SoapTextElement* title2 = new SoapTextElement(draw_area, "Simulated doublets", QFont("Arial", 20, QFont::Bold));
	draw_area->plotLayout()->addElement(2, 0, title2);

	this->draw_suite_->update(draw_area);

	G_LOG("Detect doublets by Scrublet finished.");
}

void SingleCellRnaItem::s_scrublet() {
	G_GETLOCK;

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Counts data is lost.");
			G_UNLOCK;
		return;
	}

	G_LOG("Detect doublets by Scrublet...");
	ScrubletWorker* worker = new ScrubletWorker(counts->mat_);
	G_LINK_WORKER_THREAD(ScrubletWorker, x_scrublet_ready, SingleCellRnaItem, s_receive_scrublet)
};

void SingleCellRnaItem::s_bubble_plot() {

	auto& metadata = this->data()->metadata()->mat_;
	QMap<QString, QStringList> map = metadata.get_factor_information();
	if (map.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}
	const auto& gs = this->draw_suite_->graph_settings_;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Bubble Plot",
		{ "Gene Names", "Group:Group", "Normalized" },
		{ soap::InputStyle::MultipleLineEdit, soap::InputStyle::FactorChoice, soap::InputStyle::SwitchButton},
		QList<QStringList>(), { map}
	);

	if (settings.isEmpty())return;

	auto gene_names = multiple_line_edit_to_list(settings[0]);
	if (gene_names.isEmpty())return;

	auto filter = custom::in(gene_names, this->data()->counts()->rownames_);
	if (filter.count() == 0) {
		G_NOTICE("No Gene Found in data.");
		return;
	}
	gene_names = custom::sliced(gene_names, filter);

	auto gene_index = custom::index_of(gene_names, this->data()->counts()->rownames_);
	int n_gene = gene_names.size();
	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);
	if (levels.isEmpty())return;
	int n_level = levels.size();
	QStringList factor = metadata.get_qstring(factor_name);
	bool is_normalized = switch_to_bool(settings[2]);

	auto counts = this->data()->counts();
	auto normalized = this->data()->normalized();

	if (is_normalized && (normalized == nullptr)) {
		G_WARN("Data has not been normalized.");
		return;
	}
	else if (!is_normalized && (counts == nullptr)) {
		G_WARN("Count data is missed.");
		return;
	}

	Eigen::MatrixXd values(n_level, n_gene);
	Eigen::MatrixXi proportion(n_level, n_gene);
	for (int i = 0; i < n_level; ++i) {
		filter = custom::equal(factor, levels[i]);
		for (int j = 0; j < n_gene; ++j) {
			if (is_normalized) {
				Eigen::ArrayXd exp = normalized->mat_.row(gene_index[j]);
				exp = custom::sliced(exp, filter);
				values(n_level - i - 1, j) = exp.mean();
				proportion(n_level - i - 1, j) = ceil((exp > 0).count() / (double)exp.size() * 50) + 1;
			}
			else {
				Eigen::ArrayXi exp = counts->mat_.row(gene_index[j]);
				exp = custom::sliced(exp, filter);
				values(n_level - i - 1, j) = exp.cast<double>().mean();
				proportion(n_level - i - 1, j) = ceil((exp > 0).count() / (double)exp.size() * 50) + 1;
			}
		}
	}
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
	custom_plot::bubble_plot(draw_area, axis_rect, legend_layout, custom::reversed(levels), gene_names, values, proportion, gs.get_legend_title("Expression"), gs.get_legend_title("Proportion", 1), gs);
	this->draw_suite_->update(draw_area);
};

void SingleCellRnaItem::s_distribution_plot() {

	auto& metadata = this->data()->metadata()->mat_;
	QStringList features = metadata.get_factor_name();
	if (features.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}
	const auto& gs = this->draw_suite_->graph_settings_;

	features = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Feature",
		{ "Feature", "Group", "Normalized", "Show Value" },
		{ soap::InputStyle::StringLineEdit, soap::InputStyle::ComboBox
		, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton},
		{ features}
	);
	if (features.isEmpty())return;

	bool is_normalized = switch_to_bool(features[2]);

	auto counts = this->data()->counts();
	auto normalized = this->data()->normalized();

	if (is_normalized && (normalized == nullptr)) {
		G_NOTICE("Data has not been normalized.");
		return;
	}

	bool show_value = switch_to_bool(features[3]);
	QString group_name = features[1];
	features = features[0].split(",");
	for (auto& feature : features) {
		feature = feature.trimmed();
	}

	auto& gene_names = counts->rownames_;

	QStringList valid_features;
	for (auto& feature : features) {
		if ((metadata.contains(feature) && (metadata.data_type_[feature] == CustomMatrix::DataType::IntegerNumeric ||
			metadata.data_type_[feature] == CustomMatrix::DataType::DoubleNumeric)) || gene_names.contains(feature))
			valid_features.append(feature);
	}
	if (valid_features.isEmpty()) {
		G_LOG("No valid features!");
		return;
	}
	QCustomPlot* draw_area = custom_plot::initialize_plot(gs);
	QStringList group = metadata.get_qstring(group_name);

	QStringList group_labels;
	int group_number;
	if (metadata.data_type_[group_name] == CustomMatrix::DataType::QStringFactor) {
		group_number = metadata.string_factors_[group_name].size();
		group_labels = metadata.string_factors_[group_name];
	}
	else {
		group_number = metadata.integer_factors_[group_name].size();
		group_labels = custom::cast<QString>(custom::sorted(metadata.integer_factors_[group_name]));
	}
	auto colors = gs.palette(group_labels);
	Eigen::ArrayXd feature_data;

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	for (int i = 0; i < valid_features.size(); ++i) {
		QString feature = valid_features[i];
		if (metadata.contains(feature) && (metadata.data_type_[feature] == CustomMatrix::DataType::IntegerNumeric || metadata.data_type_[feature] == CustomMatrix::DataType::DoubleNumeric)) {
			feature_data = custom::cast<Eigen::ArrayX>(metadata.get_double(feature));
		}
		else {
			if (!is_normalized) {
				feature_data = counts->get_row(feature).cast<double>();
			}
			else {
				feature_data = normalized->get_row(feature);
			}
		}
		QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);
		axis_rect->setMarginGroup(QCP::msLeft, margin_group);
		draw_area->plotLayout()->addElement(i, 0, axis_rect);

		auto [min, max] = custom_plot::patch::violin_batch(
			draw_area,
			axis_rect,
			group,
			group_labels,
			colors,
			feature_data,
			1.0,
			2.0,
			16
		);

		custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * group_number), custom_plot::utility::get_range(min, max));
		custom_plot::patch::clear_bottom_axis(axis_rect);
		custom_plot::set_left_title(axis_rect, feature, gs);
		if (!show_value) {
			axis_rect->axis(QCPAxis::atLeft)->setTicks(false);
			axis_rect->axis(QCPAxis::atLeft)->setTickLabels(false);
		}
		axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);

		if (i == valid_features.size() - 1) {
			custom_plot::set_bottom_axis_label(
				axis_rect,
				Eigen::ArrayXd::LinSpaced(group_number, 1, 2 * group_number - 1),
				group_labels,
				6,
				gs
			);
		}

	}

	this->draw_suite_->update(draw_area);
};

void SingleCellRnaItem::s_view_quality() {

	auto* metadata = &this->data()->metadata()->mat_;

	if (!(metadata->contains(METADATA_RNA_UMI_NUMBER) &&
		metadata->contains(METADATA_RNA_UNIQUE_GENE_NUMBER) &&
		metadata->contains(METADATA_RNA_MITOCHONDRIAL_CONTENT))) {
		G_WARN("Quality information has lost.");
		return;
	}
	if (!(metadata->data_type_[METADATA_RNA_UMI_NUMBER] == CustomMatrix::DataType::IntegerNumeric &&
		metadata->data_type_[METADATA_RNA_UNIQUE_GENE_NUMBER] == CustomMatrix::DataType::IntegerNumeric &&
		metadata->data_type_[METADATA_RNA_MITOCHONDRIAL_CONTENT] == CustomMatrix::DataType::DoubleNumeric)) {
		G_WARN("Quality information format has been overwritten.");
		return;
	}

	auto colors = custom_plot::utility::kmeans_palette(3);

	auto trans = metadata->get_double(METADATA_RNA_UMI_NUMBER);
	int control_point_number = 24;

	const auto& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, layout] = custom_plot::prepare_lg(gs);

	QCPAxisRect* axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 0, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[0], "", "UMI Count");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[0], "", "UMI Count");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata->get_double(METADATA_RNA_UNIQUE_GENE_NUMBER);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 1, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[1], "", "Gene Count");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[1], "", "Gene Count");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata->get_double(METADATA_RNA_MITOCHONDRIAL_CONTENT);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 2, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[2], "", "Mitochondrial Content");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[2], "", "Mitochondrial Content");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	custom_plot::add_title(draw_area, "Data Quality", gs);

	this->draw_suite_->update(draw_area);
};

void SingleCellRnaItem::__s_delete_this() {

	G_GETLOCK;
	G_UNLOCK;

	if (!YesOrNoDialog::get_response("Delete Single Cell RNA Data", "This data will be deleted.")) {
		return;
	}

	this->__remove_this();
};

void SingleCellRnaItem::s_set_random_state() {

	G_GETLOCK;
	G_UNLOCK;
	
	QStringList seed = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Random Seed",
		{ "New Random Seed:" + QString::number(this->data()->random_state_) },
		{ soap::InputStyle::StringLineEdit}
	);

	if (seed.isEmpty()) {
		return;
	}
	
	this->data()->random_state_ = seed[0].toInt();
};

void SingleCellRnaItem::s_set_species() {

	G_GETLOCK;
	G_UNLOCK;
	
	QStringList species = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Species",
		{ "New Species" },
		{ soap::InputStyle::ComboBox},
		{ { "Human", "Mouse", "Undefined" }}
	);

	if (species.isEmpty()) {
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
};

void SingleCellRnaItem::s_receive_differential_analysis(DifferentialAnalysis da, QString name) {

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

	G_LOG("Differential Expression Analysis finished");
};

void SingleCellRnaItem::s_receive_velocyto(VelocytoBase* velocyto_base) {

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

void SingleCellRnaItem::s_monocle3() {

	G_GETLOCK;

	auto embedding_names = custom::keys(DATA_SUBMODULES(Embedding));

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
		int n_cell = this->data()->counts()->cols();

		Monocle3Worker* worker = new Monocle3Worker(
			DATA_SUBMODULES(Embedding)[embedding_name],
			Eigen::ArrayX<bool>::Constant(n_cell, true)
		);

		G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellRnaItem, s_receive_monocle3);
	}
	else {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Monocle3 Settings",
			{ "Embedding" , "use part", "choose cell group"},
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
				DATA_SUBMODULES(Embedding)[embedding_name],
				custom::in(factor, levels)
			);

			G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellRnaItem, s_receive_monocle3);
		}
		else {
			int n_cell = this->data()->counts()->cols();

			Monocle3Worker* worker = new Monocle3Worker(
				DATA_SUBMODULES(Embedding)[embedding_name],
				Eigen::ArrayX<bool>::Constant(n_cell, true)
			);

			G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellRnaItem, s_receive_monocle3);
		}
	}
	
};

void SingleCellRnaItem::s_receive_monocle3(Monocle3* monocle3) {

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

void SingleCellRnaItem::s_velocyto() {

	if (this->data()->data_type_ == SingleCellRna::DataType::Integrated) {
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
	G_LINK_WORKER_THREAD(VelocytoWorker, x_velocyto_ready, SingleCellRnaItem, s_receive_velocyto);
}

void SingleCellRnaItem::s_find_deg() {

	G_GETLOCK;

	auto factor_map = this->data()->metadata()->mat_.get_factor_information(false);
	if (factor_map.isEmpty()) {
		G_LOG("No suitable feature for Find DEG");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Differential Expression Analysis settings",
		{ "Choose Group:1", "Min Percentage:0.1", "P Adjust Method" },
		{soap::InputStyle::CompareLayout, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox},
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

	QStringList metadata = this->data()->metadata()->mat_.get_qstring(factor_name);

	double minimum_percentage = settings[1].toDouble();
	if (minimum_percentage < 0.0 || minimum_percentage >= 1.0) {
		G_WARN("Invalid Percentage! Reset to 0.1.");
		minimum_percentage = 0.1;
	}

	QString p_adjust_method = settings[2];

	G_LOG("Finding DEG in " + factor_name + "...");
	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		minimum_percentage,
		p_adjust_method
	);
	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, SingleCellRnaItem, s_receive_differential_analysis);
};

void SingleCellRnaItem::s_receive_gsea(GSEA gsea) {

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

void SingleCellRnaItem::s_gsea() {

	auto normalized = this->data()->normalized();
	if (normalized == nullptr) {
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

	auto index = custom::match(metadata, comparison1);
	if (index.size() > downsample) {
		index = custom::sample(index, downsample, random_state);
	}
	selected << index;

	if (comparison2 == "REST") {
		index = custom::which(!custom::equal(metadata, comparison1));
	}
	else {
		index = custom::match(metadata, comparison2);
	}
	if (index.size() > downsample) {
		index = custom::sample(index, downsample, random_state);
	}
	selected << index;

	metadata = custom::reordered(metadata, selected);
	G_LOG("GSEAing in " + factor_name + "...");
	SparseDouble tmp = normalized->col_reordered(selected);

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
	G_LINK_WORKER_THREAD(GseaWorker, x_gsea_ready, SingleCellRnaItem, s_receive_gsea)
}

void SingleCellRnaItem::s_receive_infercnv(CNV* cnv) {

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

void SingleCellRnaItem::s_infercnv() {

	auto counts = this->data()->counts();
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
		{ "Feature:Reference", "Downsample Number", },
		{ soap::InputStyle::FactorChoice, soap::InputStyle::ComboBox },
		{ { "200", "300", "500", "1000", "ALL" } },
		{ factor_map }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	};

	auto [feature, reference] = factor_choice_to_pair(settings[0]);

	QStringList metadata = this->data()->metadata()->mat_.get_qstring(feature);
	QStringList levels = custom::unique(metadata);

	QVector<int> selected;
	if (settings[1] != "ALL") {
		int downsample = settings[1].toInt();
		for (const auto& factor : levels) {
			auto index = custom::match(metadata, factor);
			if (index.size() > downsample) {
				index = custom::sample(index, downsample, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		for (const auto& factor : levels) {
			auto index = custom::match(metadata, factor);
			selected << index;
		}
	}

	if (reference.isEmpty()) {
		reference = custom::unique(metadata);
	}
	else {
		reference = custom::unique(reference);
	}

	G_LOG("InferCnv start...");
	InferCnvWorker* worker = new InferCnvWorker(
		counts->col_reordered(selected),
		custom::reordered(metadata, selected),
		reference,
		this->data()->species_
	);

	G_LINK_WORKER_THREAD(InferCnvWorker, x_cnv_ready, SingleCellRnaItem, s_receive_infercnv);
};



void SingleCellRnaItem::s_scicnv() {
	
	auto normalized = this->data()->normalized();
	if (normalized == nullptr) {
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
		{ soap::InputStyle::FactorChoice, soap::InputStyle::ComboBox,
		soap::InputStyle::NumericLineEdit, soap::InputStyle::NumericLineEdit},
		{ { "200", "500", "1000", "All" }},
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
	QStringList factor = this->data()->metadata()->mat_.get_qstring(feature);
	QStringList levels;

	auto& metadata = this->data()->metadata()->mat_;
	if (metadata.data_type_[feature] == CustomMatrix::DataType::QStringFactor) {
		levels = metadata.string_factors_[feature];
	}
	else {
		levels = custom::cast<QString>(metadata.integer_factors_[feature]);
	}
	QVector<int> selected;
	if (settings[1] != "All") {
		int downsample = settings[1].toInt();
		for (const auto& level : levels) {
			auto index = custom::match(factor, level);
			if (index.size() > downsample) {
				index = custom::sample(index, downsample, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		for (const auto& level : levels) {
			auto index = custom::match(factor, level);
			selected << index;
		}
	}
	if (selected.isEmpty()) {
		G_WARN("No cell is selected.")
			G_UNLOCK;
		return;
	}
	G_LOG("SciCnv start...");
	ScicnvWorker* worker = new ScicnvWorker(
		normalized->col_reordered(selected),
		custom::reordered(factor, selected),
		reference,
		this->data()->species_,
		threshold,
		sharpness
	);
	G_LINK_WORKER_THREAD(ScicnvWorker, x_cnv_ready, SingleCellRnaItem, s_receive_scicnv);
};

void SingleCellRnaItem::s_receive_scicnv(CNV* cnv) {

	QString title = this->signal_emitter_->get_unique_name("SciCNV");
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

	G_LOG("SciCnv finished");
};

void SingleCellRnaItem::s_combine_existed_metadata() {

	G_GETLOCK;

	G_UNLOCK;

	QStringList metadata_names = this->data()->metadata()->mat_.colnames_;

	if (metadata_names.isEmpty()) {
		G_WARN("Metadata is Empty.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Metadata Combine Settings",
		{ "Choose Metadata", "Combine Type", "New Metadata Name", "New Metadata Type", "Join Character: " },
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::ComboBox,
		soap::InputStyle::StringLineEdit, soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ metadata_names,
		{ "Text Concatenate", "Numeric Add", "Numeric Minus", "Numeric Multiply", "Numeric Divide" },
		{ "String Factor", "String", "Numeric", "Integer", "Integer Factor" }}
	);

	if (settings.isEmpty()) {
		return;
	}

	QStringList components = simple_choice_to_list(settings[0]);

	if (components.size() < 2) {
		return;
	}

	QString combine_type = settings[1];
	QString target_data_type = settings[3];
	QString new_metadata_name = settings[2];
	QString joiner = settings[4];

	auto& metadata = this->data()->metadata()->mat_;

	if (metadata.contains(new_metadata_name)) {
		if (!YesOrNoDialog::get_response("Warning: Name Duplicated!", "Overwrite existing metadata?")) {
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
			this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<int>(res), CustomMatrix::DataType::IntegerNumeric);
		}
		else if (target_data_type == "Numeric") {
			this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<double>(res), CustomMatrix::DataType::DoubleNumeric);
		}
		else if (target_data_type == "Integer Factor") {
			this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<int>(res), CustomMatrix::DataType::IntegerFactor);
		}
		else if (target_data_type == "String Factor") {
			this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::QStringFactor);
		}
		else if (target_data_type == "String") {
			this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::QString);
		}
	}
	else {
		auto data_types = custom::sapply(components, [&metadata](auto&& name) {return metadata.data_type_[name]; });

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
			else if (combine_type == "Numeric Minus") {
				if (n_choosed != 2) {
					G_WARN("Numeric Minus only support two metadata operation.");
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
					return;
				}
				for (int i = 0; i < nrow; ++i) {
					if (choosed[1][i] == 0.0) {
						G_WARN("Meet 0.0 in divisor, process terminate");
						return;
					}
					res[i] /= choosed[1][i];
				}
			}

			if (target_data_type == "Integer") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<int>(res), CustomMatrix::DataType::IntegerNumeric);
			}
			else if (target_data_type == "Numeric") {
				this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::DoubleNumeric);
			}
			else if (target_data_type == "Integer Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<int>(res), CustomMatrix::DataType::IntegerFactor);
			}
			else if (target_data_type == "String Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<QString>(res), CustomMatrix::DataType::QStringFactor);
			}
			else if (target_data_type == "String") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<QString>(res), CustomMatrix::DataType::QString);
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
					return;
				}
				for (int i = 0; i < nrow; ++i) {
					if (choosed[1][i] == 0) {
						G_WARN("Meet 0 in divisor, process terminate");
						return;
					}
					res[i] /= choosed[1][i];
				}
			}

			if (target_data_type == "Integer") {
				this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::IntegerNumeric);
			}
			else if (target_data_type == "Numeric") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<double>(res), CustomMatrix::DataType::DoubleNumeric);
			}
			else if (target_data_type == "Integer Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, res, CustomMatrix::DataType::IntegerFactor);
			}
			else if (target_data_type == "String Factor") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<QString>(res), CustomMatrix::DataType::QStringFactor);
			}
			else if (target_data_type == "String") {
				this->data()->metadata()->mat_.update(new_metadata_name, custom::cast<QString>(res), CustomMatrix::DataType::QString);
			}
		}
	}

	this->signal_emitter_->x_update_interface();
};

void SingleCellRnaItem::s_edit_metadata() {

	G_GETLOCK;

	G_UNLOCK;
	
	LogicHandler lh(this->data());

	auto&& metadata = this->data()->metadata()->mat_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Filter Settings",
		{ "Select Cell", "Feature to edit", "New Value" },
		{soap::InputStyle::LogicLayout, soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit},
		{metadata.colnames_},
		{},
		{&lh}
	);

	if (settings.isEmpty()) {
		return;
	}

	auto filter = lh.resolve(settings[0]);
	QString feature_name = settings[1];
	QString new_val = settings[2];

	if (filter.count() == 0) {
		G_WARN("No Cell Meets Requirements.");
		return;
	}

	metadata.edit(feature_name, filter, new_val);

};

void SingleCellRnaItem::s_add_metadata() {

	G_GETLOCK;

	G_UNLOCK;

	QStringList metadata_information = CommonDialog::get_response(
		this->signal_emitter_,
		"New Metadata Setting",
		{ "Name", "Type", "Initial Value", "Initiate from existing metadata:no", "Metadata from" },
		{ soap::InputStyle::StringLineEdit, soap::InputStyle::ComboBox,
		soap::InputStyle::StringLineEdit, soap::InputStyle::SwitchButton, soap::InputStyle::ComboBox},
		{ { "String Factor", "String", "Numeric", "Integer", "Integer Factor" }, this->data()->metadata()->mat_.colnames_}
	);

	if (metadata_information.isEmpty()) {
		return;
	}

	auto& metadata = this->data()->metadata()->mat_;
	if (metadata.contains(metadata_information[0])) {
		if (!YesOrNoDialog::get_response("Warning : Name Duplicated!", "Overwrite existing metadata?")) {
			return;
		}
	}

	bool existed = switch_to_bool(metadata_information[3]);

	if (existed) {
		QString existed_name = metadata_information[4];
		CustomMatrix::DataType datatype = metadata.data_type_[existed_name];
		if (datatype == CustomMatrix::DataType::DoubleNumeric) {
			metadata.update(metadata_information[0], metadata.get_double(existed_name), CustomMatrix::DataType::DoubleNumeric);
		}
		else if (datatype == CustomMatrix::DataType::IntegerNumeric) {
			metadata.update(metadata_information[0], metadata.get_integer(existed_name), CustomMatrix::DataType::IntegerNumeric);
		}
		else if (datatype == CustomMatrix::DataType::QStringFactor) {
			metadata.update(metadata_information[0], metadata.get_qstring(existed_name), CustomMatrix::DataType::QStringFactor);
		}
		else if (datatype == CustomMatrix::DataType::QString) {
			metadata.update(metadata_information[0], metadata.get_qstring(existed_name), CustomMatrix::DataType::QString);
		}
		else if (datatype == CustomMatrix::DataType::IntegerFactor) {
			metadata.update(metadata_information[0], metadata.get_integer(existed_name), CustomMatrix::DataType::IntegerFactor);
		}
	}
	else {
		QString datatype = metadata_information[1];
		int size = metadata.rownames_.size();
		if (datatype == "Integer") {
			metadata.update(metadata_information[0], QVector<int>(size, metadata_information[2].toInt()), CustomMatrix::DataType::IntegerNumeric);
		}
		else if (datatype == "Numeric") {
			metadata.update(metadata_information[0], QVector<double>(size, metadata_information[2].toDouble()), CustomMatrix::DataType::DoubleNumeric);
		}
		else if (datatype == "Integer Factor") {
			metadata.update(metadata_information[0], QVector<int>(size, metadata_information[2].toInt()), CustomMatrix::DataType::IntegerFactor);
		}
		else if (datatype == "String Factor") {
			metadata.update(metadata_information[0], QStringList(size, metadata_information[2]), CustomMatrix::DataType::QStringFactor);
		}
		else if (datatype == "String") {
			metadata.update(metadata_information[0], QStringList(size, metadata_information[2]), CustomMatrix::DataType::QString);
		}
	}

	this->signal_emitter_->x_update_interface();
};

void SingleCellRnaItem::s_receive_integrated_data(SingleCellRna* data, QList<const SingleCellRna* > items) {

	this->signal_emitter_->unlock(custom::sapply(items, [this](auto* data) {return this->signal_emitter_->search((void*)data); }));

	this->signal_emitter_->x_data_create_soon(data, soap::VariableType::SingleCellRna, "Integrated scRNA");
};

void SingleCellRnaItem::s_integrate() {

	auto& variables = this->signal_emitter_->variable_information_;

	QStringList available_data = this->signal_emitter_->get_type_variable(soap::VariableType::SingleCellRna).keys();
	if (available_data.size() < 2) {
		G_LOG("No enough single cell data found.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Integrate SingleCellRna",
		{ "Data", "Style" },
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::ComboBox},
		{ available_data, { "Distinguish barcodes from different data", "Merge the same barcode" }}
	);
	if (settings.isEmpty())return;

	available_data = simple_choice_to_list(settings[0]);
	if (available_data.size() < 2)return;
	if (!custom::is_unique(available_data)) {
		G_WARN("Can not integrate the same data!");
		return;
	}

	auto its = this->signal_emitter_->search(available_data);

	if (!this->signal_emitter_->try_lock(its)) {
		G_WARN("Please waiting for the computation in progress.");
		return;
	}

	auto choosed = custom::sapply(its, [](IndexTree* it) {return static_cast<const SingleCellRna*>(it->data_); });

	auto species = custom::unique(custom::sapply(choosed, [](const SingleCellRna* single_cell_rna) {return single_cell_rna->species_; }));

	if (species.size() != 1) {
		G_LOG("Unmatched species!");
		return;
	}

	bool distinguish = (settings[1] == "Distinguish barcodes from different data");

	if (!distinguish) {
		G_NOTICE("Metadata will not be inherited to avoid collision.");
	}

	G_LOG("Start integrating data...");

	IntegrateWorker* worker = new IntegrateWorker(choosed, species[0], distinguish);

	G_LINK_WORKER_THREAD(IntegrateWorker, x_scrna_ready, SingleCellRnaItem, s_receive_integrated_data)
};
