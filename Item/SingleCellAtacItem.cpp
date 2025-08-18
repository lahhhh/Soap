#include "SingleCellAtacItem.h"

#include "CommonDialog.h"
#include "YesOrNoDialog.h"
#include "StatisticsDialog.h"

#include "ItemIOWorker.h"
#include "LogNormalizeWorker.h"
#include "TfidfWorker.h"
#include "SvdWorker.h"
#include "UmapWorker.h"
#include "TsneWorker.h"
#include "LeidenPartitionWorker.h"
#include "SlmWorker.h"
#include "DifferentialAnalysisWorker.h"
#include "FragmentsQualityViewWorker.h"
#include "TssPlotWorker.h"
#include "CreateCoverageTrackWorker.h"
#include "MacsCallPeakWorker.h"
#include "CalculateGeneActivityWorker.h"
#include "CellTypeAnnotationWorker.h"
#include "MotifLocateWorker.h"
#include "IntegrateWorker.h"
#include "LoadFragmentsWorker.h"
#include "Monocle3Worker.h"

void SingleCellAtacItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Metadata));
	this->check_variable(DATA_SUBMODULES(SparseInt));
	this->check_variable(DATA_SUBMODULES(SparseDouble));
	this->check_variable(DATA_SUBMODULES(Embedding));
	this->check_variable(DATA_SUBMODULES(DifferentialAnalysis));
	this->check_variable(DATA_SUBMODULES(GenomicRange));
	this->check_variable(DATA_SUBMODULES(MotifPosition));
	this->check_variable(DATA_SUBMODULES(CoverageTrack));
	this->check_variable(DATA_SUBMODULES(Fragments));
	this->check_variable(DATA_SUBMODULES(Cicero));

	auto item = new NoteItem(&this->data()->string_information_["Note"], this->data(), this->signal_emitter_);
	this->addChild(item);
};

void SingleCellAtacItem::__set_menu() {
	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	//Normalize
	ADD_MAIN_MENU("Normalize");

	ADD_MENU("Normalize | TFIDF", "TFIDF", "Normalize");

	ADD_ACTION("Default", "Normalize | TFIDF", s_tfidf_default);

	ADD_ACTION("Normalize Gene Activity", "Normalize", s_normalize_gene_activity);

	// Dimensional Reduction
	ADD_MAIN_MENU("Dimension Reduction");

	ADD_MENU("Dimension Reduction | PCA", "PCA", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | PCA", s_svd_default);
	ADD_ACTION("Custom", "Dimension Reduction | PCA", s_svd_custom);

	ADD_MENU("Dimension Reduction | UMAP", "UMAP", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | UMAP", s_umap_default);
	ADD_ACTION("Custom", "Dimension Reduction | UMAP", s_umap_custom);

	ADD_MENU("Dimension Reduction | tSNE", "tSNE", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | tSNE", s_tsne_default);
	ADD_ACTION("Custom", "Dimension Reduction | tSNE", s_tsne_custom);


	//cluster
	ADD_MAIN_MENU("Cluster");

	ADD_MENU("Cluster | Louvain", "Louvain", "Cluster");

	ADD_MENU("Cluster | Louvain | Louvain", "Louvain", "Cluster | Louvain")

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

	// DAP
	ADD_MAIN_MENU("Expression Analysis");

	ADD_ACTION("Find DAP", "Expression Analysis", s_find_dap);

	// calculate gene activity
	ADD_MAIN_ACTION("Calculate Gene Activity", s_calculate_gene_activity);

	//fragments
	ADD_MAIN_MENU("Fragments");

	ADD_ACTION("Load Fragments", "Fragments", s_load_fragments);

	// find motif

	ADD_MAIN_ACTION("Find Motifs", s_find_motifs);

	//coverage

	ADD_MAIN_MENU("Coverage");

	ADD_ACTION("Coverage Plot", "Coverage", s_coverage_plot);
	ADD_ACTION("Create Coverage Track", "Coverage", s_create_coverage_track);

	// call peak
	ADD_MAIN_MENU("Call Peaks");

	ADD_ACTION("MACS", "Call Peaks", s_call_peaks_by_macs);

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

	//filter
	ADD_MAIN_MENU("Filter");

	ADD_ACTION("By Features", "Filter", s_filter_by_features);
	ADD_ACTION("By Parameters", "Filter", s_filter_by_parameters);
	ADD_ACTION("By Quality Metrics", "Filter", s_filter_by_quality_metrics);


	// landscape plot
	ADD_MAIN_ACTION("Landscape Plot", s_show_atac_landscape);

	// integrate
	ADD_MAIN_ACTION("Integrate", s_integrate);

	// duplicate
	ADD_MAIN_ACTION("Duplicate", __s_duplicate);

	// save
	ADD_MAIN_ACTION("Save", __s_export_as_item);

	// delete
	ADD_MAIN_ACTION("Delete", __s_delete_this);
}


void SingleCellAtacItem::col_slice(const Eigen::ArrayX<bool>& col_slice) {

	this->__clear_reserve_data();

	this->data()->col_slice(col_slice);

	this->signal_emitter_->update_information(this->title_, soap::VariableType::SingleCellAtac, this->data(), true);

	this->__check_data();

	G_LOG("Filter finished.");
};

void SingleCellAtacItem::__s_rename() {
	G_GETLOCK;

	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Rename SingleCellAtac",
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

void SingleCellAtacItem::s_statistics() {

	StatisticsDialog::get_response(this->data());
};

void SingleCellAtacItem::__s_delete_this() {

	G_GETLOCK;
	G_UNLOCK;

	if (!YesOrNoDialog::get_response("Delete Single Cell ATAC Data", "This data will be deleted.")) {
		return;
	}

	this->__remove_this();
};

void SingleCellAtacItem::s_normalize_gene_activity() {

	G_UNLOCK;

	auto gene_activity_counts = this->data()->gene_activity_counts();

	if (gene_activity_counts == nullptr) {
		G_WARN("Counts Data Missed.");
		G_UNLOCK;
		return;
	}

	if (auto normalized = this->gene_activity_normalized()) {
		normalized->__remove_this();
	}

	G_LOG("Gene Activity Normalization Start...");

	auto* worker = new LogNormalizeWorker(gene_activity_counts);
	G_LINK_WORKER_THREAD(LogNormalizeWorker, x_log_normalize_ready, SingleCellAtacItem, s_receive_normalized_gene_activity);
};

void SingleCellAtacItem::s_receive_normalized_gene_activity(SparseDouble* data) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_GENE_ACTIVITY_NORMALIZED);

	data->data_type_ = SparseDouble::DataType::GeneActivity;

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
};

void SingleCellAtacItem::s_tfidf_default() {

	G_UNLOCK;

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Counts data missed.");
		G_UNLOCK;
		return;
	}

	if (auto normalized = this->normalized()) {
		normalized->__remove_this();
	}

	G_LOG("TFIDF start...");
	TfidfWorker* worker = new TfidfWorker(
		counts,
		10000
	);

	G_LINK_WORKER_THREAD(TfidfWorker, x_tfidf_ready, SingleCellAtacItem, s_receive_tfidf);
};

void SingleCellAtacItem::s_receive_tfidf(SparseDouble* data) {

	QString title = VARIABLE_NORMALIZED;

	title = this->signal_emitter_->get_unique_name(title);

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
}


void SingleCellAtacItem::s_svd_custom() {

	G_GETLOCK;

	auto normalized = this->data()->normalized();
	if (normalized == nullptr) {
		G_WARN("Normalized data missed!");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"SVD Settings",
		{ "Feature Percentile:25" },
		{soap::InputStyle::IntegerLineEdit}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	int perc = settings[0].toInt();

	if (perc < 1 || perc > 100) {
		G_WARN("Illegal Parameter.");
		G_UNLOCK;
		return;
	}


	if (auto item = this->pca()) {
		item->__remove_this();
	}

	G_LOG("SVD start...");

	SvdWorker* worker = new SvdWorker(normalized->mat_, perc, 50, this->data()->random_state_);
	G_LINK_WORKER_THREAD(SvdWorker, x_svd_ready, SingleCellAtacItem, s_receive_svd);
};

void SingleCellAtacItem::s_svd_default() {

	G_GETLOCK;

	auto normalized = this->data()->normalized();
	if (normalized == nullptr) {
		G_WARN("Normalized data missed!");
		G_UNLOCK;
		return;
	}

	if (auto item = this->pca()) {
		item->__remove_this();
	}

	G_LOG("SVD start...");

	SvdWorker* worker = new SvdWorker(normalized->mat_, 100, 50, this->data()->random_state_);
	G_LINK_WORKER_THREAD(SvdWorker, x_svd_ready, SingleCellAtacItem, s_receive_svd);
}

void SingleCellAtacItem::s_receive_svd(Eigen::MatrixXd mat, QVector<double> sd) {

	QString pca_name, sd_name;

	pca_name = VARIABLE_PCA;
	sd_name = VARIABLE_ATAC_PCA_STANDARD_DEVIATION;

	QString title = this->signal_emitter_->get_unique_name(pca_name);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Pca,
		mat,
		this->data()->counts()->colnames_,
		custom::paste(pca_name, custom::cast<QString>(custom::seq_n(1, mat.cols())), "-")
	);

	this->data()->double_vectors_[sd_name] = sd;

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(Embedding)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);
	this->set_item(item);
};


void SingleCellAtacItem::s_umap_default() {

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

	Eigen::MatrixXd umap_input;

	int start_dimension = 0;
	int last_dimension = pca->data_.mat_.cols() - 1;
	if (last_dimension > 19) {
		last_dimension = 19;
	}

	const int random_state = this->data()->random_state_;
	G_LOG("Default ATAC UMAP start, using 1~" + QString::number(last_dimension + 1) + " PCA components...");

	G_NOTICE("If you use SVD before, it is recommended to exclude the first PCA dimension in UMAP because of its high correlation with sequencing depth.");

	umap_input = pca->data_.mat_.block(0, start_dimension, pca->data_.mat_.rows(), last_dimension - start_dimension + 1);

	const QString
		metric = "Angular",
		init = "Random";
	constexpr double
		learning_rate = 1.0,
		minimum_distance = 0.3,
		spread = 1.0,
		set_op_mix_ratio = 1.0,
		repulsion_strength = 1.0;

	constexpr int
		n_neighbors = 30,
		negative_sample_rate = 5,
		n_trees = 50;

	UmapWorker* worker = new UmapWorker(
		umap_input,
		n_neighbors,
		metric,
		learning_rate,
		init,
		minimum_distance,
		spread,
		set_op_mix_ratio,
		repulsion_strength,
		negative_sample_rate,
		random_state,
		n_trees
	);
	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, SingleCellAtacItem, s_receive_umap)
};

void SingleCellAtacItem::s_umap_custom() {

	G_GETLOCK;

	Eigen::MatrixXd umap_input;
	QStringList umap_from_list;

	auto pca = this->data()->pca();
	auto harmony = this->data()->harmony();

	if (pca != nullptr) {
		umap_from_list << "PCA";
	}
	if (harmony != nullptr) {
		umap_from_list << "Harmony";
	}
	if (umap_from_list.isEmpty()) {
		G_WARN("No PCA or Harmony embedding for UMAP.");
		G_UNLOCK;
		return;
	}

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
		{ { "Angular", "Euclidean", "Manhattan" }, { "Random" }, umap_from_list}
		);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	Eigen::MatrixXd* mat;
	if (settings[12] == "PCA") {
		mat = &pca->data_.mat_;
	}
	else {
		mat = &harmony->data_.mat_;
	}

	int nrow = mat->rows(), ncol = mat->cols();
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
		G_WARN("Illegal neighbor value : " + QString::number(n_neighbors) + ". Please reset.");
		G_UNLOCK;
		return;
	}
	QString metric = settings[3];
	double learning_rate = settings[4].toDouble();
	if (learning_rate <= 0) {
		G_WARN("Learning rate must be positive. Reset to 1.0");
		learning_rate = 1.0;
	}
	QString init = settings[5];
	double minimum_distance = settings[6].toDouble();
	double spread = settings[7].toDouble();
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
	double repulsion_strength = settings[9].toDouble();
	if (repulsion_strength < 0) {
		G_LOG("Repulsion strength cannot be negative. Reset to 1.0");
		repulsion_strength = 1.0;
	}
	int negative_sample_rate = settings[10].toInt();
	if (negative_sample_rate <= 0) {
		G_LOG("Negative sample rate must be positive. Reset to 5.");
		negative_sample_rate = 5;
	}
	else if (negative_sample_rate > 10) {
		G_LOG("Larger negative sample rate may cost more time.");
	}
	int random_state = settings[11].toInt();
	umap_input = mat->block(0, start_dimension - 1, nrow, last_dimension - start_dimension + 1);

	if (auto item = this->umap()) {
		item->__remove_this();
	}

	UmapWorker* worker = new UmapWorker(
		umap_input,
		n_neighbors,
		metric,
		learning_rate,
		init,
		minimum_distance,
		spread,
		set_op_mix_ratio,
		repulsion_strength,
		negative_sample_rate,
		random_state,
		50
	);
	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, SingleCellAtacItem, s_receive_umap)
};

void SingleCellAtacItem::s_receive_umap(Eigen::MatrixXd mat) {

	QString umap_name;
	umap_name = VARIABLE_UMAP;

	QString title = this->signal_emitter_->get_unique_name(umap_name);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Umap,
		mat,
		this->data()->counts()->colnames_,
		custom::paste(umap_name, custom::cast<QString>(custom::seq_n(1, mat.cols())), "-")
	);

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		this->data()->umap(),
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);
	this->set_item(item);
};


void SingleCellAtacItem::s_tsne_default() {

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

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, SingleCellAtacItem, s_receive_tsne)
};

void SingleCellAtacItem::s_tsne_custom() {

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
			soap::InputStyle::IntegerLineEdit,
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

	if (auto item = this->tsne()) {
		item->__remove_this();
	}

	unsigned int random_state = settings[3].toUInt();

	Eigen::MatrixXd tsne_input = mat->block(0, start_dimension - 1, nrow, last_dimension - start_dimension + 1);

	G_LOG("Customized tSNE start using " + from +
		" (" + QString::number(start_dimension) + "~" + QString::number(last_dimension) + ") dimensions.");


	TsneWorker* worker = new TsneWorker(tsne_input, random_state);

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, SingleCellAtacItem, s_receive_tsne)
};

void SingleCellAtacItem::s_receive_tsne(Eigen::MatrixXd mat) {

	QString tsne_name;
	tsne_name = VARIABLE_TSNE;

	QString title = this->signal_emitter_->get_unique_name(tsne_name);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Tsne,
		mat,
		this->data()->counts()->colnames_,
		custom::paste(tsne_name, custom::cast<QString>(custom::seq_n(1, mat.cols())), "-")
	);

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		this->data()->tsne(),
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);
	this->set_item(item);
};


void SingleCellAtacItem::s_leiden_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	LeidenPartitionWorker* worker = new LeidenPartitionWorker(pca->data_.mat_, "Modularity", "Euclidean", 30, 50, 0.6);
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, SingleCellAtacItem, s_receive_leiden);
}

void SingleCellAtacItem::s_louvain_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "Louvain", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellAtacItem, s_receive_louvain);
}

void SingleCellAtacItem::s_modified_louvain_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "Modified Louvain", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellAtacItem, s_receive_modified_louvain);
}

void SingleCellAtacItem::s_smart_local_moving_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "SLM", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellAtacItem, s_receive_slm);
}

void SingleCellAtacItem::s_smart_local_moving_custom() {

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
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellAtacItem, s_receive_slm);
}

void SingleCellAtacItem::s_modified_louvain_custom() {

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
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellAtacItem, s_receive_modified_louvain);
}

void SingleCellAtacItem::s_leiden_custom() {

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
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, SingleCellAtacItem, s_receive_leiden);
}

void SingleCellAtacItem::s_louvain_custom() {

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
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, SingleCellAtacItem, s_receive_louvain);
}

void SingleCellAtacItem::s_receive_leiden(QVector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_ATAC_LEIDEN_CLUSTER, cluster, CustomMatrix::DataType::IntegerFactor);
};

void SingleCellAtacItem::s_receive_louvain(std::vector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_ATAC_LOUVAIN_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void SingleCellAtacItem::s_receive_modified_louvain(std::vector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_ATAC_MODIFIED_LOUVAIN_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void SingleCellAtacItem::s_receive_slm(std::vector<int> cluster) {


	this->data()->metadata()->mat_.update(METADATA_ATAC_SLM_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};


void SingleCellAtacItem::s_find_dap() {

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

	G_LOG("Finding DAP in " + factor_name + "...");
	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		minimum_percentage,
		p_adjust_method
	);
	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, SingleCellAtacItem, s_receive_differential_analysis);
};

void SingleCellAtacItem::s_receive_differential_analysis(DifferentialAnalysis da, QString name) {

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
};


void SingleCellAtacItem::s_filter_by_quality_metrics() {

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

	Eigen::ArrayX<bool> filter = custom::greater_equal(metadata.get_const_double_reference(METADATA_TSS_ENRICHMENT), min_tss);
	filter *= custom::less_equal(metadata.get_const_double_reference(METADATA_NUCLEOSOME_SIGNAL), max_ns);
	filter *= custom::less_equal(metadata.get_const_double_reference(METADATA_BLACKLIST_RATIO), max_black);
	filter *= custom::greater_equal(metadata.get_const_double_reference(METADATA_FRIP_SCORE), min_frip);

	if (filter.count() == 0) {
		G_WARN("No cell remain after filter.");
		G_UNLOCK;
		return;
	}
	else {

		if (in_place) {
			G_LOG("Filter Object by Quality Metrics...");
			this->col_slice(filter);
		}
		else {
			G_LOG("Filter Object and Create New Object by Quality Metrics...");
			this->signal_emitter_->x_data_create_soon(this->data()->col_sliced(filter), soap::VariableType::SingleCellAtac, "Sliced SingleCellAtac");
		}
		G_UNLOCK;
	}
};

void SingleCellAtacItem::s_filter_by_features() {
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
		this->signal_emitter_->x_data_create_soon(this->data()->col_sliced(filter), soap::VariableType::SingleCellAtac, "Sliced SingleCellAtac");
	}

	G_UNLOCK;
}

void SingleCellAtacItem::s_filter_by_parameters() {
	G_GETLOCK;

	if (this->data()->counts() == nullptr) {
		G_WARN("ATAC Count data is missed.");
		return;
	}

	auto& metadata = this->data()->metadata()->mat_;
	if (!metadata.contains(METADATA_ATAC_UMI_NUMBER) ||
		metadata.data_type_[METADATA_ATAC_UMI_NUMBER] != CustomMatrix::DataType::IntegerNumeric) {

		G_NOTICE("Data of UMI number sums is not available, this setting of minimum/maximum count will be not applied.")

	}

	if (!metadata.contains(METADATA_ATAC_UNIQUE_PEAK_NUMBER) ||
		metadata.data_type_[METADATA_ATAC_UNIQUE_PEAK_NUMBER] != CustomMatrix::DataType::IntegerNumeric) {
		G_NOTICE("Data of Peak number sums is not available, this setting" 
			"of minimum/maximum unique gene number will be not applied.")
	}

	QStringList parameters = CommonDialog::get_response(
		this->signal_emitter_,
		"Set filter parameters",
		{
		"Max Count:60000",
		"Max Unique Peak:30000",
		"Min Count:3000",
		"Min Unique Peak:2000",
		"In place:no" },
		{
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

	auto filter = this->get_parameter_filter(parameters);
	if (filter.count() == 0) {
		G_WARN("No Data remained after filtering!");
		G_UNLOCK;
		return;
	}
	bool in_place = switch_to_bool(parameters[4]);

	if (in_place) {
		G_LOG("Filter Object by Parameters...");
		this->col_slice(filter);
	}
	else {
		G_LOG("Filter Object and Create New Object by Parameters...");
		this->signal_emitter_->x_data_create_soon(this->data()->col_sliced(filter), soap::VariableType::SingleCellAtac, "Sliced SingleCellAtac");
	}

	G_UNLOCK;
}


Eigen::ArrayX<bool> SingleCellAtacItem::get_parameter_filter(const QStringList& parameters) {

	const int ncol = this->data()->counts()->mat_.cols();
	Eigen::ArrayX<bool> passed_column = Eigen::ArrayX<bool>::Constant(ncol, true);

	auto& metadata = this->data()->metadata()->mat_;

	if (!parameters[0].isEmpty() || !parameters[2].isEmpty()) {
		if (metadata.contains(METADATA_ATAC_UMI_NUMBER) &&
			metadata.data_type_[METADATA_ATAC_UMI_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {

			auto& n_umi = metadata.get_const_integer_reference(METADATA_ATAC_UMI_NUMBER);
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
		if (metadata.contains(METADATA_ATAC_UNIQUE_PEAK_NUMBER) &&
			metadata.data_type_[METADATA_ATAC_UNIQUE_PEAK_NUMBER] == CustomMatrix::DataType::IntegerNumeric) {

			auto& n_gene = metadata.get_const_integer_reference(METADATA_ATAC_UNIQUE_PEAK_NUMBER);

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

	return passed_column;
}

void SingleCellAtacItem::s_receive_qc(QMap<QString, QList<double>> qc) {

	auto&& metadata = this->data()->metadata()->mat_;

	metadata.update(METADATA_TSS_ENRICHMENT, qc[METADATA_TSS_ENRICHMENT]);
	metadata.update(METADATA_TSS_PERCENTILE, qc[METADATA_TSS_PERCENTILE]);
	metadata.update(METADATA_NUCLEOSOME_SIGNAL, qc[METADATA_NUCLEOSOME_SIGNAL]);
	metadata.update(METADATA_BLACKLIST_RATIO, qc[METADATA_BLACKLIST_RATIO]);
	metadata.update(METADATA_FRIP_SCORE, qc[METADATA_FRIP_SCORE]);

	this->data()->double_vectors_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION] = qc[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION];
};

void SingleCellAtacItem::s_calculate_quality_metrics() {
	G_GETLOCK;

	if (this->data()->counts() == nullptr) {
		G_WARN("No ATAC Counts Data Found.");
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
	G_LINK_WORKER_THREAD(FragmentsQualityViewWorker, x_qc_ready, SingleCellAtacItem, s_receive_qc);
};

void SingleCellAtacItem::s_show_fragments_length_distribution() {

	if (!this->data()->double_vectors_.contains(VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION) ||
		this->data()->double_vectors_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION].size() != 1000) {
		G_WARN("Fragments Length Distribution Data has lost. Please run Quality Control again.");
		return;
	}

	const auto& distribution = this->data()->double_vectors_[VECTOR_FRAGMENTS_LENGTH_DISTRIBUTION];

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect] = custom_plot::prepare_ar(gs);

	double max_length = distribution.size() - 1;
	double minimum_x = 0, maximum_x = max_length + 1, maximum_y = std::ranges::max(distribution) * 1.1;

	draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
	QPen pen;
	pen.setWidth(1);
	pen.setColor(custom::random_color());

	draw_area->graph()->setPen(pen);
	draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
	draw_area->graph()->setData(custom::linspaced(distribution.size(), 0, max_length), distribution);

	axis_rect->axis(QCPAxis::atLeft)->setRange(0, maximum_y);
	axis_rect->axis(QCPAxis::atBottom)->setRange(minimum_x, maximum_x);

	custom_plot::set_simple_axis(axis_rect, "Fragment Length(bp)", "Proportion", gs);

	custom_plot::add_title(draw_area, "Fragments Length Distribution", gs);

	this->draw_suite_->update(draw_area);
};

void SingleCellAtacItem::s_show_quality_matrics() {

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
	auto [draw_area, layout] = custom_plot::prepare_lg(gs);

	auto colors = gs.palette({ METADATA_TSS_ENRICHMENT , METADATA_NUCLEOSOME_SIGNAL ,
		METADATA_BLACKLIST_RATIO , METADATA_FRIP_SCORE });
	const int control_point_number = 24;

	auto trans = metadata.get_double(METADATA_TSS_ENRICHMENT);
	auto axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 0, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[0], "", "TSS Score");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[0], "", "TSS Score");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_NUCLEOSOME_SIGNAL);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 1, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[1], "", "Nucleosome Signal");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[1], "", "Nucleosome Signal");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_BLACKLIST_RATIO);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 2, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[2], "", "Blacklist Percentage");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[2], "", "Blacklist Percentage");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_FRIP_SCORE);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 3, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[3], "", "In-Peak Percentage");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[3], "", "In-Peak Percentage");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	custom_plot::add_title(draw_area, "Quality Metrics", gs);

	this->draw_suite_->update(draw_area);

};

void SingleCellAtacItem::s_tss_plot() {

	if (this->data()->counts() == nullptr) {
		G_WARN("No ATAC Counts Data Found.");
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
	G_LINK_WORKER_THREAD(TssPlotWorker, x_tss_ready, SingleCellAtacItem, s_receive_tss_plot_data);
};

void SingleCellAtacItem::s_receive_tss_plot_data(Eigen::ArrayXXd tss_matrix, QVector<double> tss_vector) {
	const auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, main_layout] = custom_plot::prepare_lg(gs);

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

	custom_plot::patch::line(draw_area, tss_line_rect, custom::linspaced(n_tss, 0, n_tss - 1), tss_vector, Qt::red, 3);
	custom_plot::patch::set_range(tss_line_rect, QCPRange(0, n_tss - 1), QCPRange(0, 1.1 * std::ranges::max(tss_vector)));

	Eigen::ArrayXd bottom_axis_location(3);
	bottom_axis_location[0] = 0;
	bottom_axis_location[1] = (n_tss - 1) / 2.0;
	bottom_axis_location[2] = n_tss - 1;

	custom_plot::set_bottom_axis_label(
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
	custom_plot::patch::remove_left_bottom_axis(color_rect);
	custom_plot::patch::set_range(color_rect, QCPRange(-0.5, n_col - 0.5), QCPRange(-0.5, n_row - 0.5));

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
	gradient.setColorStopAt(0.0, QColor(custom_plot::color::firebrick3));
	gradient.setColorStopAt(0.5, QColor(custom_plot::color::gold));
	gradient.setColorStopAt(1.0, QColor(custom_plot::color::navy));
	heatmap->setGradient(gradient);
	heatmap->rescaleDataRange();

	color_scale->setMarginGroup(QCP::msTop | QCP::msBottom, horizontal_margin);

	custom_plot::add_title(draw_area, "TSS Score", gs);

	this->draw_suite_->update(draw_area);
};


void SingleCellAtacItem::s_coverage_plot() {

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
		{ "Group", "Region", "Draw Gene:yes", "Draw Link:no", "Link Cutoff:0.5", "Draw Legend:no" },
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
	G_LINK_WORKER_THREAD(CoveragePlotWorker, x_plot_ready, SingleCellAtacItem, s_receive_coverage_plot_data);
};

void SingleCellAtacItem::s_receive_coverage_plot_data(COVERAGE_PLOT_ELEMENTS res) {

	auto& gs = this->draw_suite_->graph_settings_;

	auto draw_area = custom_plot::coverage_plot(res, gs);

	this->draw_suite_->update(draw_area);
};


void SingleCellAtacItem::s_receive_coverage_track(CoverageTrack* d) {

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

void SingleCellAtacItem::s_create_coverage_track() {

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
		fragments, 
		settings[0]);
	G_LINK_WORKER_THREAD(CreateCoverageTrackWorker, x_coverage_track_ready, SingleCellAtacItem, s_receive_coverage_track);
};

void SingleCellAtacItem::s_sample() {

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
		levels = custom::cast<QString>(this->data()->metadata()->mat_.integer_factors_[feature]);
	}
	QVector<int> selected;
	if (downsample != "ALL") {
		int downsample_number = downsample.toInt();
		for (const auto& factor : levels) {
			auto index = custom::match(metadata, factor);
			if (index.size() > downsample_number) {
				index = custom::sample(index, downsample_number, this->data()->random_state_);
			}
			selected << index;
		}
	}
	else {
		selected = custom::seq_n(0, metadata.size());
	}

	if (selected.count() == 0) {
		G_UNLOCK;
		G_NOTICE("No cell is selected.")
			return;
	}

	this->signal_emitter_->x_data_create_soon(this->data()->col_reordered(selected), soap::VariableType::SingleCellAtac, "Sampled SingleCellAtac");

	G_UNLOCK;
};


void SingleCellAtacItem::s_set_random_state() {

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

void SingleCellAtacItem::s_set_species() {

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


void SingleCellAtacItem::s_edit_metadata() {

	G_GETLOCK;

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

void SingleCellAtacItem::s_combine_existed_metadata() {
	G_GETLOCK;

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
	G_UNLOCK;
};

void SingleCellAtacItem::s_add_new_metadata() {
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

void SingleCellAtacItem::s_show_atac_landscape() {

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
		"Landscape Plot Settings",
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

	G_LINK_WORKER_THREAD(AtacLandscapePlotWorker, x_plot_ready, SingleCellAtacItem, s_receive_atac_landscape_plot);
};

void SingleCellAtacItem::s_receive_atac_landscape_plot(ATAC_LANDSCAPE_PLOT_ELEMENTS ele) {

	auto&& gs = this->draw_suite_->graph_settings_;

	int nrow = ele.mat.rows(), ncol = ele.mat.cols();

	QCustomPlot* draw_area = custom_plot::initialize_plot(gs);
	SoapTextElement* title = new SoapTextElement(draw_area, gs.get_title("ATAC Landscape"), gs.get_title_font());
	draw_area->plotLayout()->addElement(0, 0, title);

	QCPLayoutGrid* main_layout = new QCPLayoutGrid;
	draw_area->plotLayout()->addElement(1, 0, main_layout);

	QCPLayoutGrid* right_layout = new QCPLayoutGrid;

	draw_area->plotLayout()->addElement(1, 1, right_layout);

	auto legend_layout = custom_plot::patch::set_legend_layout(draw_area, right_layout);

	if (ele.scale) {
		custom_plot::add_gradient_legend(
			draw_area,
			legend_layout,
			-1.0,
			1.0,
			"Accessibiliy",			
			gs,
			custom_plot::color::navy,
			Qt::white,
			custom_plot::color::firebrick3,
			"Low",
			"High"
		);
	}
	else {

		custom_plot::add_gradient_legend(draw_area, legend_layout, ele.mat.minCoeff(), ele.mat.maxCoeff(), "Accessibility", gs);
	}

	QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	axis_rect->setMarginGroup(QCP::msLeft | QCP::msRight, margin_group);
	axis_rect->setMarginGroup(QCP::msTop | QCP::msBottom, margin_group);
	QCPAxisRect* left_bottom_legend = new QCPAxisRect(draw_area, true);

	main_layout->addElement(0, 0, left_bottom_legend);
	main_layout->addElement(0, 1, axis_rect);

	auto cluster_names = custom::sapply(ele.cell_loc, [](auto&& loc) {return std::get<0>(loc); });

	auto cluster_colors = gs.palette(cluster_names);
	int index = 0;
	for (const auto& [cluster, start, n] : ele.cell_loc) {

		double end = start + n;

		custom_plot::patch::rectangle_borderless(
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

	custom_plot::patch::remove_left_bottom_axis(left_bottom_legend);
	left_bottom_legend->setMinimumSize(custom_plot::utility::get_max_text_width(custom::sapply(ele.cell_loc, [](auto&& t) {return std::get<0>(t); }), gs.get_left_label_font()) * 1.2, 200);

	left_bottom_legend->axis(QCPAxis::atBottom)->setRange(-11, 1);
	left_bottom_legend->axis(QCPAxis::atLeft)->setRange(0, nrow);


	QCPColorMap* heatmap = new QCPColorMap(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));
	heatmap->data()->setSize(ncol, nrow);
	heatmap->data()->setRange(QCPRange(0, ncol - 1), QCPRange(0, nrow - 1));
	custom_plot::patch::remove_left_bottom_axis(axis_rect);
	custom_plot::patch::set_range(axis_rect, QCPRange(-0.5, ncol - 0.5), QCPRange(-0.5, nrow - 0.5));

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
		gradient.setColorStopAt(1.0, custom_plot::color::firebrick3);
		gradient.setColorStopAt(0.0, custom_plot::color::navy);
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

	auto chrColors = custom_plot::utility::kmeans_palette(ele.chr_loc.size());
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


void SingleCellAtacItem::s_call_peaks_by_macs() {
	G_GETLOCK;

	auto fragments = this->data()->fragments();
	if (fragments == nullptr) {
		G_WARN("No Fragments Loaded.");
		G_UNLOCK;
		return;
	}

	MacsCallPeakWorker* worker = new MacsCallPeakWorker({ fragments });
	G_LINK_WORKER_THREAD(MacsCallPeakWorker, x_genomic_range_ready, SingleCellAtacItem, s_receive_macs_peaks);
}

void SingleCellAtacItem::s_receive_macs_peaks(GenomicRange genomic_range) {

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
};


void SingleCellAtacItem::update_quality_control_information() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		return;
	}

	Eigen::ArrayXi peak_count = custom::col_sum_mt(counts->mat_);
	const int ncol_atac = counts->mat_.cols();
	Eigen::ArrayXi peak_gene(ncol_atac);

	for (std::size_t i = 0; i < ncol_atac; ++i) {
		peak_gene[i] = counts->mat_.outerIndexPtr()[i + 1] - counts->mat_.outerIndexPtr()[i];
	}

	Metadata& metadata = *this->data()->metadata();

	metadata.mat_.update(METADATA_ATAC_UMI_NUMBER, custom::cast<QVector>(peak_count));
	metadata.mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, custom::cast<QVector>(peak_gene));
}


void SingleCellAtacItem::s_calculate_gene_activity() {

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
	G_LINK_WORKER_THREAD(CalculateGeneActivityWorker, x_gene_activity_ready, SingleCellAtacItem, s_receive_gene_activity);

};

void SingleCellAtacItem::s_receive_gene_activity(SparseInt* counts) {

	auto gene_activity = this->gene_activity_counts();

	if (gene_activity != nullptr) {
		gene_activity->__remove_this();
	}

	auto name = this->signal_emitter_->get_unique_name(VARIABLE_GENE_ACTIVITY);

	counts->data_type_ = SparseInt::DataType::GeneActivity;
	DATA_SUBMODULES(SparseInt)[name] = std::move(*counts);
	
	delete counts;

	auto item = new SparseIntItem(
		name,
		this->index_tree_,
		&DATA_SUBMODULES(SparseInt)[name],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);
};


void SingleCellAtacItem::s_receive_fast_annotation(
	QStringList main_type,
	QStringList sub_type,
	QString main_type_name,
	QString sub_type_name
) {

	this->data()->metadata()->mat_.update(main_type_name, main_type, CustomMatrix::DataType::QStringFactor);
	this->data()->metadata()->mat_.update(sub_type_name, sub_type, CustomMatrix::DataType::QStringFactor);

	this->signal_emitter_->x_update_interface();
};

void SingleCellAtacItem::s_fast_annotation() {
	G_GETLOCK;

	auto gene_activity = this->data()->gene_activity_counts();
	if (gene_activity == nullptr) {
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
			gene_activity,
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
				gene_activity,
				main_type_name,
				sub_type_name,
				true,
				this->data()->metadata()->mat_.get_qstring(settings[3])
			);
		}
		else {

			worker = new CellTypeAnnotationWorker(
				gene_activity,
				main_type_name,
				sub_type_name,
				false
			);
		}

	}

	G_LINK_WORKER_THREAD(CellTypeAnnotationWorker, x_annotation_ready, SingleCellAtacItem, s_receive_fast_annotation);
};


void SingleCellAtacItem::s_find_motifs() {

	G_GETLOCK;

	auto atac_counts = this->data()->counts();
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
	G_LINK_WORKER_THREAD(MotifLocateWorker, x_motif_location_ready, SingleCellAtacItem, s_receive_motif_location);
};

void SingleCellAtacItem::s_receive_motif_location(MotifPosition motif_position) {

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
};


void SingleCellAtacItem::s_receive_integrated_data(SingleCellAtac* data, QList<const SingleCellAtac* > items) {

	this->signal_emitter_->unlock(this->signal_emitter_->search(custom::sapply(items, [](auto* data) {return (void*)data; })));

	this->signal_emitter_->x_data_create_soon(data, soap::VariableType::SingleCellAtac, "Integrated scAtac");
};

void SingleCellAtacItem::s_integrate() {

	QMap<QString, SingleCellAtac*> available_data;
	for (const auto& [variable_name, data_info] : this->signal_emitter_->variable_information_) {
		if (data_info.first == soap::VariableType::SingleCellAtac) {
			available_data[variable_name] = static_cast<SingleCellAtac*>(data_info.second);
		}
	}
	if (available_data.size() < 2) {
		G_LOG("No enough single cell atac data found.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Integrate SingleCellAtac",
		{ "Data", "Style" },
		{ soap::InputStyle::SimpleChoice, soap::InputStyle::ComboBox},
		{ available_data.keys(), 
		{ "Distinguish cells from different data", "Merge the same cell name" } }
	);

	if (settings.isEmpty())return;

	std::map<QString, SingleCellAtac*> choosed_data;

	auto choosed_name = simple_choice_to_list(settings[0]);

	if (choosed_name.size() < 2)return;

	if (!custom::is_unique(choosed_name)) {
		G_WARN("Can not integrate the same data!");
		return;
	}

	std::ranges::for_each(choosed_name, 
		[&choosed_data, &available_data](const QString& name) {choosed_data[name] = available_data[name]; });

	QList<const SingleCellAtac*> all_data;
	int velocyto_base_loaded = 0;

	QStringList fragments_missed_data;
	for (const auto& [data_name, ptr] : choosed_data) {

		auto fragments = ptr->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments of " + data_name + " has not been loaded.");
			return;
		}

		all_data << ptr;
	}

	auto species = custom::unique(custom::sapply(all_data, [](const SingleCellAtac* data) {return data->species_; }));

	if (species.size() != 1) {
		G_WARN("Unmatched species!");
		return;
	}

	if (!this->signal_emitter_->try_lock(custom::sapply(all_data, 
		[this](const SingleCellAtac* data) {return this->signal_emitter_->search((void*)data); }))) {
		G_WARN("Please waiting for computation in progress.");
		return;
	}

	bool distinguish = (settings[1] == "Distinguish cells from different data");
	if (!distinguish) {
		G_NOTICE("Metadata will not be inherited to avoid collision.");
	}

	G_LOG("Start integrating data...");
	IntegrateWorker* worker = new IntegrateWorker(all_data, distinguish, species[0]);
	G_LINK_WORKER_THREAD(IntegrateWorker, x_scatac_ready, SingleCellAtacItem, s_receive_integrated_data);
};

void SingleCellAtacItem::s_load_fragments() {
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
		barcodes = this->data()->counts()->colnames_;
	}
	LoadFragmentsWorker* worker = new LoadFragmentsWorker(barcodes, fragments_file_path);
	G_LINK_WORKER_THREAD(LoadFragmentsWorker, x_fragments_ready, SingleCellAtacItem, s_receive_loaded_fragments);
};

void SingleCellAtacItem::s_receive_loaded_fragments(Fragments* fragments) {

	fragments->cell_names_ = this->data()->counts()->colnames_;

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


void SingleCellAtacItem::s_monocle3() {

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
			{ soap::InputStyle::ComboBox },
			{ embedding_names }
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

		G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellAtacItem, s_receive_monocle3);
	}
	else {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Monocle3 Settings",
			{ "Embedding" , "use part", "choose cell group" },
			{ soap::InputStyle::ComboBox, soap::InputStyle::SwitchButton, soap::InputStyle::FactorChoice },
			{ embedding_names },
			{ factor_info }
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

			G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellAtacItem, s_receive_monocle3);
		}
		else {
			int n_cell = this->data()->counts()->cols();

			Monocle3Worker* worker = new Monocle3Worker(
				DATA_SUBMODULES(Embedding)[embedding_name],
				Eigen::ArrayX<bool>::Constant(n_cell, true)
			);

			G_LINK_WORKER_THREAD(Monocle3Worker, x_monocle3_ready, SingleCellAtacItem, s_receive_monocle3);
		}
	}

};

void SingleCellAtacItem::s_receive_monocle3(Monocle3* monocle3) {

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
};