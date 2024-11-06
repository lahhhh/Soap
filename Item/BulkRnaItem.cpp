#include "BulkRnaItem.h"

#include "YesOrNoDialog.h"
#include "CommonDialog.h"
#include "StatisticsDialog.h"

#include "ItemIOWorker.h"
#include "RnaNormalizeWorker.h"
#include "PcaWorker2.h"
#include "TsneWorker.h"
#include "UmapWorker.h"
#include "IntegrateWorker.h"
#include "LeidenPartitionWorker.h"
#include "SlmWorker.h"
#include "DifferentialAnalysisWorker.h"
#include "GseaWorker.h"

void BulkRnaItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Metadata));
	this->check_variable(DATA_SUBMODULES(DenseInt));
	this->check_variable(DATA_SUBMODULES(DenseDouble));
	this->check_variable(DATA_SUBMODULES(Embedding));
	this->check_variable(DATA_SUBMODULES(DifferentialAnalysis));
	this->check_variable(DATA_SUBMODULES(GSEA));

	auto item = new NoteItem(&this->data()->string_information_["Note"], this->data(), this->signal_emitter_);
	this->addChild(item);
};

void BulkRnaItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Normalize");
	ADD_ACTION("Normalize(log transform)", "Normalize", s_normalize1);
	ADD_ACTION("FPKM(RPKM)", "Normalize", s_fpkm);
	ADD_ACTION("TPM", "Normalize", s_tpm);

	//setting
	ADD_MAIN_MENU("Settings");

	ADD_ACTION("Random Seed", "Settings", s_set_random_state);
	ADD_ACTION("Species", "Settings", s_set_species);


	// Dimensional Reduction
	ADD_MAIN_MENU("Dimension Reduction");

	ADD_MENU("Dimension Reduction | PCA", "PCA", "Dimension Reduction");
	ADD_ACTION("Default", "Dimension Reduction | PCA", s_pca_default);
	ADD_ACTION("Custom", "Dimension Reduction | PCA", s_pca_custom);

	ADD_MENU("Dimension Reduction | t-SNE", "t-SNE", "Dimension Reduction");
	ADD_ACTION("Default", "Dimension Reduction | t-SNE", s_tsne_default);
	ADD_ACTION("Custom", "Dimension Reduction | t-SNE", s_tsne_custom);

	ADD_MENU("Dimension Reduction | UMAP", "UMAP", "Dimension Reduction");
	ADD_ACTION("Default", "Dimension Reduction | UMAP", s_umap_default);
	ADD_ACTION("Custom", "Dimension Reduction | UMAP", s_umap_custom);


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

	// DEG
	ADD_MAIN_MENU("Expression Analysis");

	ADD_ACTION("Find DEG", "Expression Analysis", s_find_deg);

	ADD_MAIN_ACTION("Integrate", s_integrate);

	ADD_MAIN_ACTION("Duplicate", __s_duplicate);

	ADD_MAIN_ACTION("Save", __s_export_as_item);

	// GSEA
	ADD_MAIN_ACTION("GSEA", s_gsea);

	// delete
	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void BulkRnaItem::__s_rename() {
	G_GETLOCK;

	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_,
		"Rename BulkRna",
		{ "New Name:" + this->title_, "Overwrite Metadata:Source:yes" },
		{ soap::InputStyle::StringLineEdit, soap::InputStyle::SwitchButton }
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

void BulkRnaItem::s_statistics() {

	StatisticsDialog::get_response(this->data());
};

void BulkRnaItem::__s_delete_this() {

	G_GETLOCK;
	G_UNLOCK;

	if (!YesOrNoDialog::get_response("Delete Bulk RNA Data", "This data will be deleted.")) {
		return;
	}

	this->__remove_this();
};

void BulkRnaItem::s_normalize1() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("No Counts Data.");
		return;
	}

	G_GETLOCK;

	if (auto normalized = this->normalized()) {
		normalized->__remove_this();
	}

	auto* worker = new RnaNormalizeWorker(
		counts,
		"Normalize1"
	);

	G_LINK_WORKER_THREAD(RnaNormalizeWorker, x_normalize_ready, BulkRnaItem, s_receive_normalize);
};

void BulkRnaItem::s_fpkm() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("No Counts Data.");
		return;
	}

	G_GETLOCK;

	if (auto normalized = this->normalized()) {
		normalized->__remove_this();
	}

	auto* worker = new RnaNormalizeWorker(
		counts,
		"FPKM"
	);

	G_LINK_WORKER_THREAD(RnaNormalizeWorker, x_normalize_ready, BulkRnaItem, s_receive_normalize);
};

void BulkRnaItem::s_tpm() {

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("No Counts Data.");
		return;
	}

	G_GETLOCK;

	if (auto normalized = this->normalized()) {
		normalized->__remove_this();
	}

	auto* worker = new RnaNormalizeWorker(
		counts,
		"TPM"
	);

	G_LINK_WORKER_THREAD(RnaNormalizeWorker, x_normalize_ready, BulkRnaItem, s_receive_normalize);
};

void BulkRnaItem::s_receive_normalize(DenseDouble normalized) {

	normalized.data_type_ = DenseDouble::DataType::Normalized;

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_NORMALIZED);

	DATA_SUBMODULES(DenseDouble)[title] = normalized;

	auto item = new DenseDoubleItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(DenseDouble)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);

	G_LOG("Normalize finished.");
};

void BulkRnaItem::s_pca_default() {

	G_GETLOCK;

	auto normalized = this->data()->normalized();
	if (normalized == nullptr) {
		G_WARN("Data must be normalized.");
		return;
	}

	if (auto item = this->pca()) {
		item->__remove_this();
	}

	auto* worker = new PcaWorker2(normalized->mat_);
	G_LINK_WORKER_THREAD(PcaWorker2, x_pca_ready, BulkRnaItem, s_receive_pca);
};

void BulkRnaItem::s_pca_custom() {

	G_GETLOCK;

	auto normalized = this->data()->normalized();
	if (normalized == nullptr) {
		G_WARN("Data must be normalized.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"PCA settings",
		{"Use Variable Features:yes", "Number of Variable Features:2000"},
		{soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit}
	);

	if (settings.isEmpty()) {
		return;
	}

	bool use_variable_features = switch_to_bool(settings[0]);
	int n_variable_features = settings[1].toInt();

	if (use_variable_features) {
		if (n_variable_features < 1000) {
			G_WARN("Two few features.");
			G_UNLOCK;
			return;
		}
	}

	if (auto item = this->pca()) {
		item->__remove_this();
	}

	auto* worker = new PcaWorker2(normalized->mat_, use_variable_features, n_variable_features);
	G_LINK_WORKER_THREAD(PcaWorker2, x_pca_ready, BulkRnaItem, s_receive_pca);
};

void BulkRnaItem::s_receive_pca(Eigen::MatrixXd emb, QVector<double> sdev, QVector<double> vp) {

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_PCA);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Pca,
		emb,
		this->data()->counts()->colnames_,
		custom::paste("PCA-", custom::cast<QString>(custom::seq_n(1, emb.cols())))
	);

	this->data()->double_vectors_[VARIABLE_PCA_STANDARD_DEVIATION] = sdev;
	this->data()->double_vectors_[VARIABLE_PCA_VARIANCE_PROPORTION] = vp;

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
};


void BulkRnaItem::s_tsne_default() {
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

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, BulkRnaItem, s_receive_tsne)
};


void BulkRnaItem::s_tsne_custom() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("No PCA Data Found.");
		return;
	}

	G_GETLOCK;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"tSNE settings",
		{
			"Dimension start:1",
			"Dimension end:20",
			"Random State:" + QString::number(this->data()->random_state_)
		},
		{
		soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit
		}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto&& mat = pca->data_.mat_;
	const int nrow = mat.rows(), ncol = mat.cols();

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

	unsigned int random_state = settings[2].toUInt();

	if (auto item = this->tsne()) {
		item->__remove_this();
	}

	Eigen::MatrixXd tsne_input = mat.block(0, start_dimension - 1, nrow, last_dimension - start_dimension + 1);

	G_LOG("Customized tSNE start using PCA (" + 
		QString::number(start_dimension) + "~" + QString::number(last_dimension) + ") dimensions.");


	TsneWorker* worker = new TsneWorker(tsne_input, random_state);

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, BulkRnaItem, s_receive_tsne);
};

void BulkRnaItem::s_receive_tsne(Eigen::MatrixXd mat) {

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

void BulkRnaItem::s_umap_default() {
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

	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, BulkRnaItem, s_receive_umap)
}

void BulkRnaItem::s_umap_custom() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("No PCA Data Found.");
		return;
	}

	G_GETLOCK;

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
			"Random state:" + QString::number(this->data()->random_state_)
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
			soap::InputStyle::IntegerLineEdit
		},
		{{ "Angular", "Euclidean", "Manhattan" },{ "Random" }}
		);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto&& mat = pca->data_.mat_;
	const int nrow = mat.rows(), ncol = mat.cols();

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
	Eigen::MatrixXd umap_input = mat.block(0, start_dimension - 1, nrow, last_dimension - start_dimension + 1);

	G_LOG("Customized UMAP start using PCA (" +
		QString::number(start_dimension) + "~" + QString::number(last_dimension) + ") dimensions.");

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

	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, BulkRnaItem, s_receive_umap)
}

void BulkRnaItem::s_receive_umap(Eigen::MatrixXd mat) {

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

void BulkRnaItem::s_receive_integrated_data(BulkRna* data, QList<const BulkRna* > items) {

	this->signal_emitter_->unlock(custom::sapply(items, [this](auto* data) {return this->signal_emitter_->search((void*)data); }));

	this->signal_emitter_->x_data_create_soon(data, soap::VariableType::BulkRna, "Integrated scRNA");
};

void BulkRnaItem::s_integrate() {

	auto& variables = this->signal_emitter_->variable_information_;

	QStringList available_data = this->signal_emitter_->get_type_variable(soap::VariableType::BulkRna).keys();
	if (available_data.size() < 2) {
		G_LOG("No enough single cell data found.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Integrate BulkRna",
		{ "Data" },
		{ soap::InputStyle::SimpleChoice },
		{ available_data}
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

	auto choosed = custom::sapply(its, [](IndexTree* it) {return static_cast<const BulkRna*>(it->data_); });

	auto species = custom::unique(custom::sapply(choosed, [](const BulkRna* rna) {return rna->species_; }));

	if (species.size() != 1) {
		G_LOG("Unmatched species!");
		return;
	}

	G_LOG("Start integrating data...");

	IntegrateWorker* worker = new IntegrateWorker(choosed, species[0]);

	G_LINK_WORKER_THREAD(IntegrateWorker, x_bulkrna_ready, BulkRnaItem, s_receive_integrated_data);
};


void BulkRnaItem::s_leiden_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	LeidenPartitionWorker* worker = new LeidenPartitionWorker(pca->data_.mat_, "Modularity", "Euclidean", 30, 50, 0.6);
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, BulkRnaItem, s_receive_leiden);
}

void BulkRnaItem::s_louvain_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "Louvain", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, BulkRnaItem, s_receive_louvain);
}

void BulkRnaItem::s_modified_louvain_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "Modified Louvain", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, BulkRnaItem, s_receive_modified_louvain);
}

void BulkRnaItem::s_smart_local_moving_default() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(pca->data_.mat_, "SLM", "Euclidean", 1, 10, 10, 1997, 30, 50, 0.6);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, BulkRnaItem, s_receive_slm);
}

void BulkRnaItem::s_smart_local_moving_custom() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("No PCA Data Found.");
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
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox
		},

		{ { "Euclidean", "Angular", "Manhattan" }, { "1", "2" } }
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

	Eigen::MatrixXd* from_matrix = &pca->data_.mat_;

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	int n_start = settings[5].toInt();
	if (n_start < 1 || n_start > 100) {
		G_WARN("Number of start should be between 1 and 100, reset to 10");
		n_start = 10;
	}

	int n_iteration = settings[6].toInt();
	if (n_iteration < 1 || n_iteration > 100) {
		G_WARN("Number of iteration should be between 1 and 100, reset to 10");
		n_iteration = 10;
	}

	int modularity_function = settings[7].toInt();

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(*from_matrix, "SLM", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, BulkRnaItem, s_receive_slm);
}

void BulkRnaItem::s_modified_louvain_custom() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("No PCA Data Found.");
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
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox
		},
		{ { "Euclidean", "Angular", "Manhattan" }, { "1", "2" } }
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

	Eigen::MatrixXd* from_matrix = &pca->data_.mat_;

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	int n_start = settings[5].toInt();
	if (n_start < 1 || n_start > 100) {
		G_WARN("Number of start should be between 1 and 100, reset to 10");
		n_start = 10;
	}

	int n_iteration = settings[6].toInt();
	if (n_iteration < 1 || n_iteration > 100) {
		G_WARN("Number of iteration should be between 1 and 100, reset to 10");
		n_iteration = 10;
	}

	int modularity_function = settings[7].toInt();

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(*from_matrix, "Modified Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, BulkRnaItem, s_receive_modified_louvain);
}

void BulkRnaItem::s_leiden_custom() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("No PCA Data Found.");
		return;
	}

	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Leiden Clustering Settings",
		{ 
			"Method", 
			"Metric",
			"Resolution:0.6", 
			"Number of Neighbors:30",
			"Number of Trees:50"
		},
		{ 
			soap::InputStyle::ComboBox, 
			soap::InputStyle::ComboBox, 
			soap::InputStyle::NumericLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit },
		{ 
			{ "Modularity", "CPM", "RBConfiguration", "RBER", "Significance", "Surprise" }, 
		{ "Euclidean", "Angular", "Manhattan" } 
		}
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

	Eigen::MatrixXd* from_matrix = &pca->data_.mat_;

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	G_LOG("Start Clustering...");

	LeidenPartitionWorker* worker = new LeidenPartitionWorker(*from_matrix, method, metric, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, BulkRnaItem, s_receive_leiden);
}

void BulkRnaItem::s_louvain_custom() {

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("No PCA Data Found.");
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
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::ComboBox
		},
		{ { "Euclidean", "Angular", "Manhattan" }, { "1", "2" } }
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

	Eigen::MatrixXd* from_matrix = &pca->data_.mat_;

	if (n_neighbors <= 0 || n_neighbors > from_matrix->cols()) {
		n_neighbors = from_matrix->cols() > 15 ? 15 : from_matrix->cols();
		G_WARN("Illegal number of neighbors! reset to " + QString::number(n_neighbors));
	}

	int n_start = settings[5].toInt();
	if (n_start < 1 || n_start > 100) {
		G_WARN("Number of start should be between 1 and 100, reset to 10");
		n_start = 10;
	}

	int n_iteration = settings[6].toInt();
	if (n_iteration < 1 || n_iteration > 100) {
		G_WARN("Number of iteration should be between 1 and 100, reset to 10");
		n_iteration = 10;
	}

	int modularity_function = settings[7].toInt();

	G_LOG("Start Clustering...");
	SlmWorker* worker = new SlmWorker(*from_matrix, "Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, BulkRnaItem, s_receive_louvain);
}

void BulkRnaItem::s_receive_leiden(QVector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_RNA_LEIDEN_CLUSTER, cluster, CustomMatrix::DataType::IntegerFactor);
};

void BulkRnaItem::s_receive_louvain(std::vector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_RNA_LOUVAIN_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void BulkRnaItem::s_receive_modified_louvain(std::vector<int> cluster) {

	this->data()->metadata()->mat_.update(METADATA_RNA_MODIFIED_LOUVAIN_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void BulkRnaItem::s_receive_slm(std::vector<int> cluster) {


	this->data()->metadata()->mat_.update(METADATA_RNA_SLM_CLUSTER, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};


void BulkRnaItem::s_find_deg() {

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
		{ soap::InputStyle::CompareLayout, soap::InputStyle::NumericLineEdit, soap::InputStyle::ComboBox },
		{ { "Bonferroni", "FDR" } },
		{ factor_map }
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
	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, BulkRnaItem, s_receive_differential_analysis);
};

void BulkRnaItem::s_receive_differential_analysis(DifferentialAnalysis da, QString name) {

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


void BulkRnaItem::s_receive_gsea(GSEA gsea) {

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

void BulkRnaItem::s_gsea() {

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
		"Permutation Type" },
		{ soap::InputStyle::CompareLayout, soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::ComboBox },
		{ { "Curated", "Ontology" }, { "Phenotype", "Geneset" } },
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

	QStringList metadata = this->data()->metadata()->mat_.get_qstring(factor_name);

	G_LOG("GSEAing in " + factor_name + "...");
	SparseDouble tmp = normalized->to_sparse();

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
	G_LINK_WORKER_THREAD(GseaWorker, x_gsea_ready, BulkRnaItem, s_receive_gsea);
}

void BulkRnaItem::s_set_species() {

	G_GETLOCK;
	G_UNLOCK;

	QStringList species = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Species",
		{ "New Species" },
		{ soap::InputStyle::ComboBox },
		{ { "Human", "Mouse", "Undefined" } }
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

void BulkRnaItem::s_set_random_state() {

	G_GETLOCK;
	G_UNLOCK;

	QStringList seed = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Random Seed",
		{ "New Random Seed:" + QString::number(this->data()->random_state_) },
		{ soap::InputStyle::StringLineEdit }
	);

	if (seed.isEmpty()) {
		return;
	}

	this->data()->random_state_ = seed[0].toInt();
};