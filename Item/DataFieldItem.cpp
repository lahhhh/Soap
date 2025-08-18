#include "DataFieldItem.h"

#include "CustomPlot.h"
#include "LogNormalizeWorker.h"
#include "TfidfWorker.h"
#include "CommonDialog.h"
#include "PcaWorker.h"
#include "UmapWorker.h"
#include "TsneWorker.h"
#include "HarmonyWorker.h"
#include "DifferentialAnalysisWorker.h"
#include "SlmWorker.h"
#include "LeidenPartitionWorker.h"
#include "SvdWorker.h"
#include "WilcoxTest.h"

#include "SingleCellMultiomeItem.h"

void DataFieldItem::__set_menu() {

	CREATE_ROOT_MENU;

	// visualize
	ADD_MAIN_MENU("Visualize");

	auto type = this->data()->data_type_;
	if (type == DataField::DataType::Rna ||	type == DataField::DataType::Atac) {

		ADD_ACTION("Quality View", "Visualize", s_view_quality);

		ADD_MAIN_ACTION("Re-estimate quality parameters", s_re_estimate_quality_parameters);
	}

	ADD_ACTION("Distribution Plot", "Visualize", s_distribution_plot);

	ADD_ACTION("Bubble Plot", "Visualize", s_bubble_plot);

	ADD_ACTION("Heatmap Plot", "Visualize", s_heatmap_plot);
	ADD_ACTION("Cell-wise Heatmap Plot", "Visualize", s_cell_heatmap_plot);

	ADD_ACTION("Correlation Heatmap Plot", "Visualize", s_correlation_heatmap_plot);

	//Normalize
	ADD_MAIN_MENU("Normalize");

	if (type == DataField::DataType::Atac) {
		ADD_MENU("Normalize | TFIDF", "TFIDF", "Normalize");

		ADD_ACTION("Default", "Normalize | TFIDF", s_tfidf_default);
	}

	ADD_MENU("Normalize | Log Normalize", "Log Normalize", "Normalize");

	ADD_ACTION("Default", "Normalize | Log Normalize", s_log_normalize_default);
	ADD_ACTION("Custom", "Normalize | Log Normalize", s_log_normalize_custom);

	// Dimensional Reduction
	ADD_MAIN_MENU("Dimension Reduction");

	if (type == DataField::DataType::Atac) {

		ADD_MENU("Dimension Reduction | SVD", "SVD", "Dimension Reduction");

		ADD_ACTION("Default", "Dimension Reduction | SVD", s_svd_default);
		ADD_ACTION("Custom", "Dimension Reduction | SVD", s_svd_custom);
	}

	ADD_MENU("Dimension Reduction | PCA", "PCA", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | PCA", s_pca_default);
	ADD_ACTION("Custom", "Dimension Reduction | PCA", s_pca_custom);

	ADD_MENU("Dimension Reduction | UMAP", "UMAP", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | UMAP", s_umap_default);
	ADD_ACTION("Custom", "Dimension Reduction | UMAP", s_umap_custom);

	ADD_MENU("Dimension Reduction | tSNE", "tSNE", "Dimension Reduction");

	ADD_ACTION("Default", "Dimension Reduction | tSNE", s_tsne_default);
	ADD_ACTION("Custom", "Dimension Reduction | tSNE", s_tsne_custom);

	// Remove Batch Effect
	ADD_MAIN_MENU("Remove Batch Effect");

	ADD_ACTION("Harmony", "Remove Batch Effect", s_harmony);

	//deg

	ADD_MAIN_MENU("Expression Analysis");

	if (type == DataField::DataType::Rna) {
		ADD_ACTION("Find DEG", "Expression Analysis", s_find_deg);
	}
	else if (type == DataField::DataType::Atac) {
		ADD_ACTION("Find DAP", "Expression Analysis", s_find_dap);
	}
	else if (type == DataField::DataType::Trans) {
		ADD_ACTION("Find DAG", "Expression Analysis", s_find_dag);
	}

	// filter
	ADD_MAIN_MENU("Filter");

	ADD_ACTION("By Features", "Filter", s_filter_by_feature);

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

}

void DataFieldItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(SparseInt));
	this->check_variable(DATA_SUBMODULES(SparseDouble));
	this->check_variable(DATA_SUBMODULES(Embedding));
	this->check_variable(DATA_SUBMODULES(DifferentialAnalysis));
};

SparseIntItem* DataFieldItem::counts() {
	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::SparseInt) {
			
			SparseInt* sparseint = static_cast<SparseInt*>(item->data_);
			if (sparseint->data_type_ == SparseInt::DataType::Counts) {
				return static_cast<SparseIntItem*>(item);
			}
		}
	}

	return nullptr;
};

Eigen::ArrayX<bool> find_variable_features(const Eigen::SparseMatrix<double>& mat) {

	const int ncol = mat.cols();
	const int nrow = mat.rows();
	double clipmax = std::sqrt(ncol);
	Eigen::ArrayXd row_mean = custom::row_mean(mat);
	Eigen::ArrayXd row_var = custom::row_var(mat);

	Eigen::ArrayX<bool> not_const = row_var > 0;
	const int not_const_row = not_const.count();
	if (not_const_row < 1000) {
		return Eigen::ArrayX<bool>();
	}
	auto not_const_index = custom::which(not_const);

	Eigen::ArrayXd not_const_row_means = custom::sliced(row_mean, not_const);
	Eigen::ArrayXd not_const_row_var = custom::sliced(row_var, not_const);

	Eigen::ArrayXd fitted = custom::loess_mt(log10(not_const_row_var), log10(not_const_row_means), 2, 0.3, 50);

	Eigen::ArrayXd expected_var = Eigen::ArrayXd::Zero(nrow);
	expected_var(not_const_index) = pow(10, fitted);

	Eigen::ArrayXd expected_sd = sqrt(expected_var);
	Eigen::ArrayXd standardized_row_var = Eigen::ArrayXd::Zero(nrow);
	Eigen::ArrayXi n_zero = Eigen::ArrayXi::Constant(nrow, ncol);

	for (int i = 0; i < ncol; ++i) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
		{
			--n_zero[it.row()];
			standardized_row_var[it.row()] += std::pow(std::min(clipmax, std::abs((it.value() - row_mean[it.row()]) / expected_sd[it.row()])), 2);
		}
	}

	for (int i = 0; i < nrow; ++i) {
		if (expected_sd[i] == 0)continue;
		standardized_row_var[i] += std::pow(row_mean[i] / expected_sd[i], 2) * n_zero[i];
	}
	standardized_row_var /= (ncol - 1);
	Eigen::ArrayXi standardized_row_var_sorted_index = custom::order(standardized_row_var, true);
	Eigen::ArrayX<bool> res = Eigen::ArrayX<bool>::Constant(nrow, false);


	int n_feature = 2000;

	for (int i = 0; i < n_feature; ++i) {
		int index = standardized_row_var_sorted_index[i];
		res[index] = true;
	}
	return res;
}

SparseDoubleItem* DataFieldItem::normalized() {
	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::SparseDouble) {
			SparseDouble* sparsedouble = static_cast<SparseDouble*>(item->data_);
			if (sparsedouble->data_type_ == SparseDouble::DataType::Normalized) {
				return static_cast<SparseDoubleItem*>(item);
			}
		}
	}

	return nullptr;
};

EmbeddingItem* DataFieldItem::pca() {
	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::Embedding) {
			Embedding* embedding = static_cast<Embedding*>(item->data_);
			auto embedding_type = embedding->data_type_;
			if (embedding_type == Embedding::DataType::Pca) {
				return static_cast<EmbeddingItem*>(item);
			}
		}
	}

	return nullptr;
};

EmbeddingItem* DataFieldItem::tsne() {
	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::Embedding) {
			Embedding* embedding = static_cast<Embedding*>(item->data_);
			auto embedding_type = embedding->data_type_;
			if (embedding_type == Embedding::DataType::Tsne) {
				return static_cast<EmbeddingItem*>(item);
			}
		}
	}

	return nullptr;
};

EmbeddingItem* DataFieldItem::umap() {
	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::Embedding) {
			Embedding* embedding = static_cast<Embedding*>(item->data_);
			auto embedding_type = embedding->data_type_;
			if (embedding_type == Embedding::DataType::Umap) {
				return static_cast<EmbeddingItem*>(item);
			}
		}
	}

	return nullptr;
};

EmbeddingItem* DataFieldItem::harmony() {
	const int child_count = this->childCount();

	for (int i = 0; i < child_count; ++i) {
		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		if (item->data_type_ == soap::VariableType::Embedding) {
			Embedding* embedding = static_cast<Embedding*>(item->data_);
			auto embedding_type = embedding->data_type_;
			if (embedding_type == Embedding::DataType::Harmony) {
				return static_cast<EmbeddingItem*>(item);
			}
		}
	}

	return nullptr;
};

void DataFieldItem::s_receive_differential_analysis(DifferentialAnalysis da, QString name) {

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

void DataFieldItem::s_find_dag() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	if (this->counts() == nullptr) {
		G_WARN("No Gene Activity Data found.");
		return;
	}
	G_GETLOCK;

	auto factor_map = single_cell_multiome->metadata()->mat_.get_factor_information(false);
	if (factor_map.isEmpty()) {
		G_LOG("No suitable feature for Find DAG");
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

	QStringList metadata = single_cell_multiome->metadata()->mat_.get_qstring(factor_name);

	double minimum_percentage = settings[1].toDouble();
	if (minimum_percentage < 0.0 || minimum_percentage >= 1.0) {
		G_WARN("Invalid Percentage! Reset to 0.1.");
		minimum_percentage = 0.1;
	}

	QString p_adjust_method = settings[2];

	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		minimum_percentage,
		p_adjust_method
	);
	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, DataFieldItem, s_receive_differential_analysis)

};

void DataFieldItem::s_find_dap() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	if (this->counts() == nullptr) {
		G_WARN("No ATAC Data found.");
		return;
	}

	G_GETLOCK;

	auto factor_map = single_cell_multiome->metadata()->mat_.get_factor_information(false);
	if (factor_map.isEmpty()) {
		G_LOG("No suitable feature for Find DAP");
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

	QStringList metadata = single_cell_multiome->metadata()->mat_.get_qstring(factor_name);

	double minimum_percentage = settings[1].toDouble();
	if (minimum_percentage < 0.0 || minimum_percentage >= 1.0) {
		G_WARN("Invalid Percentage! Reset to 0.1.");
		minimum_percentage = 0.1;
	}

	QString p_adjust_method = settings[2];

	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		minimum_percentage,
		p_adjust_method
	);
	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, DataFieldItem, s_receive_differential_analysis)
};

void DataFieldItem::s_find_deg() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	if (this->counts() == nullptr) {
		G_WARN("No RNA Data found.");
		return;
	}
	G_GETLOCK;

	auto factor_map = single_cell_multiome->metadata()->mat_.get_factor_information(false);
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

	QStringList metadata = single_cell_multiome->metadata()->mat_.get_qstring(factor_name);

	double minimum_percentage = settings[1].toDouble();
	if (minimum_percentage < 0.0 || minimum_percentage >= 1.0) {
		G_WARN("Invalid Percentage! Reset to 0.1.");
		minimum_percentage = 0.1;
	}

	QString p_adjust_method = settings[2];

	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		minimum_percentage,
		p_adjust_method
	);
	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, DataFieldItem, s_receive_differential_analysis)

};

void DataFieldItem::s_receive_umap(Eigen::MatrixXd mat) {

	QString umap_name;

	if (this->data()->data_type_ == DataField::DataType::Rna) {
		umap_name = VARIABLE_RNA_UMAP;
	}
	else if (this->data()->data_type_ == DataField::DataType::Atac) {
		umap_name = VARIABLE_ATAC_UMAP;
	}
	else if (this->data()->data_type_ == DataField::DataType::Trans) {
		umap_name = VARIABLE_TRANS_UMAP;
	}

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

void DataFieldItem::s_receive_svd(Eigen::MatrixXd mat, QVector<double> sd) {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	QString pca_name, sd_name;

	if (this->data()->data_type_ == DataField::DataType::Rna) {
		pca_name = VARIABLE_RNA_PCA;
		sd_name = VARIABLE_RNA_PCA_STANDARD_DEVIATION;
	}
	else if (this->data()->data_type_ == DataField::DataType::Atac) {
		pca_name = VARIABLE_ATAC_PCA;
		sd_name = VARIABLE_ATAC_PCA_STANDARD_DEVIATION;
	}
	else if (this->data()->data_type_ == DataField::DataType::Trans) {
		pca_name = VARIABLE_TRANS_PCA;
		sd_name = VARIABLE_TRANS_PCA_STANDARD_DEVIATION;
	}

	QString title = this->signal_emitter_->get_unique_name(pca_name);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Pca,
		mat,
		this->data()->counts()->colnames_,
		custom::paste(pca_name, custom::cast<QString>(custom::seq_n(1, mat.cols())), "-")
	);

	single_cell_multiome->double_vectors_[sd_name] = sd;

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

void DataFieldItem::s_receive_pca(Eigen::MatrixXd mat, QVector<double> sd) {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	QString pca_name, sd_name;

	if (this->data()->data_type_ == DataField::DataType::Rna) {
		pca_name = VARIABLE_RNA_PCA;
		sd_name = VARIABLE_RNA_PCA_STANDARD_DEVIATION;
	}
	else if (this->data()->data_type_ == DataField::DataType::Atac) {
		pca_name = VARIABLE_ATAC_PCA;
		sd_name = VARIABLE_ATAC_PCA_STANDARD_DEVIATION;
	}
	else if (this->data()->data_type_ == DataField::DataType::Trans) {
		pca_name = VARIABLE_TRANS_PCA;
		sd_name = VARIABLE_TRANS_PCA_STANDARD_DEVIATION;
	}

	QString title = this->signal_emitter_->get_unique_name(pca_name);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Pca,
		mat,
		this->data()->counts()->colnames_,
		custom::paste(pca_name, custom::cast<QString>(custom::seq_n(1, mat.cols())), "-")
	);

	single_cell_multiome->double_vectors_[sd_name] = sd;

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

void DataFieldItem::s_receive_leiden(QVector<int> cluster) {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	QString cluster_name;

	auto type = this->data()->data_type_;
	switch (type)
	{
	case DataField::DataType::Rna:
		cluster_name = METADATA_RNA_LEIDEN_CLUSTER;
		break;
	case DataField::DataType::Atac:
		cluster_name = METADATA_ATAC_LEIDEN_CLUSTER;
		break;
	case DataField::DataType::Trans:
		cluster_name = METADATA_TRANS_LEIDEN_CLUSTER;
		break;
	default:
		break;
	}

	single_cell_multiome->metadata()->mat_.update(cluster_name, cluster, CustomMatrix::DataType::IntegerFactor);
};

void DataFieldItem::s_receive_louvain(std::vector<int> cluster) {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	QString cluster_name;

	auto type = this->data()->data_type_;
	switch (type)
	{
	case DataField::DataType::Rna:
		cluster_name = METADATA_RNA_LOUVAIN_CLUSTER;
		break;
	case DataField::DataType::Atac:
		cluster_name = METADATA_ATAC_LOUVAIN_CLUSTER;
		break;
	case DataField::DataType::Trans:
		cluster_name = METADATA_TRANS_LOUVAIN_CLUSTER;
		break;
	default:
		break;
	}

	single_cell_multiome->metadata()->mat_.update(cluster_name, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void DataFieldItem::s_receive_modified_louvain(std::vector<int> cluster) {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	QString cluster_name;

	auto type = this->data()->data_type_;
	switch (type)
	{
	case DataField::DataType::Rna:
		cluster_name = METADATA_RNA_MODIFIED_LOUVAIN_CLUSTER;
		break;
	case DataField::DataType::Atac:
		cluster_name = METADATA_ATAC_MODIFIED_LOUVAIN_CLUSTER;
		break;
	case DataField::DataType::Trans:
		cluster_name = METADATA_TRANS_MODIFIED_LOUVAIN_CLUSTER;
		break;
	default:
		break;
	}

	single_cell_multiome->metadata()->mat_.update(cluster_name, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void DataFieldItem::s_receive_slm(std::vector<int> cluster) {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	QString cluster_name;

	auto type = this->data()->data_type_;
	switch (type)
	{
	case DataField::DataType::Rna:
		cluster_name = METADATA_RNA_SLM_CLUSTER;
		break;
	case DataField::DataType::Atac:
		cluster_name = METADATA_ATAC_SLM_CLUSTER;
		break;
	case DataField::DataType::Trans:
		cluster_name = METADATA_TRANS_SLM_CLUSTER;
		break;
	default:
		break;
	}

	single_cell_multiome->metadata()->mat_.update(cluster_name, custom::cast<QVector>(cluster), CustomMatrix::DataType::IntegerFactor);
};

void DataFieldItem::s_smart_local_moving_custom() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto pca = this->data()->pca();
	auto harmony = this->data()->harmony();

	QStringList from_list;
	if (pca) {
		from_list << "PCA";
	}
	if (harmony) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No valid data.");
		return;
	}

	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Louvain Clustering Settings",
		{ 
			"Resolution:0.6",
			"Random State:" + QString::number(single_cell_multiome->random_state_), 
			"Number of Neighbors:30", "Number of Trees:50", 
			"Metric", 
			"From", 
			"Number of Start:10", 
			"Number of Iterations:10", 
			"Modularity Function", 
			"Start Dimension:1",
			"End Dimension:50" 
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
			soap::InputStyle::ComboBox, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit
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
	int start_dim = settings[9].toInt();
	int end_dim = settings[10].toInt();

	int ncol = from_matrix->cols();
	if (start_dim < 0 || end_dim > ncol || end_dim <= start_dim) {
		G_WARN("Illegal dimension selection");
		G_UNLOCK;
		return;
	}

	G_LOG("Start Clustering...");

	SlmWorker* worker;
	if (start_dim == 1 && end_dim == ncol) {
		worker = new SlmWorker(*from_matrix, "SLM", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	}
	else {
		worker = new SlmWorker(from_matrix->operator()(Eigen::all, Eigen::seqN(start_dim - 1, end_dim - start_dim + 1)), "SLM", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	}

	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, DataFieldItem, s_receive_slm);
};

void DataFieldItem::s_modified_louvain_custom() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto pca = this->data()->pca();
	auto harmony = this->data()->harmony();

	QStringList from_list;
	if (pca) {
		from_list << "PCA";
	}
	if (harmony) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No valid data.");
		return;
	}

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Louvain Clustering Settings",
		{ 
			"Resolution:0.6", 
			"Random State:" + QString::number(single_cell_multiome->random_state_), 
			"Number of Neighbors:30",
			"Number of Trees:50", 
			"Metric", 
			"From", 
			"Number of Start:10", 
			"Number of Iterations:10", 
			"Modularity Function", 
			"Start Dimension:1", 
			"End Dimension:50" 
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
			soap::InputStyle::ComboBox, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit},
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

	int start_dim = settings[9].toInt();
	int end_dim = settings[10].toInt();

	int ncol = from_matrix->cols();
	if (start_dim < 0 || end_dim > ncol || end_dim <= start_dim) {
		G_WARN("Illegal dimension selection");
		G_UNLOCK;
		return;
	}

	G_LOG("Start Clustering...");

	SlmWorker* worker;
	if (start_dim == 1 && end_dim == ncol) {
		worker = new SlmWorker(*from_matrix, "Modified Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	}
	else {
		worker = new SlmWorker(from_matrix->operator()(Eigen::all, Eigen::seqN(start_dim - 1, end_dim - start_dim + 1)), "Modified Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	}

	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, DataFieldItem, s_receive_modified_louvain);
};

void DataFieldItem::s_leiden_custom() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto pca = this->data()->pca();
	auto harmony = this->data()->harmony();

	QStringList from_list;
	if (pca != nullptr) {
		from_list << "PCA";
	}
	if (harmony != nullptr) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No valid data.");
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
			"Number of Trees:50", 
			"From", 
			"Start Dimension:1", 
			"End Dimension:50" 
		},
		{ 
			soap::InputStyle::ComboBox, 
			soap::InputStyle::ComboBox, 
			soap::InputStyle::NumericLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::ComboBox,
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit
		},
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

	int start_dim = settings[6].toInt();
	int end_dim = settings[7].toInt();

	int ncol = from_matrix->cols();
	if (start_dim < 0 || end_dim > ncol || end_dim <= start_dim) {
		G_WARN("Illegal dimension selection");
		G_UNLOCK;
		return;
	}

	G_LOG("Start Clustering...");

	LeidenPartitionWorker* worker;
	if (start_dim == 1 && end_dim == ncol) {
		worker = new LeidenPartitionWorker(*from_matrix, method, metric, n_neighbors, n_trees, resolution);
	}
	else {
		worker = new LeidenPartitionWorker(
			from_matrix->operator()(Eigen::all, Eigen::seqN(start_dim - 1, end_dim - start_dim + 1)),
			method, 
			metric, 
			n_neighbors, 
			n_trees, 
			resolution);
	}

	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, DataFieldItem, s_receive_leiden);
}

void DataFieldItem::s_louvain_custom() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto pca = this->data()->pca();
	auto harmony = this->data()->harmony();

	QStringList from_list;
	if (pca) {
		from_list << "PCA";
	}
	if (harmony) {
		from_list << "Harmony";
	}

	if (from_list.isEmpty()) {
		G_WARN("No valid data.");
		return;
	}

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Louvain Clustering Settings",
		{ 
			"Resolution:0.6", 
			"Random State:" + QString::number(single_cell_multiome->random_state_), 
			"Number of Neighbors:30", 
			"Number of Trees:50", 
			"Metric", 
			"From", 
			"Number of Start:10",
			"Number of Iterations:10",
			"Modularity Function", 
			"Start Dimension:1", 
			"End Dimension:50" 
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
			soap::InputStyle::ComboBox, 
			soap::InputStyle::IntegerLineEdit, 
			soap::InputStyle::IntegerLineEdit
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
		from_matrix = &this->pca()->data()->data_.mat_;
	}
	else {
		from_matrix = &this->harmony()->data()->data_.mat_;
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

	int start_dim = settings[9].toInt();
	int end_dim = settings[10].toInt();

	int ncol = from_matrix->cols();
	if (start_dim < 0 || end_dim > ncol || end_dim <= start_dim) {
		G_WARN("Illegal dimension selection");
		G_UNLOCK;
		return;
	}

	G_LOG("Start Clustering...");

	SlmWorker* worker;
	if (start_dim == 1 && end_dim == ncol) {
		worker = new SlmWorker(*from_matrix, "Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	}
	else {
		worker = new SlmWorker(from_matrix->operator()(Eigen::all, Eigen::seqN(start_dim - 1, end_dim - start_dim + 1)), "Louvain", metric, modularity_function, n_start, n_iteration, random_state, n_neighbors, n_trees, resolution);
	}
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, DataFieldItem, s_receive_louvain);
};

void DataFieldItem::s_leiden_default() {

	auto pca = this->data()->pca();

	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}

	G_GETLOCK;

	G_LOG("Start Clustering...");

	LeidenPartitionWorker* worker = new LeidenPartitionWorker(
		pca->data_.mat_,
		"Modularity",
		"Euclidean",
		30,
		50,
		0.6
	);
	G_LINK_WORKER_THREAD(LeidenPartitionWorker, x_leiden_ready, DataFieldItem, s_receive_leiden);
};

void DataFieldItem::s_louvain_default() {

	auto pca = this->data()->pca();

	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}
	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(
		pca->data_.mat_,
		"Louvain",
		"Euclidean",
		1,
		10,
		10,
		1997,
		30,
		50,
		0.6
	);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, DataFieldItem, s_receive_louvain);
};

void DataFieldItem::s_modified_louvain_default() {

	auto pca = this->data()->pca();

	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}
	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(
		pca->data_.mat_,
		"Modified Louvain",
		"Euclidean",
		1,
		10,
		10,
		1997,
		30,
		50,
		0.6
	);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, DataFieldItem, s_receive_modified_louvain);
};

void DataFieldItem::s_smart_local_moving_default() {

	auto pca = this->data()->pca();

	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		return;
	}
	G_GETLOCK;

	G_LOG("Start Clustering...");

	SlmWorker* worker = new SlmWorker(
		pca->data_.mat_,
		"SLM",
		"Euclidean",
		1,
		10,
		10,
		1997,
		30,
		50,
		0.6
	);
	G_LINK_WORKER_THREAD(SlmWorker, x_cluster_ready, DataFieldItem, s_receive_slm);
};

void DataFieldItem::s_harmony() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto pca = this->data()->pca();
	if (pca == nullptr) {
		G_WARN("PCA has not been conducted!");
		G_UNLOCK;
		return;
	}

	QStringList factor_names = single_cell_multiome->metadata()->mat_.get_factor_name(false);

	if (factor_names.isEmpty()) {
		G_LOG("No suitable feature for Harmony");
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Harmony Settings",
		{ "Batch Factor" , "Dimension Start:1", "Dimension End:50"},
		{soap::InputStyle::SimpleChoice, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit },
		{factor_names}
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
		factor_list << single_cell_multiome->metadata()->mat_.get_qstring(factor_name);
	}

	int dim_start = settings[1].toInt();
	int dim_end = settings[2].toInt();

	int ndim = pca->data_.mat_.cols();

	if (dim_start < 1 || dim_start > ndim || dim_end < 1 || dim_end > ndim || dim_start >= dim_end) {
		G_WARN("Illegal Dimension Setting.");
		G_UNLOCK;
		return;
	}

	if (auto item = this->harmony()) {
		item->__remove_this();
	}

	HarmonyWorker* worker = new HarmonyWorker(
		pca->data_.mat_.block(0, dim_start - 1, pca->data_.mat_.rows(), dim_end - dim_start + 1),
		factor_list
	);

	G_LOG("Harmony start...");

	G_LINK_WORKER_THREAD(HarmonyWorker, x_harmony_ready, DataFieldItem, s_receive_harmony)
};

void DataFieldItem::s_svd_custom() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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

	SvdWorker* worker = new SvdWorker(normalized->mat_, perc, 50, single_cell_multiome->random_state_);
	G_LINK_WORKER_THREAD(SvdWorker, x_svd_ready, DataFieldItem, s_receive_svd);
};

void DataFieldItem::s_svd_default() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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

	SvdWorker* worker = new SvdWorker(normalized->mat_, 25, 50, single_cell_multiome->random_state_);
	G_LINK_WORKER_THREAD(SvdWorker, x_svd_ready, DataFieldItem, s_receive_svd);
}

void DataFieldItem::s_pca_default() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Counts data missed!");
		G_UNLOCK;
		return;
	}

	if (auto item = this->pca()) {
		item->__remove_this();
	}

	PcaWorker* worker = new PcaWorker(&counts->mat_, 2000, 50, single_cell_multiome->random_state_);
	G_LINK_WORKER_THREAD(PcaWorker, x_pca_ready, DataFieldItem, s_receive_pca);
};

void DataFieldItem::s_pca_custom() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Counts data missed!");
		G_UNLOCK;
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"PCA settings",
		{
			"Number of dimension:50",
			"Random State:" + QString::number(single_cell_multiome->random_state_),
			"Variable Feature Proportion:0.1",
			"Fixed Variable Feature Number?:yes",
			"Fixed Variable Feature Number:2000"
		},
		{
		soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::IntegerLineEdit,
			soap::InputStyle::NumericLineEdit,
			soap::InputStyle::SwitchButton,
			soap::InputStyle::IntegerLineEdit
	}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	int n_dimension = settings[0].toInt();
	if (n_dimension < 10 || n_dimension > 100) {
		G_WARN("Number of dimension should be between 10 and 100.");
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
			G_WARN("Number of Variable Features shold be more than 1000.");
			G_UNLOCK;
			return;
		}
	}
	else {
		variable_proportion = settings[2].toDouble();
		if (variable_proportion < 0.01 || variable_proportion > 1) {
			G_WARN("Number of Variable Features Proportion should be between 0.01 and 1.");
			G_UNLOCK;
			return;
		}
	}

	if (auto item = this->pca()) {
		item->__remove_this();
	}

	PcaWorker* worker = new PcaWorker(
		&counts->mat_,
		variable_feature_number_fixed ? variable_number : variable_proportion,
		n_dimension,
		random_state
	);
	G_LINK_WORKER_THREAD(PcaWorker, x_pca_ready, DataFieldItem, s_receive_pca)
};

void DataFieldItem::s_filter_by_feature() {
	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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

	auto item = static_cast<SingleCellMultiomeItem*>(this->signal_emitter_->get_item(single_cell_multiome));
	item->col_slice(filter, in_place);

	G_UNLOCK;
};

void DataFieldItem::s_receive_normalize(SparseDouble* data) {
	
	QString title;

	switch (this->data()->data_type_)
	{
	case DataField::DataType::Rna:
		title = VARIABLE_RNA_NORMALIZED;
		break;
	case DataField::DataType::Atac:
		title = VARIABLE_ATAC_NORMALIZED;
		break;
	case DataField::DataType::Trans:
		title = VARIABLE_TRANS_NORMALIZED;
		break;
	default:
		break;
	}

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
};

void DataFieldItem::s_receive_harmony(Eigen::MatrixXd mat) {

	QString harmony_name;

	if (this->data()->data_type_ == DataField::DataType::Rna) {
		harmony_name = VARIABLE_RNA_HARMONY;
	}
	else if (this->data()->data_type_ == DataField::DataType::Atac) {
		harmony_name = VARIABLE_ATAC_HARMONY;
	}
	else if (this->data()->data_type_ == DataField::DataType::Trans) {
		harmony_name = VARIABLE_TRANS_HARMONY;
	}

	QString title = this->signal_emitter_->get_unique_name(harmony_name);

	DATA_SUBMODULES(Embedding)[title] =
		Embedding(
			Embedding::DataType::Harmony,
			mat,
			this->data()->counts()->colnames_,
			custom::paste(harmony_name, custom::cast<QString>(custom::seq_n(1, mat.cols())), "-")
		);

	auto item = new EmbeddingItem(
		title,
		this->index_tree_,
		this->data()->harmony(),
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);
	this->set_item(item);
};

void DataFieldItem::s_tsne_default() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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

	TsneWorker* worker = new TsneWorker(tsne_input, single_cell_multiome->random_state_);

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, DataFieldItem, s_receive_tsne)
};

void DataFieldItem::s_tsne_custom() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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
			"Random State:" + QString::number(single_cell_multiome->random_state_)
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

	G_LINK_WORKER_THREAD(TsneWorker, x_tsne_ready, DataFieldItem, s_receive_tsne)
};

void DataFieldItem::s_receive_tsne(Eigen::MatrixXd mat) {

	QString tsne_name;

	if (this->data()->data_type_ == DataField::DataType::Rna) {
		tsne_name = VARIABLE_RNA_TSNE;
	}
	else if (this->data()->data_type_ == DataField::DataType::Atac) {
		tsne_name = VARIABLE_ATAC_TSNE;
	}
	else if (this->data()->data_type_ == DataField::DataType::Trans) {
		tsne_name = VARIABLE_TRANS_TSNE;
	}

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

void DataFieldItem::s_umap_default() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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

	const int random_state = single_cell_multiome->random_state_;
	G_LOG("Default RNA UMAP start, using 1~" + QString::number(last_dimension + 1) + " PCA components...");

	if (this->data()->data_type_ == DataField::DataType::Atac) {
		G_NOTICE("If you use SVD before, it is recommended to exclude the first PCA dimension in UMAP because of its high correlation with sequencing depth.");
	}

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
	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, DataFieldItem, s_receive_umap)
};

void DataFieldItem::s_umap_custom() {

	G_GETLOCK;

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

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
			"Random state:" + QString::number(single_cell_multiome->random_state_), 
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
	G_LINK_WORKER_THREAD(UmapWorker, x_umap_ready, DataFieldItem, s_receive_umap)
};

void DataFieldItem::s_tfidf_default() {

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
	G_LINK_WORKER_THREAD(TfidfWorker, x_tfidf_ready, DataFieldItem, s_receive_tfidf);
};

void DataFieldItem::s_receive_tfidf(SparseDouble* data) {


	QString title;

	switch (this->data()->data_type_)
	{
	case DataField::DataType::Rna:
		title = VARIABLE_RNA_NORMALIZED;
		break;
	case DataField::DataType::Atac:
		title = VARIABLE_ATAC_NORMALIZED;
		break;
	case DataField::DataType::Trans:
		title = VARIABLE_TRANS_NORMALIZED;
		break;
	default:
		break;
	}

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

void DataFieldItem::s_log_normalize_default() {

	G_GETLOCK;

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Counts data missed.");
		G_UNLOCK;
		return;
	}

	if (auto item = this->normalized()) {
		item->__remove_this();
	}

	G_LOG("Log normalize start...");
	LogNormalizeWorker* worker = new LogNormalizeWorker(
		counts,
		10000.0
	);
	G_LINK_WORKER_THREAD(LogNormalizeWorker, x_log_normalize_ready, DataFieldItem, s_receive_normalize)

};

void DataFieldItem::s_log_normalize_custom() {

	G_GETLOCK;

	auto counts = this->data()->counts();
	if (counts == nullptr) {
		G_WARN("Counts data missed.");
		G_UNLOCK;
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set the Scale Factor",
		{ "Scale Factor(>0):10000" },
		{ soap::InputStyle::IntegerLineEdit}
	);
	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	double scale_factor = settings[0].toDouble();
	if (scale_factor <= 0.0) {
		G_NOTICE("Scale factor must be positive.");
		G_UNLOCK;
		return;
	}

	if (auto item = this->normalized()) {
		item->__remove_this();
	}

	G_LOG("Log normalize start...");
	LogNormalizeWorker* worker = new LogNormalizeWorker(
		counts,
		scale_factor
	);
	G_LINK_WORKER_THREAD(LogNormalizeWorker, x_log_normalize_ready, DataFieldItem, s_receive_normalize)
};

void DataFieldItem::atac_view_quality() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	const auto& metadata = single_cell_multiome->metadata()->mat_;

	if (!metadata.contains(METADATA_ATAC_UMI_NUMBER, CustomMatrix::DataType::IntegerNumeric)
		|| !metadata.contains(METADATA_ATAC_UNIQUE_PEAK_NUMBER, CustomMatrix::DataType::IntegerNumeric)
		)
	{
		G_WARN("Quality information has been lost.");
		return;
	}

	auto trans = metadata.get_double(METADATA_ATAC_UMI_NUMBER);
	int control_point_number = 24;

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, layout] = custom_plot::prepare_lg(gs);
	
	auto colors = gs.palette({ METADATA_ATAC_UMI_NUMBER , METADATA_ATAC_UNIQUE_PEAK_NUMBER });

	QCPAxisRect* axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 0, axis_rect);

	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[0], "", "UMI Count");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[0], "", "UMI Count");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_ATAC_UNIQUE_PEAK_NUMBER);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 1, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[1], "", "Peak Count");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[1], "", "Peak Count");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	custom_plot::add_title(draw_area, "ATAC Data Quality", gs);

	this->draw_suite_->update(draw_area);
}

void DataFieldItem::rna_view_quality() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		G_UNLOCK;
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	const auto& metadata = single_cell_multiome->metadata()->mat_;

	if (!(metadata.contains(METADATA_RNA_UMI_NUMBER) &&
		metadata.contains(METADATA_RNA_UNIQUE_GENE_NUMBER) &&
		metadata.contains(METADATA_RNA_MITOCHONDRIAL_CONTENT))) {
		G_WARN("Quality information has lost.");
		return;
	}
	if (!(metadata.data_type_.at(METADATA_RNA_UMI_NUMBER) == CustomMatrix::DataType::IntegerNumeric &&
		metadata.data_type_.at(METADATA_RNA_UNIQUE_GENE_NUMBER) == CustomMatrix::DataType::IntegerNumeric &&
		metadata.data_type_.at(METADATA_RNA_MITOCHONDRIAL_CONTENT) == CustomMatrix::DataType::DoubleNumeric)) {
		G_WARN("Quality information format has been overwritten.");
		return;
	}

	auto colors = custom_plot::utility::kmeans_palette(3);

	auto trans = metadata.get_double(METADATA_RNA_UMI_NUMBER);
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

	trans = metadata.get_double(METADATA_RNA_UNIQUE_GENE_NUMBER);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 1, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[1], "", "Gene Count");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[1], "", "Gene Count");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	trans = metadata.get_double(METADATA_RNA_MITOCHONDRIAL_CONTENT);
	axis_rect = custom_plot::patch::new_axis_rect(draw_area);
	layout->addElement(0, 2, axis_rect);
	if (gs.use_boxplot()) {
		custom_plot::patch::single_box_plot(draw_area, axis_rect, trans, colors[2], "", "Mitochondrial Count");
	}
	else {
		custom_plot::patch::single_violin_plot(draw_area, axis_rect, trans, control_point_number, colors[2], "", "Mitochondrial Count");
	}
	custom_plot::patch::clear_bottom_axis(axis_rect);

	custom_plot::add_title(draw_area, "RNA Data Quality", gs);

	this->draw_suite_->update(draw_area);
}

void DataFieldItem::s_view_quality() {

	auto type = this->data()->data_type_;

	switch (type)
	{
	case DataField::DataType::Rna:
		this->rna_view_quality();
		break;
	case DataField::DataType::Atac:
		this->atac_view_quality();
		break;
	default:
		break;
	}
};

void DataFieldItem::s_correlation_heatmap_plot() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	const auto& metadata = single_cell_multiome->metadata()->mat_;

	auto factor_info = metadata.get_factor_information(false);

	if (factor_info.isEmpty()) {
		G_WARN("No Factor has more than one level.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Correlation Heatmap Plot Settings",
		{"Normalize:yes", "Choose Factor"},
		{soap::InputStyle::SwitchButton, soap::InputStyle::FactorChoice},
		{},
		{factor_info}
	);

	if (settings.isEmpty()) {
		return;
	}

	bool normalize = switch_to_bool(settings[0]);

	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);

	if (levels.isEmpty()) {
		return;
	}
	int n_level = levels.size();
	QStringList factor = metadata.get_qstring(factor_name);
	int n_cell = factor.size();

	Eigen::MatrixXd data(n_cell, n_level);
	if (normalize) {
		auto normalized = this->data()->normalized();

		if (normalized == nullptr) {
			G_WARN("No Normalized Data.");
			return;
		}

		auto&& mat = normalized->mat_;

		for (int i = 0; i < n_level; ++i) {
			data.col(i) = custom::column_slice_and_row_mean(mat, custom::equal(factor, levels[i]));
		}
	}
	else {
		auto counts = this->data()->counts();

		if (counts == nullptr) {
			G_WARN("No Counts Data.");
			return;
		}

		auto&& mat = counts->mat_;

		for (int i = 0; i < n_level; ++i) {
			data.col(i) = custom::column_slice_and_row_mean(mat, custom::equal(factor, levels[i]));
		}
	}

	auto cor = custom::cor(data);

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	custom_plot::heatmap_plot2(
		draw_area,
		axis_rect,
		10,
		10,
		2,
		levels,
		levels,
		cor,
		gs
	);

	custom_plot::add_gradient_legend(
		draw_area,
		legend_layout,
		cor.minCoeff(),
		cor.maxCoeff(),
		"Correlation",
		gs);

	this->draw_suite_->update(draw_area);
};

void DataFieldItem::s_distribution_plot() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	const auto& metadata = single_cell_multiome->metadata()->mat_;

	auto factor_info = metadata.get_factor_information(false);

	if (factor_info.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Feature",
		{ "Feature", "Group", "Normalized", "Show Value" },
		{ soap::InputStyle::MultipleLineEdit, soap::InputStyle::FactorChoice,
		soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton},
		{},
		{ factor_info }
	);

	if (settings.isEmpty())return;

	bool normalized = switch_to_bool(settings[2]);

	bool show_value = switch_to_bool(settings[3]);
	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);
	QStringList factors = metadata.get_qstring(factor_name);
	
	if (levels.isEmpty()) {
		return;
	}

	auto features = multiple_line_edit_to_list(settings[0]);
	if (features.isEmpty()) return;

	FeatureHandler handler(single_cell_multiome);

	auto data = custom::sapply(
		features,
		[&handler, normalized, this](auto&& s) {
		return handler.get_data(QUERY_INFO{ s, normalized, this->data()->data_type_ == DataField::DataType::Trans });
	}
	);

	auto filter = custom::sapply(data, [](auto&& d) {return d.is_continuous(); });
	data = custom::sliced(data, filter);

	if (data.isEmpty()) {
		G_WARN("No valid data.");
		return;
	}

	const auto& gs = this->draw_suite_->graph_settings_;
	QCustomPlot* draw_area = custom_plot::initialize_plot(gs);
	auto colors = gs.palette(levels);

	QCPMarginGroup* margin_group = new QCPMarginGroup(draw_area);

	int n_valid_feature = data.size();
	int group_number = levels.size();

	for (int i = 0; i < n_valid_feature; ++i) {

		auto feature_data = data[i].get_continuous();

		QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);
		axis_rect->setMarginGroup(QCP::msLeft, margin_group);
		draw_area->plotLayout()->addElement(i, 0, axis_rect);

		auto [min, max] = custom_plot::patch::violin_batch(
			draw_area,
			axis_rect,
			factors,
			levels,
			colors,
			custom::cast<Eigen::ArrayX>(feature_data),
			1.0,
			2.0,
			16
		);

		custom_plot::patch::set_range(axis_rect, QCPRange(0, 2 * group_number), custom_plot::utility::get_range(min, max));
		custom_plot::patch::remove_bottom_axis(axis_rect);
		custom_plot::patch::clear_left_axis(axis_rect);
		axis_rect->axis(QCPAxis::atLeft)->setBasePen(QPen(Qt::black, 3));
		custom_plot::set_left_title(axis_rect, data[i].name, gs);
		if (show_value) {
			axis_rect->axis(QCPAxis::atLeft)->setTicks(true);
			axis_rect->axis(QCPAxis::atLeft)->setTickLabels(true);
		}
		axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);
		draw_area->plotLayout()->setRowStretchFactor(i, 1.0);
	}

	QCPAxisRect* axis_rect = new QCPAxisRect(draw_area, true);
	custom_plot::patch::set_fixed_height(axis_rect, 1);
	axis_rect->setMarginGroup(QCP::msLeft, margin_group);
	draw_area->plotLayout()->addElement(n_valid_feature, 0, axis_rect);
	draw_area->plotLayout()->setRowStretchFactor(n_valid_feature, 0.1);
	axis_rect->axis(QCPAxis::atBottom)->setRange(QCPRange(0, 2 * group_number));
	custom_plot::patch::clear_bottom_axis(axis_rect);
	custom_plot::patch::remove_left_axis(axis_rect);
	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(QPen(Qt::black, 3));
	custom_plot::set_bottom_axis_label(
		axis_rect,
		Eigen::ArrayXd::LinSpaced(group_number, 1, 2 * group_number - 1),
		levels,
		6,
		gs
	);

	this->draw_suite_->update(draw_area);
};

void DataFieldItem::s_cell_heatmap_plot() {
	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto& metadata = single_cell_multiome->metadata()->mat_;
	QMap<QString, QStringList> map = metadata.get_factor_information();

	if (map.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Heatmap Plot",
		{ "Feature Names", "Group(Ordered By):Group", "Normalized:yes", "Groups(Annotation)",
		"Show Row Name:no", "Show Column Name:no", "Show Annotation Name:no","Show Annotation Legend:yes"},
		{ soap::InputStyle::MultipleLineEdit, soap::InputStyle::FactorChoice,
		soap::InputStyle::SwitchButton, soap::InputStyle::SimpleChoice, soap::InputStyle::SwitchButton,
		soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton},
		{ map.keys()}, {map}
	);

	if (settings.isEmpty()) {
		return;
	}

	auto feature_names = multiple_line_edit_to_list(settings[0]);
	if (feature_names.isEmpty())return;

	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);
	if (levels.isEmpty())return;

	bool normalized = switch_to_bool(settings[2]);
	auto annotation_names = simple_choice_to_list(settings[3]);

	HEATMAP_PLOT_ELEMENTS ele;

	ele.show_row_names = switch_to_bool(settings[4]);
	ele.show_column_names = switch_to_bool(settings[5]);
	ele.show_column_annotation_name = switch_to_bool(settings[6]);
	ele.show_annotation_legend = switch_to_bool(settings[7]);

	QStringList factor_data = metadata.get_qstring(factor_name);

	QVector<int> index;
	for (auto&& level : levels) {
		auto sub_index = custom::match(factor_data, level);
		index << sub_index;
	}

	int n_cell_use = index.size();

	if (normalized) {

		auto normalized = this->data()->normalized();
		if (normalized == nullptr) {
			G_WARN("Data has not been normalized.");
			return;
		}

		auto filter = custom::in(feature_names, normalized->rownames_);
		if (filter.count() == 0) {
			G_NOTICE("No Feature Found in data.");
			return;
		}

		feature_names = custom::sliced(feature_names, filter);

		auto feature_index = custom::index_of(feature_names, normalized->rownames_);

		int n_feature = feature_names.size();


		Eigen::MatrixXd values(n_feature, n_cell_use);
		for (int i = 0; i < n_feature; ++i) {
			Eigen::ArrayXd exp = normalized->mat_.row(feature_index[i]);
			values.row(i) = custom::reordered(exp, index);
		}

		const auto& gs = this->draw_suite_->graph_settings_;

		ele.mat = values;
		ele.legend_title = "Expression";
		ele.column_names = custom::reordered(normalized->colnames_, index);
		ele.row_names = feature_names;
		for (auto&& annotation : annotation_names) {
			ele.column_annotations << qMakePair(annotation, custom::reordered(metadata.get_qstring(annotation), index));
		}
		
		auto draw_area = custom_plot::heatmap_plot(ele, gs);

		this->draw_suite_->update(draw_area);
	}
	else {
		auto counts = this->data()->counts();
		if (counts == nullptr) {
			G_WARN("No Counts Data.");
			return;
		}

		auto filter = custom::in(feature_names, counts->rownames_);
		if (filter.count() == 0) {
			G_NOTICE("No Feature Found in data.");
			return;
		}

		feature_names = custom::sliced(feature_names, filter);

		auto feature_index = custom::index_of(feature_names, counts->rownames_);

		int n_feature = feature_names.size();


		Eigen::MatrixXd values(n_feature, n_cell_use);
		for (int i = 0; i < n_feature; ++i) {
			Eigen::ArrayXd exp = counts->mat_.row(feature_index[i]).cast<double>();
			values.row(i) = custom::reordered(exp, index);
		}

		const auto& gs = this->draw_suite_->graph_settings_;

		ele.mat = values;
		ele.legend_title = "Expression";
		ele.column_names = custom::reordered(counts->colnames_, index);
		ele.row_names = feature_names;
		for (auto&& annotation : annotation_names) {
			ele.column_annotations << qMakePair(annotation, custom::reordered(metadata.get_qstring(annotation), index));
		}

		auto draw_area = custom_plot::heatmap_plot(ele, gs);

		this->draw_suite_->update(draw_area);
	}

};

void DataFieldItem::s_heatmap_plot() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto& metadata = single_cell_multiome->metadata()->mat_;
	QMap<QString, QStringList> map = metadata.get_factor_information();

	if (map.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Heatmap Plot",
		{ "Feature Names", "Group:Group", "Normalized:yes", "Cell Width:10",
		"Cell Height:10", "Border Width:0", "Annotation at top:no", "Normalize by row:yes" },
		{ soap::InputStyle::MultipleLineEdit, soap::InputStyle::FactorChoice,
		soap::InputStyle::SwitchButton, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton	},
		{}, { map}
	);

	if (settings.isEmpty())return;

	auto feature_names = multiple_line_edit_to_list(settings[0]);
	if (feature_names.isEmpty())return;

	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);
	if (levels.isEmpty())return;

	int n_level = levels.size();
	QStringList factor = metadata.get_qstring(factor_name);
	bool is_normalized = switch_to_bool(settings[2]);

	int cell_width = settings[3].toInt();
	if (cell_width <= 0) {
		G_WARN("Cell Width must be positive.");
		return;
	}

	int cell_height = settings[4].toInt();
	if (cell_height <= 0) {
		G_WARN("Cell Height must be positive.");
		return;
	}

	int border_width = settings[5].toInt();
	if (border_width < 0) {
		G_WARN("Border Width must not be negative.");
		return;
	}

	bool annotation_at_top = switch_to_bool(settings[6]);
	bool normalize_by_row = switch_to_bool(settings[7]);

	if (is_normalized) {

		auto normalized = this->data()->normalized();
		if(normalized == nullptr) {
			G_WARN("Data has not been normalized.");
			return;
		}

		auto filter = custom::in(feature_names, normalized->rownames_);
		if (filter.count() == 0) {
			G_NOTICE("No Feature Found in data.");
			return;
		}

		feature_names = custom::sliced(feature_names, filter);

		auto feature_index = custom::index_of(feature_names, normalized->rownames_);

		int n_feature = feature_names.size();

		Eigen::MatrixXd values(n_feature, n_level);
		for (int i = 0; i < n_level; ++i) {
			filter = custom::equal(factor, levels[i]);
			for (int j = 0; j < n_feature; ++j) {
				Eigen::ArrayXd exp = normalized->mat_.row(feature_index[j]);
				exp = custom::sliced(exp, filter);
				values(j, i) = exp.mean();
			}
		}

		if (normalize_by_row) {
			for (int j = 0; j < n_feature; ++j) {
				auto [min_exp, max_exp] = std::ranges::minmax(values.row(j));
				if (min_exp == max_exp) {
					values.row(j).setConstant(0.0);
					continue;
				}
				else {
					double mean = values.row(j).mean();
					double std = custom::sd(values.row(j));
					values.row(j) = (values.row(j).array() - mean) / std;
				}
			}
		}

		const auto& gs = this->draw_suite_->graph_settings_;
		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::heatmap_plot(
			draw_area,
			axis_rect,
			legend_layout,
			cell_width,
			cell_height,
			border_width,
			feature_names,
			levels,
			values,
			normalize_by_row,
			gs.get_legend_title(factor_name, 1),
			annotation_at_top,
			gs);

		this->draw_suite_->update(draw_area);
	}
	else {

		auto counts = this->data()->counts();
		if(counts == nullptr) {
			G_WARN("Count data is missed.");
			return;
		}

		auto filter = custom::in(feature_names, counts->rownames_);
		if (filter.count() == 0) {
			G_NOTICE("No Feature Found in data.");
			return;
		}

		feature_names = custom::sliced(feature_names, filter);

		auto feature_index = custom::index_of(feature_names, counts->rownames_);

		int n_feature = feature_names.size();

		Eigen::MatrixXd values(n_level, n_feature);

		for (int i = 0; i < n_level; ++i) {
			filter = custom::equal(factor, levels[i]);
			for (int j = 0; j < n_feature; ++j) {
				Eigen::ArrayXi exp = counts->mat_.row(feature_index[j]);
				exp = custom::sliced(exp, filter);
				values(i, j) = exp.cast<double>().mean();
			}
		}

		if (normalize_by_row) {
			for (int j = 0; j < n_feature; ++j) {
				auto [min_exp, max_exp] = std::ranges::minmax(values.row(j));
				if (min_exp == max_exp) {
					values.row(j).setConstant(0.0);
					continue;
				}
				else {
					double mean = values.row(j).mean();
					double std = custom::sd(values.row(j));
					values.row(j) = (values.row(j).array() - mean) / std;
				}
			}
		}

		const auto& gs = this->draw_suite_->graph_settings_;
		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::heatmap_plot(
			draw_area,
			axis_rect,
			legend_layout,
			cell_width,
			cell_height,
			border_width,
			levels,
			feature_names,
			values,
			normalize_by_row,
			gs.get_legend_title(factor_name, 1),
			annotation_at_top,
			gs);

		this->draw_suite_->update(draw_area);

	}
};

void DataFieldItem::s_bubble_plot() {

	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This Field is not a part of Single Cell Multiome Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto& metadata = single_cell_multiome->metadata()->mat_;
	QMap<QString, QStringList> map = metadata.get_factor_information();

	if (map.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Bubble Plot",
		{ "Feature Names", "Group:Group", "Normalized" },
		{ soap::InputStyle::MultipleLineEdit, soap::InputStyle::FactorChoice, soap::InputStyle::SwitchButton},
		{}, { map}
	);
	if (settings.isEmpty())return;

	auto feature_names = multiple_line_edit_to_list(settings[0]);
	if (feature_names.isEmpty())return;

	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);
	if (levels.isEmpty())return;
	int n_level = levels.size();
	QStringList factor = metadata.get_qstring(factor_name);
	bool is_normalized = switch_to_bool(settings[2]);


	if (is_normalized) {

		auto normalized = this->data()->normalized();
		if (normalized == nullptr) {
			G_WARN("Data has not been normalized.");
			return;
		}

		auto filter = custom::in(feature_names, normalized->rownames_);
		if (filter.count() == 0) {
			G_NOTICE("No Feature Found in data.");
			return;
		}

		feature_names = custom::sliced(feature_names, filter);

		auto feature_index = custom::index_of(feature_names, normalized->rownames_);

		int n_feature = feature_names.size();

		Eigen::MatrixXd values(n_level, n_feature);
		Eigen::MatrixXi proportion(n_level, n_feature);
		for (int i = 0; i < n_level; ++i) {
			filter = custom::equal(factor, levels[i]);
			for (int j = 0; j < n_feature; ++j) {
				Eigen::ArrayXd exp = normalized->mat_.row(feature_index[j]);
				exp = custom::sliced(exp, filter);
				values(i, j) = exp.mean();
				proportion(i, j) = ceil((exp > 0).count() / (double)exp.size() * 50) + 1;
			}
		}

		const auto& gs = this->draw_suite_->graph_settings_;
		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::bubble_plot(
			draw_area,
			axis_rect,
			legend_layout,
			levels,
			feature_names,
			values,
			proportion,
			gs.get_legend_title("Expression"),
			gs.get_legend_title("Proportion", 1),
			gs);

		this->draw_suite_->update(draw_area);
	}
	else {

		auto counts = this->data()->counts();
		if (counts == nullptr) {
			G_WARN("Count data is missed.");
			return;
		}

		auto filter = custom::in(feature_names, counts->rownames_);
		if (filter.count() == 0) {
			G_NOTICE("No Feature Found in data.");
			return;
		}

		feature_names = custom::sliced(feature_names, filter);

		auto feature_index = custom::index_of(feature_names, counts->rownames_);

		int n_feature = feature_names.size();

		Eigen::MatrixXd values(n_level, n_feature);
		Eigen::MatrixXi proportion(n_level, n_feature);

		for (int i = 0; i < n_level; ++i) {
			filter = custom::equal(factor, levels[i]);
			for (int j = 0; j < n_feature; ++j) {
				Eigen::ArrayXi exp = counts->mat_.row(feature_index[j]);
				exp = custom::sliced(exp, filter);
				values(i, j) = exp.cast<double>().mean();
				proportion(i, j) = ceil((exp > 0).count() / (double)exp.size() * 50) + 1;
			}
		}

		const auto& gs = this->draw_suite_->graph_settings_;
		auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

		custom_plot::bubble_plot(
			draw_area,
			axis_rect,
			legend_layout,
			custom::reversed(levels),
			feature_names,
			values,
			proportion,
			gs.get_legend_title("Expression"),
			gs.get_legend_title("Proportion", 1),
			gs);

		this->draw_suite_->update(draw_area);

	}
};

bool DataFieldItem::atac_re_estimate_quality_parameters() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		return false;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto atac_counts = this->data()->counts();

	if (atac_counts == nullptr) {
		return false;
	}

	const int n_cell = atac_counts->cols();

	Eigen::ArrayXi peak_count = custom::col_sum_mt(atac_counts->mat_);
	Eigen::ArrayXi peak_gene(n_cell);

	for (std::size_t i = 0; i < n_cell; ++i) {
		peak_gene[i] = atac_counts->mat_.outerIndexPtr()[i + 1] - atac_counts->mat_.outerIndexPtr()[i];
	}

	auto metadata = single_cell_multiome->metadata();

	metadata->mat_.update(METADATA_ATAC_UMI_NUMBER, custom::cast<QVector>(peak_count));
	metadata->mat_.update(METADATA_ATAC_UNIQUE_PEAK_NUMBER, custom::cast<QVector>(peak_gene));

	return true;
};

bool DataFieldItem::rna_re_estimate_quality_parameters() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		return false;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();
	auto species = single_cell_multiome->species_;

	auto rna_counts = this->data()->counts();

	if (rna_counts == nullptr) {
		return false;
	}

	const int n_cell = rna_counts->cols();

	Eigen::ArrayXi col_count = custom::col_sum_mt(rna_counts->mat_);
	Eigen::ArrayXi col_gene(n_cell);
	Eigen::ArrayXd mitochondrial_content = Eigen::ArrayXd::Zero(n_cell);
	Eigen::ArrayXd ribosomal_content = Eigen::ArrayXd::Zero(n_cell);

	for (int i = 0; i < n_cell; ++i) {
		col_gene[i] = rna_counts->mat_.outerIndexPtr()[i + 1] - rna_counts->mat_.outerIndexPtr()[i];
	}
	QList<int> mitochondrial_location, ribosomal_location;
	QStringList gene_symbols = rna_counts->rownames_;
	const int n_gene = gene_symbols.size();
	if (species == soap::Species::Human) {

		for (int i = 0; i < n_gene; ++i) {
			if (gene_symbols.at(i).startsWith("MT-")) {
				mitochondrial_location << i;
			}

			if (gene_symbols.at(i).startsWith("RPS") || gene_symbols.at(i).startsWith("RPL")) {
				ribosomal_location << i;
			}
		}
	}
	else if (species == soap::Species::Mouse) {

		for (int i = 0; i < n_gene; ++i) {
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
		mitochondrial_content /= col_count.cast<double>();
	}
	if (ribosomal_location.length() > 0) {
		for (int& i : ribosomal_location) {
			ribosomal_content += rna_counts->mat_.row(i).cast<double>();
		}
		ribosomal_content /= col_count.cast<double>();
	}

	custom::remove_na(mitochondrial_content);
	custom::remove_na(ribosomal_content);

	auto metadata = single_cell_multiome->metadata();
	metadata->mat_.update(METADATA_RNA_UMI_NUMBER, custom::cast<QVector>(col_count));
	metadata->mat_.update(METADATA_RNA_UNIQUE_GENE_NUMBER, custom::cast<QVector>(col_gene));
	metadata->mat_.update(METADATA_RNA_MITOCHONDRIAL_CONTENT, custom::cast<QVector>(mitochondrial_content));
	metadata->mat_.update(METADATA_RNA_RIBOSOMAL_CONTENT, custom::cast<QVector>(ribosomal_content));

	return true;
};

void DataFieldItem::s_re_estimate_quality_parameters() {

	if (this->data()->data_type_ == DataField::DataType::Rna) {
		if(!this->rna_re_estimate_quality_parameters()){
			G_WARN("Quality Parameter Estimation failed.");
		}
		else {
			G_NOTICE("Quality Parameter Estimation finished.");
		}

		return;
	}

	if (this->data()->data_type_ == DataField::DataType::Atac) {
		if (!this->atac_re_estimate_quality_parameters()) {
			G_WARN("Quality Parameter Estimation failed.");
		}
		else {
			G_NOTICE("Quality Parameter Estimation finished.");
		}

		return;
	}
};