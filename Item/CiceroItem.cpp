#include "CiceroItem.h"

#include "MatrixWindow.h"
#include "CustomPlot.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "DifferentialAnalysisWorker.h"
#include "SimpleUmapWorker.h"
#include "GenomeUtility.h"

void CiceroItem::__set_menu() {
	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

	ADD_MAIN_ACTION("Differential Analysis", s_differential_analysis);

	ADD_MAIN_ACTION("Show Regulation Group", s_show_group);
	ADD_MAIN_ACTION("Show Transcription Factor Possibility", s_show_group_transcription_factor_possibility);

	ADD_MAIN_ACTION("Run UMAP", s_umap);

	ADD_MAIN_ACTION("Show CCAN Coverage", s_show_ccan_coverage);

	ADD_MAIN_ACTION("Scan Downstream Ccan", s_scan_downstream_ccan);

	ADD_MAIN_ACTION("Ccan Embedding Feature Plot", s_ccan_embedding_feature_plot);

	ADD_MAIN_ACTION("Reset CCAN Threshold", s_reset_threshold);

};

void CiceroItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Embedding));
	this->check_variable(DATA_SUBMODULES(DifferentialAnalysis));
};

void CiceroItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->regulation_groups_.size()) + " ]");
};

void CiceroItem::s_show_ccan_coverage() {


	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("No SingleCell Data Found.");
		return;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto fragments = single_cell_multiome->fragments();
	if (fragments == nullptr) {
		G_WARN("No Fragments Loaded.");
		return;
	}

	const auto& metadata = single_cell_multiome->metadata()->mat_;

	auto group_info = metadata.get_factor_information();

	if (group_info.isEmpty()) {
		G_WARN("No Suitable Metadata for Visualization.");
		return;
	}
	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Coverage Plot Settings",
		{ "Group", "ccan ID", "Draw Gene:yes", "Draw Link:yes", "Link Cutoff:0.5", "Draw Legend:no"},
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

	QString ccan_id = settings[1];

	int id = ccan_id.toInt() - 1;

	int n_group = this->data()->regulation_groups_.size();
	if (id < 0 || id >= n_group) {
		G_WARN("Illegal ID.");
		return;
	}

	QStringList peak_names = _Cs reordered(this->data()->peak_names_, this->data()->regulation_groups_[id]);
	QVector<std::pair<int, int>> ccan_locs;

	auto [chr, start, end, success] = _Cs string_to_peak(peak_names[0]);

	if (!success) {
		G_WARN("Illegal peak name");
		return;
	}

	ccan_locs << std::make_pair(start, end);

	for (auto&& p : peak_names) {
		auto [c, s, e, suc] = _Cs string_to_peak(p);
		if (!suc) {
			G_WARN("Illegal peak name");
			return;
		}

		ccan_locs << std::make_pair(s, e);

		if (start > s) {
			start = s;
		}

		if (end < e) {
			end = e;
		}
	}

	start -= 10000;

	if (start < 1) {
		start = 1;
	}

	end += 10000;

	QString region = chr + ":" + QString::number(start) + "-" + QString::number(end);

	bool draw_gene = switch_to_bool(settings[2]);
	bool draw_link = switch_to_bool(settings[3]);
	double link_cutoff = settings[4].toDouble();
	bool draw_legend = switch_to_bool(settings[5]);

	CoveragePlotWorker* worker = new CoveragePlotWorker(
		single_cell_multiome,
		fragments,
		factors,
		levels,
		region, 
		"CCAN " + ccan_id,
		ccan_locs,
		draw_gene,
		draw_link, 
		link_cutoff,
		draw_legend);
	G_LINK_WORKER_THREAD(CoveragePlotWorker, x_plot_ready, CiceroItem, s_receive_coverage_plot_data);
};

void CiceroItem::s_receive_coverage_plot_data(COVERAGE_PLOT_ELEMENTS res) {

	auto& gs = this->draw_suite_->graph_settings_;

	auto draw_area = _Cp coverage_plot(res, gs);

	this->draw_suite_->update(draw_area);
};

void CiceroItem::s_show_group() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Show Regulation Group Content",
		{ "Group ID" },
		{soap::InputStyle::StringLineEdit}
	);

	if (settings.isEmpty()) {
		return;
	}

	int id = settings[0].toInt() - 1;

	int n_group = this->data()->regulation_groups_.size();
	if (id < 0 || id >= n_group) {
		G_WARN("Illegal ID.");
		return;
	}

	QStringList res = _Cs reordered(this->data()->peak_names_, this->data()->regulation_groups_[id]);

	MatrixWindow::show_matrix(&res, "Regulation Group " + QString::number(id + 1), this->signal_emitter_);
};

void CiceroItem::s_show_group_transcription_factor_possibility() {


	if (!this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("Illegal Data.");
		return;
	}

	SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

	auto motif_position = single_cell_multiome->motif_position();

	if (!motif_position) {
		G_WARN("Motif Position is not determined.");
		return;
	}

	if (motif_position->peak_names_ != this->data()->peak_names_) {
		G_WARN("Unmatched peak data.");
		return;
	}

	auto rna_normalized = single_cell_multiome->rna_normalized();
	if (!rna_normalized) {
		G_WARN("Rna Data is not normalized.");
		return;
	}

	auto factor_info = single_cell_multiome->metadata()->mat_.get_factor_information(false);

	if (factor_info.isEmpty()) {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Show Transcription Factor Possibility",
			{ "Group ID" },
			{soap::InputStyle::StringLineEdit}
		);

		if (settings.isEmpty()) {
			return;
		}

		int id = settings[0].toInt() - 1;

		int n_group = this->data()->regulation_groups_.size();
		if (id < 0 || id >= n_group) {
			G_WARN("Illegal ID.");
			return;
		}

		Eigen::ArrayXd group_activity = this->data()->regulation_group_normalized_.mat_.row(id);

		auto matches = motif_position->get_motif_match_times(this->data()->regulation_groups_[id]);

		auto motif_names = motif_position->motif_names_;

		auto valid_motif_index = _Cs which(_Cs in(motif_names, rna_normalized->rownames_));
		if (valid_motif_index.isEmpty()) {
			G_WARN("No Motif in RNA Data.");
			return;
		}

		matches = _Cs reordered(matches, valid_motif_index);
		motif_names = _Cs reordered(motif_names, valid_motif_index);

		int n_motif = motif_names.size();
		QVector<double> possibility(n_motif);

		for (int i = 0; i < n_motif; ++i) {
			possibility[i] = _Cs correlation_pearson(group_activity, rna_normalized->get_row(motif_names[i]));

			possibility[i] *= (int)(matches[i] > 0);
		}

		auto nan_filter = _Cs sapply(possibility, [](auto d) {return !std::isnan(d); });

		motif_names = _Cs sliced(motif_names, nan_filter);
		possibility = _Cs sliced(possibility, nan_filter);

		auto order = _Cs order(possibility, true);

		CustomMatrix res;
		res.set_rownames(_Cs reordered(motif_names, order));
		res.update("Possibility", _Cs reordered(possibility, order));

		MatrixWindow::show_matrix(&res, "Trancription Factor Possibility of Regulation Group " + QString::number(id + 1), this->signal_emitter_);
	}
	else {
		auto settings = CommonDialog::get_response(
			this->signal_emitter_,
			"Show Transcription Factor Possibility",
			{ "Group ID" , "Show part?:no", "Select part" },
			{soap::InputStyle::StringLineEdit, soap::InputStyle::SwitchButton,
			soap::InputStyle::FactorChoice},
			{},
			{factor_info}
		);

		if (settings.isEmpty()) {
			return;
		}

		int id = settings[0].toInt() - 1;
		bool show_part = switch_to_bool(settings[1]);
		Eigen::ArrayX<bool> filter;
		if (show_part) {
			auto [factor, levels] = factor_choice_to_pair(settings[2]);

			if (levels.isEmpty()) {
				return;
			}

			filter = _Cs in(single_cell_multiome->metadata()->mat_.get_qstring(factor), levels);
		}

		int n_group = this->data()->regulation_groups_.size();
		if (id < 0 || id >= n_group) {
			G_WARN("Illegal ID.");
			return;
		}

		auto matches = motif_position->get_motif_match_times(this->data()->regulation_groups_[id]);

		auto motif_names = motif_position->motif_names_;
		auto valid_motif_index = _Cs which(_Cs in(motif_names, rna_normalized->rownames_));
		if (valid_motif_index.isEmpty()) {
			G_WARN("No Motif in RNA Data.");
			return;
		}

		matches = _Cs reordered(matches, valid_motif_index);
		motif_names = _Cs reordered(motif_names, valid_motif_index);
		int n_motif = motif_names.size();

		Eigen::ArrayXd group_activity = this->data()->regulation_group_normalized_.mat_.row(id);
		if (show_part) {
			group_activity = _Cs sliced(group_activity, filter);
		}

		QVector<double> possibility(n_motif);

		for (int i = 0; i < n_motif; ++i) {

			if (show_part) {

				

				possibility[i] = _Cs correlation_pearson(group_activity, _Cs sliced(rna_normalized->get_row(motif_names[i]), filter));
			
				possibility[i] *= (int)(matches[i] > 0);
			}
			else {

				possibility[i] = _Cs correlation_pearson(group_activity, rna_normalized->get_row(motif_names[i]));

				possibility[i] *= (int)(matches[i] > 0);
			}
		}

		auto nan_filter = _Cs sapply(possibility, [](auto d) {return !std::isnan(d); });

		motif_names = _Cs sliced(motif_names, nan_filter);
		possibility = _Cs sliced(possibility, nan_filter);

		auto order = _Cs order(possibility, true);

		CustomMatrix res;
		res.set_rownames(_Cs reordered(motif_names, order));
		res.update("Possibility", _Cs reordered(possibility, order));

		MatrixWindow::show_matrix(&res, "Trancription Factor Possibility of Regulation Group " + QString::number(id + 1), this->signal_emitter_);

	}
};

void CiceroItem::s_differential_analysis() {

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
		{ "Choose Group:1", "Minimum Percentage:0.1", "P Adjust Method"},
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

	G_LOG("Differential Analysis in " + factor_name + "...");
	DifferentialAnalysisWorker* worker = new DifferentialAnalysisWorker(
		this->data(),
		factor_name,
		metadata,
		comparison,
		p_adjust_method,
		minimum_percentage
	);

	G_LINK_WORKER_THREAD(DifferentialAnalysisWorker, x_differential_analysis_ready, CiceroItem, s_receive_differential_analysis);
};

void CiceroItem::s_receive_differential_analysis(DifferentialAnalysis da, QString name) {

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

void CiceroItem::s_umap() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("This data is not part of Single Cell Multiome Data.");
		return;
	}

	G_GETLOCK;

	SimpleUmapWorker* worker = new SimpleUmapWorker(this->data()->regulation_group_normalized_.mat_.transpose());
	G_LINK_WORKER_THREAD(SimpleUmapWorker, x_umap_ready, CiceroItem, s_receive_umap);
};

void CiceroItem::s_receive_umap(Eigen::MatrixXd emb) {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		return;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto rna_counts = single_cell_multiome->rna_counts();

	if (rna_counts == nullptr) {
		return;
	}

	QString title = this->signal_emitter_->get_unique_name(VARIABLE_CICERO_UMAP);

	DATA_SUBMODULES(Embedding)[title] = Embedding(
		Embedding::DataType::Umap,
		emb,
		rna_counts->colnames_,
		_Cs paste("Cicero UMAP-", _Cs cast<QString>(_Cs seq_n(1, emb.cols()))));

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

void CiceroItem::s_reset_threshold() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("Cannot reset threshold without SingleCellMultiome Data.");
		return;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto atac_counts = single_cell_multiome->atac_counts();

	if (atac_counts == nullptr) {
		G_WARN("No Atac Count Data Found.");
		return;
	}

	G_GETLOCK;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Reset Threshold for CCAN",
		{ "Threshold(0~1)" },
		{soap::InputStyle::NumericLineEdit}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	double threshold = settings[0].toDouble();

	if (threshold <= 0.0 || threshold >= 1.0) {
		G_WARN("Illegal Threshold Value");
		G_UNLOCK;
		return;
	}

	G_LOG("Calculating New CCANs...");

	QVector<QVector<int>> new_regulation_groups;

	Eigen::ArrayXi cluster = _Cs cluster_louvain_igraph(this->data()->connections_, threshold);

	int n_vertice = cluster.size();

	auto unique_cluster = _Cs unique(cluster);

	QVector<int> peak_index;

	for (auto c : unique_cluster) {
		peak_index.clear();

		for (int i = 0; i < n_vertice; ++i) {
			if (cluster[i] == c) {
				peak_index << i;
			}
		}

		if (peak_index.size() > 2) {
			new_regulation_groups << peak_index;
		}
	}

	if (new_regulation_groups.isEmpty()) {
		G_WARN("No Result for new threshold");
		G_UNLOCK;
		return;
	}

	this->__clear_reserve_data();

	DATA_SUBMODULES(Embedding).clear();
	DATA_SUBMODULES(DifferentialAnalysis).clear();

	this->data()->regulation_groups_ = new_regulation_groups;

	int nrow = new_regulation_groups.size();

	G_LOG(QString::number(nrow) + " CCAN Found.");

	int n_cell = atac_counts->cols();

	this->data()->regulation_group_counts_.mat_.resize(nrow, n_cell);
	this->data()->regulation_group_counts_.mat_.setZero();

#pragma omp parallel for
	for (int i = 0; i < nrow; ++i) {
		for (auto ind : this->data()->regulation_groups_[i]) {
			this->data()->regulation_group_counts_.mat_.row(i) += atac_counts->mat_.row(ind);
		}
	}

	int n_group = this->data()->regulation_group_counts_.rows();

	this->data()->regulation_group_counts_.colnames_ = atac_counts->colnames_;
	this->data()->regulation_group_counts_.rownames_ = _Cs cast<QString>(_Cs seq_n(1, n_group));

	this->data()->regulation_group_normalized_.colnames_ = this->data()->regulation_group_counts_.colnames_;
	this->data()->regulation_group_normalized_.rownames_ = this->data()->regulation_group_counts_.rownames_;
	this->data()->regulation_group_normalized_.mat_ = this->data()->regulation_group_counts_.mat_.cast<double>();

	double mean = this->data()->regulation_group_normalized_.mat_.sum() / n_cell;
	Eigen::ArrayXd depth = this->data()->regulation_group_normalized_.mat_.colwise().sum();
	for (int i = 0; i < n_cell; ++i) {
		if (depth[i] != 0.0) {
			this->data()->regulation_group_normalized_.mat_.col(i) *= mean / depth[i];
		}
	}

	this->data()->regulation_group_normalized_.mat_ = log(this->data()->regulation_group_normalized_.mat_.array() + 1.0).eval();

	this->data()->peak_names_ = atac_counts->rownames_;

	this->signal_emitter_->update_information(this->title_, soap::VariableType::Cicero, this->data());

	this->__check_data();

	G_UNLOCK;
};

void CiceroItem::s_ccan_embedding_feature_plot() {

	auto emb_names = this->index_tree_->search(soap::VariableType::Embedding);

	if (emb_names.isEmpty()) {
		G_WARN("No Embedding for visualisation.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Settings",
		{ "Ccan ID", "Embedding" },
		{soap::InputStyle::IntegerLineEdit, soap::InputStyle::ComboBox},
		{emb_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	int id = settings[0].toInt() - 1;

	int n_group = this->data()->regulation_groups_.size();
	if (id < 0 || id >= n_group) {
		G_WARN("Illegal ID.");
		return;
	}

	Eigen::ArrayXd group_activity = this->data()->regulation_group_normalized_.mat_.row(id);

	auto item = (EmbeddingItem*)this->signal_emitter_->get_item(settings[1]);

	item->show(group_activity, "CCAN " + settings[0], false, "Activity");
};

void CiceroItem::s_scan_downstream_ccan() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("No enough Data.");
		return;
	}

	auto single_cell_multiome = this->get_root<SingleCellMultiome>();

	auto motif_position = single_cell_multiome->motif_position();

	if (motif_position == nullptr) {
		G_WARN("Motif Position is not located.");
		return;
	}

	auto rna_normalized = single_cell_multiome->rna_normalized();

	if (rna_normalized == nullptr) {
		G_WARN("RNA Data should be normalized.");
		return;
	}

	auto valid_names = _Cs intersect(motif_position->motif_names_, rna_normalized->rownames_);

	if (valid_names.isEmpty()) {
		G_WARN("No Valid Trancription Factor Names Found.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Scanning Settings",
		{ "Transcription Factor Name" },
		{soap::InputStyle::LineEditWithCompleter},
		{valid_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString tf_name = settings[0];
	if (!valid_names.contains(tf_name)) {
		G_WARN("Invalid Transcription Factor Name");
		return;
	}

	Eigen::ArrayXd tf_exp = rna_normalized->get_row(tf_name);

	int tf_ind = motif_position->motif_names_.indexOf(tf_name);

	int n_ccan = this->data()->regulation_groups_.size();

	Eigen::ArrayX<bool> bind(n_ccan);

#pragma omp parallel for
	for (int i = 0; i < n_ccan; ++i) {
		bind(i) = motif_position->get_match_count(this->data()->regulation_groups_[i], tf_ind) > 0;
	}

	if (bind.count() == 0) {
		G_NOTICE("No ccan has binding motif.");
		return;
	}

	Eigen::ArrayXd correlation = Eigen::ArrayXd::Zero(n_ccan);

#pragma omp parallel for
	for (int i = 0; i < n_ccan; ++i) {
		if (!bind(i)) {
			continue;
		}
		correlation[i] = _Cs correlation_pearson(this->data()->regulation_group_normalized_.mat_.row(i), tf_exp);
	}

	for (auto&& c : correlation) {
		if (std::isnan(c)) {
			c = 0.0;
		}
	}

	auto order = _Cs order(correlation, true);

	order = _Cs sliced(order, bind);

	CustomMatrix res;
	res.set_rownames(_Cs paste("CCAN", _Cs cast<QString>(_Cs seq_n(1, n_ccan)), " "));
	res.update("Correlation", _Cs cast<QVector>(correlation));

	res.row_reorder(order);

	MatrixWindow::show_matrix(&res, "Potential downstream ccan of " + tf_name, this->signal_emitter_);

};