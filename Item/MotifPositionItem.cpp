#include "MotifPositionItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "FileWritingWorker.h"
#include "ItemIOWorker.h"
#include "CustomPlot.h"
#include "TranscriptionalFactorFootprintingWorker.h"
#include "ChromVARWorker.h"
#include "GenomeUtility.h"

void MotifPositionItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->n_peak()) + " | " + QString::number(this->data()->n_motif()) + " ]");
}

void MotifPositionItem::__check_data() {

	this->check_variable(DATA_SUBMODULES(Footprint));
	this->check_variable(DATA_SUBMODULES(ChromVAR));

}

void MotifPositionItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Show Motif", s_show_motif);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);	

	ADD_MAIN_ACTION("Show Motif Peak Location", s_show_motif_location);
	ADD_MAIN_ACTION("Show Motif Binding Sites", s_show_motif_binding_location);

	ADD_MAIN_ACTION("Show Motif Weight", s_show_motif_matrix);

	ADD_MAIN_ACTION("Show Motif Overlapping Rate", s_show_motif_overlapping_rate);

	ADD_MAIN_ACTION("Show Binding Site Sequence", s_show_typical_binding_sequence);

	if (this->attached_to(soap::VariableType::SingleCellMultiome) || this->attached_to(soap::VariableType::SingleCellAtac)) {

		ADD_MAIN_ACTION("Run ChromVAR", s_chromvar);

		ADD_MAIN_ACTION("Show Difference Footprinting", s_difference_footprinting);

		ADD_MAIN_ACTION("Show Transcriptional Factor Footprinting", s_show_transcriptional_factor_footprinting);

		ADD_MAIN_ACTION("Batch Footprinting Estimation", s_batch_task);

		ADD_MAIN_ACTION("Show Multiple Footprint", s_multiple_footprint_plot);
	}

	ADD_MAIN_ACTION("Delete", __s_delete_this);
};

void MotifPositionItem::s_show_motif_matrix() {

	if (this->data()->motifs_.empty()) {
		G_WARN("Database is Embty.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Motif Name",
		{ "Motif Name" },
		{ soap::InputStyle::LineEditWithCompleter},
		{ this->data()->motif_names_}
	);
	if (settings.isEmpty()) {
		return;
	}
	QString motif_name = settings[0];
	if (!this->data()->motifs_.contains(motif_name)) {
		G_WARN(motif_name + " is not found in database.");
		return;
	}

	auto&& pwm = this->data()->motifs_[motif_name];
	int motif_size = pwm.weight_.mat_.cols();

	MatrixWindow::show_matrix(&pwm.weight_.mat_, {"A", "C", "G", "T"},
		_Cs cast<QString>(_Cs linspaced(motif_size, 1, motif_size)), motif_name, this->signal_emitter_);
}

void MotifPositionItem::s_show_motif_binding_location() {

	auto motif_names = this->data()->motif_names_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Motif Settings",
		{ "Motif Name" },
		{soap::InputStyle::LineEditWithCompleter},
		{motif_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString motif_name = settings[0];

	if (!motif_names.contains(motif_name)) {
		G_WARN("Illegal motif name.");
		return;
	}

	auto loc = this->data()->get_position(motif_name);

	if (loc.is_empty()) {
		G_WARN("No Binding Sites.");
		return;
	}

	auto sites = loc.get_range_names();
	std::ranges::sort(sites);

	MatrixWindow::show_matrix(&sites, "Binding Sites of " + motif_name, this->signal_emitter_);
};

void MotifPositionItem::s_show_motif_location() {

	auto motif_names = this->data()->motif_names_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Motif Settings",
		{ "Motif Name" },
		{soap::InputStyle::LineEditWithCompleter},
		{motif_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString motif_name = settings[0];

	if (!motif_names.contains(motif_name)) {
		G_WARN("Illegal motif name.");
		return;
	}

	auto res = this->data()->get_match_peak_position(settings[0]);

	if (res.isEmpty()) {
		G_WARN("No Motif Location");
	}
	else {

		std::ranges::sort(res);

		MatrixWindow::show_matrix(&res, "Location of " + settings[0], this->signal_emitter_);
	}
};

void MotifPositionItem::__show_this() {

	MatrixWindow::show_matrix(
		this->data(), 
		this->data()->peak_names_, 
		this->data()->motif_names_,
		this->title_, 
		this->signal_emitter_, 
		false, 
		this->data_
	);
};

void MotifPositionItem::s_show_typical_binding_sequence() {

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {
		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

		if (single_cell_multiome->species_ != soap::Species::Human) {
			G_WARN("Now only support human genome.");
			return;
		}

	}
	else {
		G_NOTICE("Can not find species information. use human genome.")
	}

	QStringList tf_names = this->data()->motif_names_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose Transcription Factor",
		{ "Transcription Factor"},
		{ soap::InputStyle::LineEditWithCompleter},
		{tf_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString tf_name = settings[0];

	if (!tf_names.contains(tf_name)) {
		G_WARN("Illegal Transcription Factor Name.");
	}

	GenomicRange loc = this->data()->get_position(tf_name);

	int n_loc = loc.size();

	if (n_loc > 20) {
		loc.row_slice(0, 20);
	}

	n_loc = loc.size();

	auto seqs = _Cs get_human_grch38_sequence(loc);

	int count{ 0 };

	for (int i = 0; i < n_loc; ++i) {

		if (count > 9) {
			break;
		}

		if (seqs[i].isEmpty()) {
			continue;
		}

		++count;

		auto [name, start, end, strand] = loc.at(i);

		G_LOG("Sequence [" + QString::number(count) + "] " + name + ":" + QString::number(start) + "-" + QString::number(end)
			+ " " + strand + " " + seqs[i]);
	}

};

void MotifPositionItem::s_difference_footprinting() {

	const Metadata* metadata{ nullptr };
	Bias bias;
	soap::Species species;
	const Fragments* fragments{ nullptr };

	QMap<QString, QStringList> map;

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);
		metadata = single_cell_multiome->metadata();
		bias = single_cell_multiome->insertion_bias_;
		species = single_cell_multiome->species_;

		fragments = single_cell_multiome->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments is not loaded.");
			return;
		}

		map = metadata->mat_.get_factor_information();

		if (map.isEmpty()) {
			G_LOG("No suitable metadata detected.");
			return;
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);
		metadata = single_cell_atac->metadata();
		bias = single_cell_atac->insertion_bias_;
		species = single_cell_atac->species_;

		fragments = single_cell_atac->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments is not loaded.");
			return;
		}

		map = metadata->mat_.get_factor_information();

		if (map.isEmpty()) {
			G_LOG("No suitable metadata detected.");
			return;
		}
	}
	else {
		G_WARN("Illegal data status.");
		return;
	}

	G_GETLOCK;

	auto pic_name = QFileDialog::getSaveFileName(nullptr, "Set Picture Name", "", "PNG(*.png)");

	if (pic_name.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QStringList tf_names = this->data()->motif_names_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Footprinting Settings",
		{ "Transcription Factor 1(draw)", "Transcription Factor 2(exclude)", 
		"Factor:Factor"},
		{ soap::InputStyle::LineEditWithCompleter, soap::InputStyle::LineEditWithCompleter, 
		soap::InputStyle::FactorChoice},
		{tf_names, tf_names}, { map }
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString tf1 = settings[0];
	QString tf2 = settings[1];

	if (!tf_names.contains(tf1) || !tf_names.contains(tf2)) {
		G_UNLOCK;
		G_WARN("Illegal Transcription Factor Name.");
		return;
	}

	GenomicRange g1 = this->data()->get_position(tf1);
	GenomicRange g2 = this->data()->get_position(tf2);

	g1.subtract(g2);

	if (g1.is_empty()) {
		G_WARN("No Location after removing overlapping regions.");
		G_UNLOCK;
		return;
	}

	G_LOG(QString::number(g1.size()) + " binding sites found after removing binding sites of " + tf2 + " in " + tf1 + " binding sites");

	auto [factor_name, levels] = factor_choice_to_pair(settings[2]);
	if (levels.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QStringList factors = metadata->mat_.get_qstring(factor_name);

	TranscriptionalFactorFootprintingWorker* worker = new TranscriptionalFactorFootprintingWorker(
		metadata,
		bias,
		species,
		g1,
		pic_name,
		fragments,
		factor_name,
		levels,
		factors,
		this->draw_suite_->graph_settings_,
		700,
		900
	);

	G_LINK_WORKER_THREAD_NO_RESPONSE(TranscriptionalFactorFootprintingWorker, MotifPositionItem);
};

void MotifPositionItem::s_batch_task() {


	const Metadata* metadata{ nullptr };
	Bias bias;
	soap::Species species;
	const Fragments* fragments{ nullptr };

	QMap<QString, QStringList> map;

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);
		metadata = single_cell_multiome->metadata();
		bias = single_cell_multiome->insertion_bias_;
		species = single_cell_multiome->species_;

		fragments = single_cell_multiome->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments is not loaded.");
			return;
		}

		map = metadata->mat_.get_factor_information();

		if (map.isEmpty()) {
			G_LOG("No suitable metadata detected.");
			return;
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);
		metadata = single_cell_atac->metadata();
		bias = single_cell_atac->insertion_bias_;
		species = single_cell_atac->species_;

		fragments = single_cell_atac->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments is not loaded.");
			return;
		}

		map = metadata->mat_.get_factor_information();

		if (map.isEmpty()) {
			G_LOG("No suitable metadata detected.");
			return;
		}
	}
	else {
		G_WARN("Illegal data status.");
		return;
	}

	G_GETLOCK;

	QString directory = QFileDialog::getExistingDirectory(nullptr, "Choose Output Directory");

	if (directory.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Batch Task Settings",
		{ "Factor:Factor" },
		{ soap::InputStyle::FactorChoice},
		QList<QStringList>(), { map }
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

	QStringList factors = metadata->mat_.get_qstring(factor_name);

	TranscriptionalFactorFootprintingWorker* worker = new TranscriptionalFactorFootprintingWorker(
		metadata,
		bias,
		species,
		this->data(),
		this->data()->motif_names_,
		fragments,
		directory,
		factor_name,
		levels,
		factors,
		this->draw_suite_->graph_settings_,
		700,
		900
	);

	G_LINK_WORKER_THREAD_NO_RESPONSE(TranscriptionalFactorFootprintingWorker, MotifPositionItem);
};

void MotifPositionItem::s_show_transcriptional_factor_footprinting() {
	

	const Metadata* metadata{ nullptr };
	Bias bias;
	soap::Species species;
	const Fragments* fragments{ nullptr };

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {

		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);
		metadata = single_cell_multiome->metadata();
		bias = single_cell_multiome->insertion_bias_;
		species = single_cell_multiome->species_;

		fragments = single_cell_multiome->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments is not loaded.");
			return;
		}
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {

		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);
		metadata = single_cell_atac->metadata();
		bias = single_cell_atac->insertion_bias_;
		species = single_cell_atac->species_;

		fragments = single_cell_atac->fragments();
		if (fragments == nullptr) {
			G_WARN("Fragments is not loaded.");
			return;
		}
	}
	else {
		G_WARN("Illegal data status.");
		return;
	}

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Transcriptional Factor Name",
		{ "Transcriptional Factor Name" },
		{ soap::InputStyle::LineEditWithCompleter},
		{ this->data()->motif_names_}
	);
	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	QString transcriptional_factor_name = settings[0];
	if (!this->data()->contains(transcriptional_factor_name)) {
		G_WARN("Transcriptional factor : " + transcriptional_factor_name + " is not found in database.");
		G_UNLOCK;
		return;
	}

	TranscriptionalFactorFootprintingWorker* worker = new TranscriptionalFactorFootprintingWorker(
		metadata,
		bias,
		species,
		this->data(),
		{ transcriptional_factor_name },
		fragments
	);
	G_LINK_WORKER_THREAD(TranscriptionalFactorFootprintingWorker, x_footprint_ready, MotifPositionItem, s_receive_transcriptional_factor_footprinting);
}

void MotifPositionItem::s_receive_chromvar(ChromVAR* chrom_var) {

	QString title = this->signal_emitter_->get_unique_name("ChromVAR");

	DATA_SUBMODULES(ChromVAR)[title] = std::move(*chrom_var);

	ChromVARItem* item = new ChromVARItem(
		title,
		this->index_tree_,
		&DATA_SUBMODULES(ChromVAR)[title],
		this->draw_suite_,
		this->information_area_,
		this->signal_emitter_
	);

	this->set_item(item);
};

void MotifPositionItem::s_receive_transcriptional_factor_footprinting(Footprint footprint) 
{
	QString transcriptional_factor_name = footprint.motif_.motif_name_;
	QString title = this->signal_emitter_->get_unique_name(transcriptional_factor_name);
	
	DATA_SUBMODULES(Footprint)[title] = footprint;

	FootprintItem* item = new FootprintItem(
		title, 
		this->index_tree_, 
		&DATA_SUBMODULES(Footprint)[title], 
		this->draw_suite_, 
		this->information_area_, 
		this->signal_emitter_
	);

	this->set_item(item);
};

void MotifPositionItem::s_show_motif() {

	if (this->data()->motifs_.empty()) {
		G_WARN("Database is Embty.");
		return;
	}
	if (this->data()->motifs_.cbegin()->second.data_type_ != PatternWeightMatrix::DataType::Frequency) {
		G_NOTICE("Only Frequency is supported in showing motif.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Set Motif Name",
		{ "Motif Name" },
		{ soap::InputStyle::LineEditWithCompleter},
		{ this->data()->motif_names_}
	);
	if (settings.isEmpty()) {
		return;
	}
	QString motif_name = settings[0];
	if (!this->data()->motifs_.contains(motif_name)) {
		G_WARN(motif_name + " is not found in database.");
		return;
	}

	const PatternWeightMatrix& pwm = this->data()->motifs_[motif_name];

	if (pwm.data_type_ != PatternWeightMatrix::DataType::Frequency) {
		return;
	}

	int length = pwm.weight_.cols();

	auto&& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect] = _Cp prepare_ar(gs);

	_CpPatch set_range(axis_rect, QCPRange(0, length * 100), QCPRange(0, 400));
	_CpPatch remove_left_bottom_axis(axis_rect);

	QStringList characters = pwm.weight_.rownames_;

	QVector<QPixmap> acgt(4);
	QImage img(FILE_CHARACTER_A_PNG);
	acgt[0] = QPixmap::fromImage(img);

	img.load(FILE_CHARACTER_C_PNG);
	acgt[1] = QPixmap::fromImage(img);

	img.load(FILE_CHARACTER_G_PNG);
	acgt[2] = QPixmap::fromImage(img);

	img.load(FILE_CHARACTER_T_PNG);
	acgt[3] = QPixmap::fromImage(img);

	for (int i = 0; i < length; ++i) {
		Eigen::ArrayXd col = pwm.weight_.mat_.col(i); // frequency
		auto col_order = _Cs order(col);
		double now = 0;
		double sum = col.sum();
		double uncertainty{ 0.0 };
		for (auto d : col) {
			if (d != 0.0) {
				uncertainty -= d / sum * log2(d / sum);
			}
		}

		for (auto c : col_order) {
			double frequency = col[c] / sum;
			if (frequency > 0) {
				double height = frequency * (2 - uncertainty);

				if (height > 0.01) {
					QCPItemPixmap* pic_t = new QCPItemPixmap(draw_area);
					pic_t->setPixmap(acgt[c].scaled(100, height * 400));
					pic_t->setVisible(true);
					pic_t->setScaled(true, Qt::IgnoreAspectRatio);
					pic_t->bottomRight->setCoords(i * 100 + 100, now);
					pic_t->topLeft->setCoords(i * 100, now + height * 400);
					now += height * 400;
				}
			}
		}

	}
	_CpPatch set_range(axis_rect, QCPRange(0, length * 100), QCPRange(0, 800));
	_Cp add_title(draw_area, motif_name, gs);
	this->draw_suite_->update(draw_area);
};

void MotifPositionItem::footprint_plot_patch(
	QCustomPlot* draw_area,
	QCPLayoutGrid* layout,
	const QString& transcriptional_factor_name,
	const QString& factor_name,
	const QStringList& levels,
	const QStringList& factors,
	bool show_in_group,
	const Eigen::ArrayX<bool>& show_filter,
	const QList<QColor>& colors
) {
	QCPAxisRect* axis_rect = new QCPAxisRect(draw_area);
	_CpPatch set_border_only(axis_rect, Qt::black, 2);

	layout->addElement(0, 0, axis_rect);
	draw_area->addGraph(axis_rect->axis(QCPAxis::atBottom), axis_rect->axis(QCPAxis::atLeft));

	const Footprint* data{ nullptr };
	for (auto&& fp : DATA_SUBMODULES(Footprint)) {
		if (fp.second.motif_.motif_name_ == transcriptional_factor_name) {
			data = &fp.second;
		}
	}

	if (data == nullptr) {
		return;
	}

	const int nrow = data->insertion_matrix_.mat_.rows(), ncol = data->insertion_matrix_.mat_.cols();
	if (nrow != factors.size() || ncol <= 495) {
		return ;
	}
	Eigen::ArrayXd flank_mean = data->insertion_matrix_.mat_.block(0, 0, nrow, 50).rowwise().sum().cast<double>() + data->insertion_matrix_.mat_.block(0, ncol - 50, nrow, 50).rowwise().sum().cast<double>();
	flank_mean /= 100;
	double all_mean = flank_mean.mean();
	if (all_mean == 0) {
		return;
	}
	for (int i = 0; i < nrow; ++i) {
		if (flank_mean[i] == 0) {
			flank_mean[i] = all_mean;
		}
	}
	Eigen::MatrixXd normalized = data->insertion_matrix_.mat_.cast<double>();
	normalized.array().colwise() /= flank_mean;

	QList<Eigen::ArrayXd> locations;
	QStringList levels_use;
	QList<QColor> color_use;
	const int n_level = levels.size();
	for (int i = 0; i < n_level; ++i) {

		auto&& factor = levels[i];

		auto index = show_in_group ? 
			_Cs which(_Cs equal(factors, factor) * show_filter) : _Cs match(factors, factor);

		if (index.isEmpty()) {
			continue;
		}

		Eigen::MatrixXd sub_matrix = normalized(index, Eigen::all);

		Eigen::ArrayXd means = sub_matrix.colwise().mean();

		levels_use << factor;

		color_use << colors[i];

		locations << means - _Cs cast<Eigen::ArrayX>(data->expected_insertions_);
	}

	if (locations.isEmpty()) {
		return;
	}

	QVector<double> x_axis = _Cs cast<double>(data->insertion_matrix_.colnames_);

	const int ncolor = color_use.size();
	for (int i = 0; i < ncolor; ++i) {
		_CpPatch line(draw_area, axis_rect, x_axis, _Cs cast<QVector>(locations[i]), color_use[i], 2);
	}
	double min_ob = std::ranges::min(_Cs sapply(locations, [](auto&& arr) {return arr.minCoeff(); }));
	double max_ob = std::ranges::max(_Cs sapply(locations, [](auto&& arr) {return arr.maxCoeff(); }));

	auto [min_x, max_x] = std::ranges::minmax(x_axis);

	_CpPatch set_range(axis_rect, QCPRange(min_x, max_x), QCPRange(1.1 * min_ob - 0.1 * max_ob, 1.1 * max_ob - 0.1 * min_ob));

	_CpPatch add_title(draw_area, layout, transcriptional_factor_name, this->draw_suite_->graph_settings_.get_title_font());
};

void MotifPositionItem::s_chromvar() {

	G_GETLOCK;

	if (this->attached_to(soap::VariableType::SingleCellMultiome)) {
		SingleCellMultiome* single_cell_multiome = this->trace_back<SingleCellMultiome>(1);

		auto atac_counts = single_cell_multiome->atac_counts();
		if (atac_counts == nullptr) {
			G_WARN("No ATAC Counts Data.");
			return;
		}

		ChromVARWorker* worker = new ChromVARWorker(single_cell_multiome->species_, this->data(), atac_counts);
		G_LINK_WORKER_THREAD(ChromVARWorker, x_chromvar_ready, MotifPositionItem, s_receive_chromvar);
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
		SingleCellAtac* single_cell_atac = this->trace_back<SingleCellAtac>(1);

		auto atac_counts = single_cell_atac->counts();
		if (atac_counts == nullptr) {
			G_WARN("No ATAC Counts Data.");
			return;
		}

		ChromVARWorker* worker = new ChromVARWorker(single_cell_atac->species_, this->data(), atac_counts);
		G_LINK_WORKER_THREAD(ChromVARWorker, x_chromvar_ready, MotifPositionItem, s_receive_chromvar);
	}
	else {
		G_WARN("No Counts Data Found.");
		G_UNLOCK;
		return;
	}
};

void MotifPositionItem::s_multiple_footprint_plot() {

	if (!this->stem_from(soap::VariableType::SingleCellMultiome)) {
		G_WARN("Footprint can only be visualized in a single-cell multiome data object.");
		return;
	}
	SingleCellMultiome* single_cell_multiome = this->get_root<SingleCellMultiome>();
	const auto& metadata = single_cell_multiome->metadata()->mat_;

	QMap<QString, QStringList> map = metadata.get_factor_information();

	if (map.isEmpty()) {
		G_LOG("No suitable metadata detected.");
		return;
	}

	if (DATA_SUBMODULES(Footprint).size() < 2) {
		G_WARN("No Enough Footprint.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Footprint Plot Setting",
		{ "Footprint", "Factor:Factor", "Number of row:1", "Show in group:no", "Show Group" },
		{ soap::InputStyle::MultipleLineEditWithCompleter, soap::InputStyle::FactorChoice,
		soap::InputStyle::IntegerLineEdit, soap::InputStyle::SwitchButton, soap::InputStyle::FactorChoice},
		{ _Cs keys(DATA_SUBMODULES(Footprint))},
		{ map, map}
	);

	if (settings.isEmpty())return;

	QStringList footprints = multiple_line_edit_with_completer_to_list(settings[0]);

	if (footprints.isEmpty()) {
		return;
	}

	auto [factor_name, levels] = factor_choice_to_pair(settings[1]);
	if (levels.isEmpty())return;

	QStringList factors = metadata.get_qstring(factor_name);

	int nrow = settings[2].toInt();
	if (nrow < 0) {
		G_WARN("Number of rows can not be less than 1.");
		return;
	}

	bool show_in_group = switch_to_bool(settings[3]);

	auto [show_group_name, show_levels] = factor_choice_to_pair(settings[4]);
	if (show_in_group) {
		if (show_levels.isEmpty()) {
			return;
		}
	}

	QStringList show_group = metadata.get_qstring(show_group_name);
	auto show_filter = _Cs in(show_group, show_levels);

	QStringList valid_footprints;

	auto computed = _Cs sapply(DATA_SUBMODULES(Footprint),
		[](auto&& fp) {return fp.second.motif_.motif_name_; }
	);

	const int n_footprint = footprints.size();
	for (int i = 0; i < n_footprint; ++i) {
		if (computed.contains(footprints[i])) {
			valid_footprints << footprints[i];
		}
		else {
			G_WARN(footprints[i] + " is not found in Footprint data.");
		}
	}

	auto& gs = this->draw_suite_->graph_settings_;

	auto [draw_area, left_layout, legend_layout] = _Cp prepare_lg_lg(gs);

	const int ncolor = levels.size();
	auto colors = gs.palette(levels);

	const int n_valid_footprint = valid_footprints.size();
	for (int i = 0; i < n_valid_footprint; ++i) {
		QCPLayoutGrid* sub_layout = new QCPLayoutGrid;
		int col = i / nrow;
		int row = i - col * nrow;
		left_layout->addElement(row, col, sub_layout);

		this->footprint_plot_patch(
			draw_area,
			sub_layout,
			valid_footprints[i],
			factor_name,
			levels,
			factors,
			show_in_group,
			show_filter,
			colors
		);
	}

	_Cp add_round_legend(
		draw_area, 
		legend_layout, 
		levels, 
		colors,
		factor_name, 
		gs
	);

	this->draw_suite_->update(draw_area);
};

void MotifPositionItem::s_show_motif_overlapping_rate() {

	QStringList tf_names = this->data()->motif_names_;

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Overlap Calculation Settings",
		{ "Transcription Factor 1", "Transcription Factor 2" },
		{soap::InputStyle::LineEditWithCompleter, soap::InputStyle::LineEditWithCompleter},
		{tf_names, tf_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString tf1 = settings[0];
	QString tf2 = settings[1];

	if (!tf_names.contains(tf1) || !tf_names.contains(tf2)) {
		G_WARN("Illegal Transcription Factor Name.");
		return;
	}

	GenomicRange g1 = this->data()->get_position(tf1);
	GenomicRange g2 = this->data()->get_position(tf2);

	int overlap = g1.count_overlap(g2);
	int overlap2 = g2.count_overlap(g1);

	G_NOTICE(QString::number(overlap) + " / " + QString::number(g1.size()) + " match position in " + tf1 + " and " 
		+ QString::number(overlap2) + " / " + QString::number(g2.size()) + " match position in " + tf2 + " are overlapped.");


};
