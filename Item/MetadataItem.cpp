#include "MetadataItem.h"

#include "ComparePlotDialog.h"

#include "FileIO.h"
#include "FileWritingWorker.h"
#include "MetadataViewWindow.h"
#include "ItemIOWorker.h"
#include "CommonDialog.h"
#include "MatrixWindow.h"
#include "CustomPlot.h"
#include "WilcoxTest.h"

void MetadataItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->mat_,
		this->title_,
		this->signal_emitter_);
}

void MetadataItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("Compare Plot", s_compare_plot);

	ADD_MAIN_ACTION("View", s_view);

	ADD_MAIN_ACTION("Extract", s_extract);

	ADD_MAIN_MENU("Edit");

	ADD_ACTION("Delete Metadata", "Edit", s_delete_metadata);
	ADD_ACTION("Change Data Type", "Edit", s_change_data_type);
	ADD_ACTION("Change Data Name", "Edit", s_change_data_name);
	ADD_ACTION("Slice Data", "Edit", s_slice_data);
	ADD_ACTION("Remove First", "Edit", s_remove_first);
	ADD_ACTION("Remove Last", "Edit", s_remove_last);
	ADD_ACTION("Remove All", "Edit", s_remove_all);
	ADD_ACTION("Remove Prefix", "Edit", s_remove_prefix);
	ADD_ACTION("Remove Suffix", "Edit", s_remove_suffix);
	ADD_ACTION("Remove From First", "Edit", s_remove_from_first);
	ADD_ACTION("Remove From Last", "Edit", s_remove_from_last);
	ADD_ACTION("Remove Until First", "Edit", s_remove_until_first);
	ADD_ACTION("Remove Until Last", "Edit", s_remove_until_last);
	ADD_ACTION("Row Name as Column", "Edit", s_rownames_as_column);
	ADD_ACTION("Add Prefix", "Edit", s_add_prefix);
	ADD_ACTION("Add Suffix", "Edit", s_add_suffix);


	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

}


void MetadataItem::s_rownames_as_column() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Set Name",
		{ "Name" },
		{ soap::InputStyle::StringLineEdit }
	);

	if (settings.isEmpty()) {
		return;
	}

	auto& mat = this->data()->mat_;
	mat.update(settings[0], mat.rownames_);

	this->__s_update_interface();
};

void MetadataItem::s_delete_metadata() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Select deleted metadata",
		{ "Deleted" },
		{ soap::InputStyle::ComboBox },
		{ this->data()->mat_.colnames_ }
	);

	if (settings.isEmpty()) {
		return;
	}

	this->data()->mat_.remove(settings[0]);

	this->__s_update_interface();


};

void MetadataItem::s_view() {

	if (this->attached_to(soap::VariableType::SingleCellRna)) {
		MetadataViewWindow::view(this->trace_back<SingleCellRna>(1), this->signal_emitter_);
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
		MetadataViewWindow::view(this->trace_back<SingleCellAtac>(1), this->signal_emitter_);
	}
	else if (this->attached_to(soap::VariableType::SingleCellMultiome)) {
		MetadataViewWindow::view(this->trace_back<SingleCellMultiome>(1), this->signal_emitter_);
	}
	else if (this->is_atomic()) {
		G_NOTICE("This metadata is not attached to any project.");
		return;
	}
};

void MetadataItem::s_extract() {

	G_GETLOCK;

	G_UNLOCK;

	DataFrame* dataFrame = new DataFrame(this->data()->mat_);

	this->signal_emitter_->x_data_create_soon(dataFrame, soap::VariableType::DataFrame, this->title_);

}

void MetadataItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rownames_.size()) + " | " + QString::number(this->data()->mat_.colnames_.size()) + " ]");

}

void MetadataItem::__s_delete_this() {

	if (this->attached_to(soap::VariableType::SingleCellRna)) {
		G_WARN("Metadata Item cannot be deleted from SingleCellRna.");
		return;
	}
	else if (this->attached_to(soap::VariableType::SingleCellMultiome)) {
		G_WARN("Metadata Item cannot be deleted from SingleCellMultiome.");
		return;
	}
	else if (this->attached_to(soap::VariableType::SingleCellAtac)) {
		G_WARN("Metadata Item cannot be deleted from SingleCellAtac.");
		return;
	}

	G_GETLOCK;

	G_UNLOCK;

	this->__remove_this();
};

void MetadataItem::s_compare_plot() {

	const auto& metadata = this->data()->mat_;

	auto factor_info = metadata.get_factor_information(false);

	if (factor_info.isEmpty()) {
		G_WARN("No Factor has more than one level.");
		return;
	}

	FeatureHandler handler;
	if (this->stem_from(soap::VariableType::SingleCellRna)) {
		handler.set(this->get_root<SingleCellRna>());
	}
	else if (this->stem_from(soap::VariableType::SingleCellAtac)) {
		handler.set(this->get_root<SingleCellAtac>());
	}
	else if (this->stem_from(soap::VariableType::SingleCellMultiome)) {
		handler.set(this->get_root<SingleCellMultiome>());
	}
	else if (this->stem_from(soap::VariableType::BulkRna)) {
		handler.set(this->get_root<BulkRna>());
	}
	else {
		handler.set(this->data());
	}

	auto features = handler.get_feature_names();

	if (features.numeric_names.isEmpty()) {
		G_WARN("No Valid Feature.");
		return;
	}

	auto settings = ComparePlotDialog::get_response(
		factor_info,
		features.numeric_names
	);

	auto&& [feature, normalize, use_gene_activity, show_p_value, group, type, comparisons] = settings;

	if (feature.isEmpty()) {
		return;
	}

	if (type == "Custom" && comparisons.isEmpty()) {
		return;
	}

	auto feature_data = handler.get_data({ feature, normalize, use_gene_activity });
	if (!feature_data.is_continuous()) {
		G_WARN("Invalid Feature");
		return;
	}

	auto factor_data = handler.get_data({ group });
	if (!factor_data.is_factor()) {
		G_WARN("Invalid Factor");
		return;
	}

	auto&& gs = this->draw_suite_->graph_settings_;

	auto values = feature_data.get_continuous();
	auto factors = factor_data.get_factor();
	auto levels = factor_data.get_levels();
	int n_level = levels.size();
	auto colors = gs.palette(levels);
	QList<Eigen::ArrayXd> vals = _Cs sapply(levels, [&values, &factors](auto&& level) {
		return _Cs cast<Eigen::ArrayX>(_Cs sliced(values, _Cs equal(factors, level)));
	});
	QList<std::tuple<int, int, double>> sigs;
	if (type == "All") {
		for (int i = 1; i < n_level; ++i) {
			for (int j = 0; j < n_level - i; ++j) {
				double p = wilcox_test(vals[j], vals[j + i]);
				sigs << std::make_tuple(j, j + i, p);
			}
		}
	}
	else if (type == "Significant Only") {
		for (int i = 1; i < n_level; ++i) {
			for (int j = 0; j < n_level - i; ++j) {
				double p = wilcox_test(vals[j], vals[j + i]);
				if (p >= 0.05) {
					continue;
				}
				sigs << std::make_tuple(j, j + i, p);
			}
		}
	}
	else {
		for (auto&& [c1, c2] : comparisons) {
			auto ind1 = levels.indexOf(c1);
			auto ind2 = levels.indexOf(c2);

			if (ind1 == -1 || ind2 == -1 || ind1 == ind2) {
				G_WARN("Illegal Comparison : " + c1 + " vs. " + c2);
				return;
			}

			if (ind1 > ind2) {
				std::swap(ind1, ind2);
			}

			double p = wilcox_test(vals[ind1], vals[ind2]);
			sigs << std::make_tuple(ind1, ind2, p);
		}
	}

	auto [draw_area, axis_rect, legend_layout] = _Cp prepare(gs);

	double min{ 0.0 }, max{ 0.0 };
	Eigen::ArrayXd bottom_ind;

	if (gs.use_boxplot()) {
		auto [m, n] = std::ranges::minmax(values);
		min = m;
		max = n;

		_CpPatch box_batch(
			draw_area,
			axis_rect,
			factors,
			levels,
			colors,
			_Cs cast<Eigen::ArrayX>(values),
			1.0,
			2.0,
			gs.boxplot_draw_outlier(),
			gs.get_scatter_point_size()
		);
		axis_rect->axis(QCPAxis::atBottom)->setRange(0, 2 * n_level);

		bottom_ind = Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1);

		_Cp set_bottom_axis_label(
			axis_rect,
			bottom_ind,
			levels,
			6,
			gs
		);
		
	}
	else {

		if (n_level == 2 && gs.use_facet_violin_plot()) {
			auto [m, n] = _CpPatch violin_facet(
				draw_area,
				axis_rect,
				factors,
				levels,
				colors,
				_Cs cast<Eigen::ArrayX>(values),
				1.0);

			min = m;
			max = n;

			bottom_ind.resize(2);
			bottom_ind[0] = 0.6;
			bottom_ind[1] = 1.4;

			axis_rect->axis(QCPAxis::atBottom)->setRange(0, 2);
			_Cp set_bottom_title(axis_rect, group, gs, true);
		}
		else {
			auto [m, n] = _CpPatch violin_batch(
				draw_area,
				axis_rect,
				factors,
				levels,
				colors,
				_Cs cast<Eigen::ArrayX>(values),
				1.0,
				2.0);

			min = m;
			max = n;

			axis_rect->axis(QCPAxis::atBottom)->setRange(0, 2 * n_level);

			bottom_ind = Eigen::ArrayXd::LinSpaced(n_level, 1, 2 * n_level - 1);
			_Cp set_bottom_axis_label(
				axis_rect,
				bottom_ind,
				levels,
				6,
				gs
			);
		}
	}

	_Cp set_simple_axis_no_title(axis_rect, gs);

	if (!sigs.isEmpty()) {

		QFont label_font = gs.get_scatter_label_font();

		double sep = (max - min) / 30;

		double anno_start = max + sep;
		auto [w, h] = this->draw_suite_->current_plot_->size();
		auto [anno_width, anno_height] = _CpUtility calculate_text_size("n.s.", label_font,
			axis_rect->axis(QCPAxis::atBottom)->range().size(), max - min, w, h);
		anno_height *= 2.0;
		sep = anno_height * 0.3;

		for (auto&& [i1, i2, p] : sigs) {
			_CpPatch line(
				draw_area,
				axis_rect,
				QVector<double>{bottom_ind[i1], bottom_ind[i1], bottom_ind[i2], bottom_ind[i2]},
				QVector<double>{anno_start, anno_start + sep, anno_start + sep, anno_start},
				Qt::black,
				2
			);

			QString l{ "n.s." };

			if (p < 0.001) {
				l = "***";
			}
			else if (p < 0.01) {
				l = "**";
			}
			else if (p < 0.05) {
				l = "*";
			}

			if (show_p_value) {
				l += " p = ";
				l += QString::number(p, 103, 4);
			}

			_CpPatch add_label(draw_area, axis_rect, l, 
				(bottom_ind[i1] + bottom_ind[i2]) / 2, anno_start + sep,
				label_font,
				Qt::AlignHCenter | Qt::AlignBottom);

			anno_start += anno_height * 2.0; 
		}

		axis_rect->axis(QCPAxis::atLeft)->setRange(_CpUtility get_range(min, anno_start));
	}
	else {
		axis_rect->axis(QCPAxis::atLeft)->setRange(_CpUtility get_range(min, max));
	}

	_Cp add_round_legend(draw_area, legend_layout, levels, colors, group, gs);
	_Cp set_left_title(axis_rect, feature, gs, true);
	_Cp add_title(draw_area, feature, gs);
	this->draw_suite_->update(draw_area);
};

void MetadataItem::s_remove_prefix() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "To delete" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString to_delete = settings[1];

	if (to_delete.isEmpty()) {
		return;
	}
	int delete_size = to_delete.size();

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		if (d.startsWith(to_delete)) {
			d = d.sliced(delete_size);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_suffix() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "To delete" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString to_delete = settings[1];

	if (to_delete.isEmpty()) {
		return;
	}
	int delete_size = to_delete.size();

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		if (d.endsWith(to_delete)) {
			d = d.sliced(0, d.size() - delete_size);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_until_first() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "Until First" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString ut = settings[1];

	if (ut.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		int ind = d.indexOf(ut);

		if (ind != -1) {
			d = d.sliced(ind);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_from_first() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "From First" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString ff = settings[1];

	if (ff.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		int ind = d.indexOf(ff);

		if (ind != -1) {
			d = d.sliced(0, ind);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_until_last() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "Until Last" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString ut = settings[1];

	if (ut.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		int ind = d.lastIndexOf(ut);

		if (ind != -1) {
			d = d.sliced(ind);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_from_last() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "From Last" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString fl = settings[1];

	if (fl.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		int ind = d.lastIndexOf(fl);

		if (ind != -1) {
			d = d.sliced(0, ind);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_all() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "To delete" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString to_delete = settings[1];

	if (to_delete.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		d.remove(to_delete);
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_first() {

	G_GETLOCK;

	G_UNLOCK;

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "To delete" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString to_delete = settings[1];

	if (to_delete.isEmpty()) {
		return;
	}

	int delete_size = to_delete.size();

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		auto ind = d.indexOf(to_delete);

		if (ind != -1) {
			d = d.sliced(0, ind) + d.sliced(ind + delete_size);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_add_suffix() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "Suffix" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString suffix = settings[1];

	if (suffix.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		d += suffix;
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_add_prefix() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "Prefix" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString prefix = settings[1];

	if (prefix.isEmpty()) {
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		d = prefix + d;
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_remove_last() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "To delete" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ valid_names }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];
	QString to_delete = settings[1];

	if (to_delete.isEmpty()) {
		return;
	}

	int delete_size = to_delete.size();

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		auto ind = d.lastIndexOf(to_delete);

		if (ind != -1) {
			d = d.sliced(0, ind) + d.sliced(ind + delete_size);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_slice_data() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid Metadata.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Edit Settings",
		{ "Feature", "Type", "Start:1","end:-1" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox, soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit },
		{ valid_names, { "Reserve", "Delete" } }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString feature_name = settings[0];

	bool is_delete = settings[1] == "Delete";

	int start = settings[2].toInt();
	if (start == 0) {
		G_WARN("Illegal Start Location.");
		return;
	}

	int end = settings[3].toInt();

	if (end == 0) {
		G_WARN("Illegal End Location.");
		return;
	}

	if (start > 0 && end > 0 && start > end) {
		G_WARN("Illegal Slice");
		return;
	}

	if (start < 0 && end < 0 && start > end) {
		G_WARN("Illegal Slice");
		return;
	}

	auto data = metadata.get_qstring(feature_name);

	for (auto&& d : data) {

		if (d.isEmpty()) {
			G_WARN("Slice Can not be performed.");
			return;
		}

		int size = d.size();

		int s{ 0 }, e{ 0 };// [s,e)

		if (start > 0) {
			s = start - 1;
		}
		else {
			s = size + start;
		}

		if (end > 0) {
			e = end;
		}
		else {
			e = size + end + 1;
		}

		if (s >= size || s < 0) {
			G_WARN("Slice Can not be performed in " + d);
			return;
		}

		if (e > size || e < 1) {
			G_WARN("Slice Can not be performed in " + d);
			return;
		}

		if (s >= e) {
			G_WARN("Slice Can not be performed in " + d);
			return;
		}

		if (is_delete) {

			if (e == size) {
				if (s == 0) {
					d = "";
				}
				else {
					d = d.sliced(0, s);
				}
			}
			else {
				if (s == 0) {
					d = d.sliced(e);
				}
				else {
					d = d.sliced(0, s) + d.sliced(e);
				}
			}
		}
		else {
			d = d.sliced(s, e - s);
		}
	}

	metadata.update(feature_name, data, metadata.data_type_[feature_name]);
};

void MetadataItem::s_change_data_name() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Change Data Name",
		{ "Original Name", "New Name" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit },
		{ this->data()->mat_.colnames_ }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString original_name = settings[0];
	QString new_name = settings[1];
	if (new_name.isEmpty()) {
		return;
	}
	if (this->data()->mat_.colnames_.contains(new_name)) {
		G_WARN(new_name + " is already existed in data.");
		return;
	}

	this->data()->mat_.change_name(original_name, new_name);
}

void MetadataItem::s_change_data_type() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Change Data Type",
		{ "Data", "New Type" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox },
		{ this->data()->mat_.colnames_, { "Integer", "Numeric", "Integer Factor", "String", "String Factor" } }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString name = settings[0];
	QString new_type = settings[1];
	if (new_type == "Integer") {
		this->data()->mat_.to_integer_numeric(name);
	}
	else if (new_type == "Numeric") {
		this->data()->mat_.to_double_numeric(name);
	}
	else if (new_type == "Integer Factor") {
		this->data()->mat_.to_integer_factor(name);
	}
	else if (new_type == "String") {
		this->data()->mat_.to_string(name);
	}
	else if (new_type == "String Factor") {
		this->data()->mat_.to_string_factor(name);
	}

	G_LOG("Data type of " + name + " changed.");
};

