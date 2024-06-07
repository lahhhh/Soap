#include "MetadataItem.h"

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

	ADD_MAIN_ACTION("Compare Violin Plot", s_compare_violin_plot);

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
	ADD_ACTION("Remove If Start With", "Edit", s_remove_if_start_with);
	ADD_ACTION("Remove If End With", "Edit", s_remove_if_end_with);
	ADD_ACTION("Remove From First", "Edit", s_remove_from_first);
	ADD_ACTION("Remove From Last", "Edit", s_remove_from_last);
	ADD_ACTION("Remove Until First", "Edit", s_remove_until_first);
	ADD_ACTION("Remove Until Last", "Edit", s_remove_until_last);
	ADD_ACTION("Add Prefix", "Edit", s_add_prefix);
	ADD_ACTION("Add Suffix", "Edit", s_add_suffix);


	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

}

void MetadataItem::s_delete_metadata() {

	G_GETLOCK;

	G_UNLOCK;

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

void MetadataItem::s_compare_violin_plot() {

	const auto& metadata = this->data()->mat_;

	auto factor_info = metadata.get_factor_information(false);

	if (factor_info.isEmpty()) {
		G_WARN("No Factor has more than one level.");
		return;
	}

	auto valid_features = metadata.get_type_names(CustomMatrix::DataType::IntegerNumeric) <<
		metadata.get_type_names(CustomMatrix::DataType::DoubleNumeric);

	if (valid_features.isEmpty()) {
		G_WARN("No Valid Feature.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Facet Violin Plot Settings",
		{ "Features", "Factor", "Show Significance:yes", "facet:yes" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::FactorChoice,
		 soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton },
		{ valid_features },
		{ factor_info }
	);

	if (settings.isEmpty()) {
		return;
	}

	auto feature = settings[0];

	auto [factor, levels] = factor_choice_to_pair(settings[1]);
	if (levels.size() != 2) {
		return;
	}

	bool show_significance = switch_to_bool(settings[2]);
	bool facet = switch_to_bool(settings[3]);

	QStringList group = metadata.get_qstring(factor);

	auto& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = _Cp prepare(gs);

	auto feat = _Cs cast<Eigen::ArrayX>(metadata.get_double(feature));
	_Cp set_simple_axis_no_title(axis_rect, gs);

	int control_point_number = 16;

	auto colors = gs.palette(levels);

	if (facet) {

		auto [min, max] = _CpPatch violin_facet(
			draw_area,
			axis_rect,
			group,
			levels,
			colors,
			feat,
			1.0,
			control_point_number
		);

		if (show_significance) {

			Eigen::ArrayXd e1 = _Cs sliced(feat, _Cs equal(group, levels[0]));
			Eigen::ArrayXd e2 = _Cs sliced(feat, _Cs equal(group, levels[1]));

			double p = wilcox_test(e1, e2);

			double val_span = max - min;
			double marker_loc = max + 0.1 * val_span;
			_CpPatch line(draw_area, axis_rect,
				QVector<double>{ 0.4, 1.6},
				QVector<double>{marker_loc, marker_loc},
				Qt::black,
				2);

			QString sig{ "n.s." };
			if (p < 0.001) {
				sig = "***";
			}
			else if (p < 0.01) {
				sig = "**";
			}
			else if (p < 0.05) {
				sig = "*";
			}

			_CpPatch add_label(draw_area, axis_rect, sig, 1.0, marker_loc,
				gs.get_scatter_label_font(),
				Qt::AlignHCenter | Qt::AlignBottom);

			max += 0.2 * val_span;
		}

		_CpPatch set_range(axis_rect, QCPRange(0.0, 2.0), _CpUtility get_range(min, max));
		_Cp set_bottom_axis_label(
			axis_rect,
			_Cs cast<Eigen::ArrayX>(QVector<double>{1.0}),
			{ feature },
			6,
			gs
		);
	}
	else {

		auto [min, max] = _CpPatch violin_batch(
			draw_area,
			axis_rect,
			group,
			levels,
			colors,
			feat,
			1.0,
			2.0,
			16
		);

		if (show_significance) {

			Eigen::ArrayXd e1 = _Cs sliced(feat, _Cs equal(group, levels[0]));
			Eigen::ArrayXd e2 = _Cs sliced(feat, _Cs equal(group, levels[1]));

			double p = wilcox_test(e1, e2);

			double val_span = max - min;
			double marker_loc = max + 0.1 * val_span;
			_CpPatch line(draw_area, axis_rect,
				QVector<double>{ 0.9, 3.1},
				QVector<double>{marker_loc, marker_loc},
				Qt::black,
				2);

			QString sig{ "n.s." };
			if (p < 0.001) {
				sig = "***";
			}
			else if (p < 0.01) {
				sig = "**";
			}
			else if (p < 0.05) {
				sig = "*";
			}

			_CpPatch add_label(draw_area, axis_rect, sig, 2.0, marker_loc,
				gs.get_scatter_label_font(),
				Qt::AlignHCenter | Qt::AlignBottom);

			max += 0.2 * val_span;
		}

		_CpPatch set_range(axis_rect, QCPRange(0.0, 4.0), _CpUtility get_range(min, max));

		_Cp set_bottom_axis_label(
			axis_rect,
			_Cs cast<Eigen::ArrayX>(QVector<double>{1.0, 3.0}),
			levels,
			6,
			gs
		);
	}

	_Cp add_round_legend(draw_area, legend_layout, levels, colors, factor, gs);
	_Cp set_left_title(axis_rect, feature, gs, true);
	_Cp add_title(draw_area, feature, gs);
	this->draw_suite_->update(draw_area);
};

void MetadataItem::s_remove_if_start_with() {

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

void MetadataItem::s_remove_if_end_with() {

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

