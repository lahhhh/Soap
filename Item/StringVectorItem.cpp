#include "StringVectorItem.h"

#include "CustomMatrix.h"
#include "qCustomPlot.h"

#include "ItemIOWorker.h"
#include "FileWritingWorker.h"
#include "CommonDialog.h"
#include "MatrixWindow.h"
#include "TextEditWindow.h"

void StringVectorItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_MENU("Edit");
	ADD_ACTION("Manually", "Edit", s_edit_manually);
	ADD_ACTION("Sort", "Edit", s_edit_sort);
	ADD_ACTION("Remove Empty Elements", "Edit", s_edit_remove_empty_elements);

	ADD_MAIN_MENU("Convert");
	ADD_ACTION("Gene Name", "Convert", s_convert_to_genenames);

	ADD_MAIN_MENU("Cooperate");
	ADD_ACTION("Intersect", "Cooperate", s_intersect_with);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void StringVectorItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->data_.size()) + " ]");
};

void StringVectorItem::__show_this() {
	MatrixWindow::show_matrix(
		&this->data()->data_,
		this->title_,
		this->signal_emitter_
	);
}

void StringVectorItem::s_intersect_with() {

	auto& variables = this->signal_emitter_->variable_information_;

	QStringList available_data = this->signal_emitter_->get_type_variable(soap::VariableType::StringVector).keys();
	if (available_data.size() < 2) {
		G_LOG("No enough string vector found.");
		return;
	}

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"String Vector Intersect",
		{ "Data"},
		{ soap::InputStyle::SimpleChoice },
		{ available_data }
	);
	if (settings.isEmpty())return;

	available_data = simple_choice_to_list(settings[0]);
	if (available_data.size() < 2)return;
	if (!_Cs is_unique(available_data)) {
		G_WARN("Can not intersect the same data!");
		return;
	}

	auto its = this->signal_emitter_->search(available_data);

	auto choosed = _Cs sapply(its, [](IndexTree* it) {return static_cast<const StringVector*>(it->data_); });

	QStringList res = choosed[0]->data_;

	for (int i = 1; i < choosed.size(); ++i) {
		res = _Cs intersect(res, choosed[i]->data_);
	}

	this->signal_emitter_->x_data_create_soon(new StringVector(res), soap::VariableType::StringVector, "Intersected");
};

void StringVectorItem::s_edit_manually() {
	G_GETLOCK;

	TextEditWindow::view(
		&this->data()->data_,
		this->data_,
		this->signal_emitter_
	);

	G_UNLOCK;
};

void StringVectorItem::s_edit_sort() {
	G_GETLOCK;

	this->data()->data_ = _Cs sorted(this->data()->data_);	

	G_UNLOCK;
};

void StringVectorItem::s_convert_to_genenames() {

	this->signal_emitter_->x_data_create_soon(new GeneName(this->data()->data_), soap::VariableType::GeneName, "Converted From StringVector");
};

void StringVectorItem::s_edit_remove_empty_elements() {

	G_GETLOCK;

	this->data()->data_.removeAll("");

	this->__s_update_interface();

	G_UNLOCK;
};