#include "VariableItem.h"

#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "FileWritingWorker.h"

VariableItem::VariableItem(	
	void* data,
	soap::VariableType data_type,
	const QString& title, 
	IndexTree* index_tree,
	const QStringList& item,
	PlotsSuite* draw_suite, 
	InformationTextBrowser* information_area, 
	SignalEmitter* signal_emitter
) :
	QTreeWidgetItem(item),
	data_(data),
	data_type_(data_type),
	title_(title),
	index_tree_(index_tree),
	draw_suite_(draw_suite),
	information_area_(information_area),
	signal_emitter_(signal_emitter)	
{};

void VariableItem::__duplicate_this() {};

bool VariableItem::attached_to(soap::VariableType type) const{

	if (this->is_atomic()) {

		return false;
	}

	return this->index_tree_->parent_->type_ == type;
};

bool VariableItem::stem_from(soap::VariableType type) const {

	return this->index_tree_->root()->type_ == type;
};

bool VariableItem::is_atomic() const {

	return this->index_tree_->parent_ == nullptr;
};

void VariableItem::__remove_information() {

	int count = this->childCount();

	for (int i = 0; i < count; ++i) {

		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		
		item->__remove_information();
	}

	this->signal_emitter_->remove_information(this->title_);
};

void VariableItem::__change_data_name(const QString& new_name) {};

void VariableItem::__check_data() {};

void VariableItem::__clear_reserve_data() {

	this->index_tree_->children_.clear();

	this->__data_delete_soon();

	this->__remove_information();

	while (this->childCount() > 0) {
		
		this->removeChild(this->child(0));
	}
};

void VariableItem::__set_menu() {}

void VariableItem::set_item(VariableItem* item) {

	this->addChild(item);

	this->signal_emitter_->update_information(item->title_, item->data_type_, item->data_);
}

void VariableItem::__data_delete_soon() {

	int count = this->childCount();

	for (int i = 0; i < count; ++i) {

		VariableItem* item = static_cast<VariableItem*>(this->child(i));
		
		item->__data_delete_soon();
	}

	this->signal_emitter_->x_data_delete_soon(this->data_, this->data_type_, this);
}

void VariableItem::__remove_data() {};

void VariableItem::__remove_this() {

	this->__data_delete_soon();
	
	this->__remove_information();

	this->signal_emitter_->remove_index(this->index_tree_);

	this->__remove_data();

	delete this;
};

bool VariableItem::__check_lock() {

	if (!this->signal_emitter_->try_lock(this->index_tree_)) {

		G_WARN("Please waiting for the computation in progress.");

		return false;
	}

	return true;
}

void VariableItem::__show_this(){}

void VariableItem::__s_update_interface() {}

void VariableItem::__s_rename() {

	G_GETLOCK;

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Set New Name for this item",
		{ "New Name:" + this->title_ },
		{ soap::InputStyle::StringLineEdit}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QString new_name = settings[0];
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
};

void VariableItem::__s_duplicate() {

	G_GETLOCK;

	G_UNLOCK;
	
	this->__duplicate_this();
};

void VariableItem::__s_delete_this() {

	G_GETLOCK;

	G_UNLOCK;

	this->__remove_this();
};

void VariableItem::__s_export_as_item() {

	G_GETLOCK;

	QString export_item_path = QFileDialog::getSaveFileName(nullptr, "Set File Name for this item", "", "soap item(*.sif)");
	if (export_item_path.isEmpty()) {

		G_UNLOCK;
		return;
	}
	
	ItemIOWorker* worker = new ItemIOWorker(export_item_path, this->data_type_, this->data_, this->title_);
	G_LINK_WORKER_THREAD_NO_RESPONSE(ItemIOWorker, VariableItem);
};

void VariableItem::__s_export_as_csv() {

	G_GETLOCK;

	QString export_csv_path = QFileDialog::getSaveFileName(nullptr, "Please Set File Name", "", "CSV(*.csv)");

	if (export_csv_path.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"CSV Settings",
		{ "quote:yes", "keep rownames:yes", "keep colnames:yes" },
		{soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QStringList write_settings;

	bool quote = switch_to_bool(settings[0]);
	if (!quote) {
		settings << CSV_NO_QUOTE;
	}

	bool keep_rownames = switch_to_bool(settings[1]);
	if (!keep_rownames) {
		settings << CSV_NO_ROWNAMES;
	}

	bool keep_colnames = switch_to_bool(settings[2]); 
	if (!keep_colnames) {
		settings << CSV_NO_COLNAMES;
	}

	FileWritingWorker* worker = new FileWritingWorker(
		this->data_, 
		this->data_type_, 
		export_csv_path, 
		soap::FileType::CSV, 
		settings
	);
	G_LINK_WORKER_THREAD_NO_RESPONSE(FileWritingWorker, VariableItem);
};

void VariableItem::__s_export_as_tsv() {
	G_GETLOCK;

	QString export_tsv_path = QFileDialog::getSaveFileName(nullptr, "Please Set File Name", "", "TSV(*.tsv)");

	if (export_tsv_path.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"TSV Settings",
		{ "quote:yes", "keep rownames:yes", "keep colnames:yes" },
		{soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	QStringList write_settings;

	bool quote = switch_to_bool(settings[0]);
	if (!quote) {
		settings << TSV_NO_QUOTE;
	}

	bool keep_rownames = switch_to_bool(settings[1]);
	if (!keep_rownames) {
		settings << TSV_NO_ROWNAMES;
	}

	bool keep_colnames = switch_to_bool(settings[2]);
	if (!keep_colnames) {
		settings << TSV_NO_COLNAMES;
	}

	FileWritingWorker* worker = new FileWritingWorker(
		this->data_,
		this->data_type_,
		export_tsv_path,
		soap::FileType::TSV,
		settings
	);
	G_LINK_WORKER_THREAD_NO_RESPONSE(FileWritingWorker, VariableItem);
};

void VariableItem::__s_receive_message(QString message, int mode) {

	if (mode == 0) {

		G_LOG(message);
	}
	else if (mode == 1) {

		G_WARN(message);
	}
	else {

		G_NOTICE(message);
	}

};

void VariableItem::__s_work_finished() {

	this->signal_emitter_->x_update_interface();

	G_UNLOCK;
};