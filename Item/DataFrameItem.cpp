#include "DataFrameItem.h"

#include "NumericMatrix.h"

#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "FileWritingWorker.h"
#include "MetadataViewWindow.h"
#include "MatrixWindow.h"

void DataFrameItem::__s_update_interface() {

	this->setText(2, "[ " + QString::number(this->data()->mat_.rows()) + " | " + QString::number(this->data()->mat_.cols()) + " ]");
}

void DataFrameItem::__show_this() {
	if (this->data()->mat_.rows() > 100000) {
		MatrixWindow::show_matrix(
			&this->data()->mat_,
			this->title_,
			this->signal_emitter_,
			false, this->data_);
	}
	else {
		MatrixWindow::show_matrix(
			&this->data()->mat_,
			this->title_,
			this->signal_emitter_);
	}
}

void DataFrameItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Rename", __s_rename);

	ADD_MAIN_ACTION("View", s_view);

	// edit

	ADD_MAIN_MENU("Edit");

	ADD_MENU("Edit | Adjust Shape", "Adjust Shape", "Edit");

	ADD_ACTION("Filter", "Edit | Adjust Shape", s_filter_by_feature);
	ADD_ACTION("Choose Columns Remained", "Edit | Adjust Shape", s_choose_columns_remained);
	ADD_ACTION("First Row as Column Names", "Edit | Adjust Shape", s_first_row_as_colname);
	ADD_ACTION("First Column as Row Names", "Edit | Adjust Shape", s_first_column_as_rowname);
	ADD_ACTION("Delete First Column", "Edit | Adjust Shape", s_delete_first_column);
	ADD_ACTION("Delete First Row", "Edit | Adjust Shape", s_delete_first_row);
	ADD_ACTION("Column as Row Names", "Edit | Adjust Shape", s_column_as_rowname);
	ADD_ACTION("Row Names as Column", "Edit | Adjust Shape", s_rownames_as_column);
	ADD_ACTION("Delete Column", "Edit | Adjust Shape", s_delete_column);

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
	ADD_ACTION("Add Prefix", "Edit", s_add_prefix);
	ADD_ACTION("Add Suffix", "Edit", s_add_suffix);

	ADD_MAIN_MENU("Promote to");

	ADD_MENU("Promote to | Metadata", "Metadata", "Promote to");
	ADD_ACTION("by Row Name", "Promote to | Metadata", s_promote_to_metadata_by_rowname);
	ADD_ACTION("by Metadata", "Promote to | Metadata", s_promote_to_metadata_by_metadata);

	ADD_MAIN_MENU("Convert to");

	ADD_ACTION("Genome Annotation", "Convert to", s_convert_to_genome_annotation);
	ADD_ACTION("Numeric Matrix", "Convert to", s_convert_to_numeric_matrix);
	ADD_ACTION("DenseDouble", "Convert to", s_convert_to_dense_double);
	ADD_ACTION("RNA Seq", "Convert to", s_convert_to_bulkrna);

	ADD_MAIN_MENU("Export");
	ADD_ACTION("as Item", "Export", __s_export_as_item);
	ADD_ACTION("as CSV", "Export", __s_export_as_csv);
	ADD_ACTION("as TSV", "Export", __s_export_as_tsv);

	ADD_MAIN_ACTION("Delete", __s_delete_this);

};

void DataFrameItem::s_view() {

	MetadataViewWindow::view(this->data(), this->signal_emitter_);
};

void DataFrameItem::s_change_data_name() {
	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Change Data Name",
		{ "Original Name", "New Name" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit},
		{ this->data()->mat_.colnames_}
	);
	if (settings.isEmpty()) {
		G_UNLOCK;
		return;
	}
	QString original_name = settings[0];
	QString new_name = settings[1];
	if (new_name.isEmpty()) {
		G_UNLOCK;
		return;
	}
	if (this->data()->mat_.colnames_.contains(new_name)) {
		G_WARN(new_name + " is already existed in data.");
		G_UNLOCK;
		return;
	}
	this->data()->mat_.change_name(original_name, new_name);
	G_UNLOCK;
}

void DataFrameItem::s_change_data_type() {
	G_GETLOCK;
	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Change Data Type",
		{ "Data", "New Type" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox},
		{ this->data()->mat_.colnames_, { "Integer", "Numeric", "Integer Factor", "String", "String Factor" }}
	);

	if (settings.isEmpty()) {
		G_UNLOCK;
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
	G_UNLOCK;
	G_LOG("Data type of " + name + " changed.");
};

void DataFrameItem::s_convert_to_dense_double() {

	G_GETLOCK;
	G_UNLOCK;

	if (this->data()->mat_.is_empty()) {
		G_WARN("This is an empty data frame.");
		return;
	}

	auto&& data = this->data()->mat_;

	int nrow = data.rows(), ncol = data.cols();

	DenseDouble* dd = new DenseDouble();

	dd->rownames_ = data.rownames_;
	dd->colnames_ = data.colnames_;

	auto&& mat = dd->mat_;
	mat.resize(nrow, ncol);

	for (int i = 0; i < ncol; ++i) {

		QString colname = data.colnames_[i];

		auto data_type = data.data_type_[colname];

		switch (data_type)
		{
		case CustomMatrix::DataType::DoubleNumeric:
		{
			auto&& d = data.get_const_double_reference(colname);
			for (int j = 0; j < nrow; ++j) {
				mat(j, i) = d[j];
			}
			break;
		}
		case CustomMatrix::DataType::IntegerNumeric:
		case CustomMatrix::DataType::IntegerFactor:
		{
			auto&& d = data.get_const_integer_reference(colname);
			for (int j = 0; j < nrow; ++j) {
				mat(j, i) = d[j];
			}
			break;
		}
		case CustomMatrix::DataType::QString:
		case CustomMatrix::DataType::QStringFactor:
		{
			auto&& d = data.get_const_qstring_reference(colname);
			for (int j = 0; j < nrow; ++j) {
				mat(j, i) = d[j].toDouble();
			}
			break;
		}
		default:
			break;
		}

	}

	this->signal_emitter_->x_data_create_soon(dd, soap::VariableType::DenseDouble, this->title_);
}

void DataFrameItem::s_convert_to_numeric_matrix() {

	G_GETLOCK;
	G_UNLOCK;

	if (this->data()->mat_.is_empty()) {
		G_WARN("This is an empty data frame.");
		return;
	}

	auto&& data = this->data()->mat_;

	int nrow = data.rows(), ncol = data.cols();

	NumericMatrix* nm = new NumericMatrix();
	
	auto&& mat = nm->data_;
	mat.resize(nrow, ncol);

	for (int i = 0; i < ncol; ++i) {

		QString colname = data.colnames_[i];

		auto data_type = data.data_type_[colname];

		switch (data_type)
		{
		case CustomMatrix::DataType::DoubleNumeric:
		{
			auto&& d = data.get_const_double_reference(colname);
			for (int j = 0; j < nrow; ++j) {
				mat(j, i) = d[j];
			}
			break;
		}
		case CustomMatrix::DataType::IntegerNumeric:
		case CustomMatrix::DataType::IntegerFactor:
		{
			auto&& d = data.get_const_integer_reference(colname);
			for (int j = 0; j < nrow; ++j) {
				mat(j, i) = d[j];
			}
			break;
		}
		case CustomMatrix::DataType::QString:
		case CustomMatrix::DataType::QStringFactor:
		{
			auto&& d = data.get_const_qstring_reference(colname);
			for (int j = 0; j < nrow; ++j) {
				mat(j, i) = d[j].toDouble();
			}
			break;
		}
		default:
			break;
		}

	}

	this->signal_emitter_->x_data_create_soon(nm, soap::VariableType::NumericMatrix, this->title_);
};

void DataFrameItem::s_convert_to_bulkrna() {

	G_GETLOCK;
	G_UNLOCK;

	auto types = custom::unique(custom::values(this->data()->mat_.data_type_));

	if (types.size() != 1 || types[0] != CustomMatrix::DataType::IntegerNumeric) {
		G_WARN("Data should contain only integer numeric data.");
		return;
	}

	auto&& data = this->data()->mat_;

	int nrow = data.rows(), ncol = data.cols();

	Eigen::MatrixXi mat(nrow, ncol);

	for (int i = 0; i < ncol; ++i) {

		mat.col(i) = custom::cast<Eigen::ArrayX>(data.get_const_integer_reference(data.colnames_[i]));
	}

	BulkRna* rna = new BulkRna();
	SUBMODULES(*rna, Metadata)[VARIABLE_METADATA] = Metadata(CustomMatrix(data.colnames_));
	SUBMODULES(*rna, DenseInt)[VARIABLE_COUNTS] = DenseInt{ DenseInt::DataType::Counts, mat, data.rownames_, data.colnames_ };
	this->signal_emitter_->x_data_create_soon(rna, soap::VariableType::BulkRna, "From " + this->title_);

};

void DataFrameItem::s_convert_to_genome_annotation() {

	G_GETLOCK;
	G_UNLOCK;

	if (!this->data()->mat_.contains({ "Seq Names", "Range Start", "Range End", "Strand" })) {
		G_WARN("If this dataframe needs to be converted to a genome annotation, its column names should contain \"Seq Names\", \"Range Start\", \"Range End\", and \"Strand\". You can \
		rename the proper data by menu [Edit >> Change Data Name]");
		return;
	}

	if (this->data()->mat_.data_type_["Seq Names"] != CustomMatrix::DataType::QStringFactor)
	{
		G_WARN("The type of [Seq Names] should be string factor.");
		return;
	}
	if (this->data()->mat_.data_type_["Range Start"] != CustomMatrix::DataType::IntegerNumeric)
	{
		G_WARN("The type of [Range Start] should be numeric integer.");
		return;
	}
	if (this->data()->mat_.data_type_["Range End"] != CustomMatrix::DataType::IntegerNumeric)
	{
		G_WARN("The type of [Range End] should be numeric integer.");
		return;
	}
	if (this->data()->mat_.data_type_["Strand"] != CustomMatrix::DataType::QStringFactor)
	{
		G_WARN("The type of [Strand] should be string factor.");
		return;
	}
	QStringList strand_level = this->data()->mat_.string_factors_["Strand"];
	if (strand_level.size() > 3) {
		G_WARN("Strand should only contain three level : +, -, *");
		return;
	}
	QStringList strand_normal_level;
	strand_normal_level << "+" << "-" << "*";
	for (const auto& level : strand_level) {
		if (!strand_normal_level.contains(level)) {
			G_WARN("Unexpected level in Strand : " + level + " .");
			return;
		}
	}

	GenomicRange* genomic_range = new GenomicRange(
		RunLengthEncoding<QString>::from_container(
			custom::sapply(
				this->data()->mat_.get_const_qstring_reference("Seq Names"), [](auto&& name) {return custom::standardize_chromosome_name(name); }
			)
		),
		IRange(this->data()->mat_.get_const_integer_reference("Range Start"),
			custom::minus(
				this->data()->mat_.get_const_integer_reference("Range End"), this->data()->mat_.get_const_integer_reference("Range Start")
			)
		),
		RunLengthEncoding<char>::from_container(
			custom::sapply(
				this->data()->mat_.get_const_qstring_reference("Strand"), [](const QString& str)->char {	return str[0].toLatin1();}
			)
		)
	);


	genomic_range->metadata_ = this->data()->mat_;
	genomic_range->metadata_.remove("Seq Names");
	genomic_range->metadata_.remove("Range Start");
	genomic_range->metadata_.remove("Range End");
	genomic_range->metadata_.remove("Strand");
	genomic_range->finalize();
	this->signal_emitter_->x_data_create_soon(genomic_range, soap::VariableType::GenomicRange, "From " + this->title_);
};

QString DataFrameItem::metadata_insert(
	const CustomMatrix& from,
	CustomMatrix& to,
	int style,
	int type,
	const QStringList& from_index,
	const QStringList& to_index
) {
	int match{ 0 };
	if (style == 0) {
		match = custom::intersect_length(from_index, to_index);

		if (match == 0) {
			return "No Content Matched!";
		}
	}
	else {
		match = std::min(from.rows(), to.rows());
	}

	if (match == to.rows()) {

		if (style == 0) {

			auto order = custom::index_of(to_index, from_index);
			for (auto&& [name, data_type] : from.data_type_) {
				if (data_type == CustomMatrix::DataType::IntegerNumeric) {
					to.update(name, custom::reordered(from.get_const_integer_reference(name), order), CustomMatrix::DataType::IntegerNumeric);
				}
				else if (data_type == CustomMatrix::DataType::DoubleNumeric) {
					to.update(name, custom::reordered(from.get_const_double_reference(name), order), CustomMatrix::DataType::DoubleNumeric);
				}
				else if (data_type == CustomMatrix::DataType::QStringFactor) {
					to.update(name, custom::reordered(from.get_const_qstring_reference(name), order), CustomMatrix::DataType::QStringFactor);
				}
				else if (data_type == CustomMatrix::DataType::QString) {
					to.update(name, custom::reordered(from.get_const_qstring_reference(name), order), CustomMatrix::DataType::QString);
				}
				else if (data_type == CustomMatrix::DataType::IntegerFactor) {
					to.update(name, custom::reordered(from.get_const_integer_reference(name), order), CustomMatrix::DataType::IntegerFactor);
				}
			}
		}
		else {
			for (auto&& [name, data_type] : from.data_type_) {
				if (data_type == CustomMatrix::DataType::IntegerNumeric) {
					to.update(name, from.get_const_integer_reference(name).sliced(0, match), CustomMatrix::DataType::IntegerNumeric);
				}
				else if (data_type == CustomMatrix::DataType::DoubleNumeric) {
					to.update(name, from.get_const_double_reference(name).sliced(0, match), CustomMatrix::DataType::DoubleNumeric);
				}
				else if (data_type == CustomMatrix::DataType::QStringFactor) {
					to.update(name, from.get_const_qstring_reference(name).sliced(0, match), CustomMatrix::DataType::QStringFactor);
				}
				else if (data_type == CustomMatrix::DataType::QString) {
					to.update(name, from.get_const_qstring_reference(name).sliced(0, match), CustomMatrix::DataType::QString);
				}
				else if (data_type == CustomMatrix::DataType::IntegerFactor) {
					to.update(name, from.get_const_integer_reference(name).sliced(0, match), CustomMatrix::DataType::IntegerFactor);
				}
			}
		}

		return "DataFrame Insertion Finished";
	}
	else if (type == 0) {
		return "Missing value deteced!";
	}

	QVector<int> order, from_order;
	int nrow = to.rows();
	if (style == 0) {
		order = custom::index_of(from_index, to_index);
		from_order = custom::match(order, [](auto t) {return t >= 0; });
		order.removeAll(-1);
	}
	else {
		order = from_order = custom::seq_n(0, match);
	}

	for (auto&& [name, data_type] : from.data_type_) {
		if (data_type == CustomMatrix::DataType::IntegerNumeric) {
			if (type == 2 && to.contains(name) && to.data_type_[name] == CustomMatrix::DataType::IntegerNumeric) {
				to.update(name, custom::assigned(to.get_const_integer_reference(name), order, from.get_const_integer_reference(name), from_order), CustomMatrix::DataType::IntegerNumeric);
			}
			else {
				to.update(name, custom::assigned(QVector<int>(nrow, -1), order, from.get_const_integer_reference(name), from_order), CustomMatrix::DataType::IntegerNumeric);
			}
		}
		else if (data_type == CustomMatrix::DataType::DoubleNumeric) {
			if (type == 2 && to.contains(name) && to.data_type_[name] == CustomMatrix::DataType::DoubleNumeric) {
				to.update(name, custom::assigned(to.get_const_double_reference(name), order, from.get_const_double_reference(name), from_order), CustomMatrix::DataType::DoubleNumeric);
			}
			else {
				to.update(name, custom::assigned(QVector<double>(nrow, -0.0), order, from.get_const_double_reference(name), from_order), CustomMatrix::DataType::DoubleNumeric);
			}
		}
		else if (data_type == CustomMatrix::DataType::QStringFactor) {
			if (type == 2 && to.contains(name) && to.data_type_[name] == CustomMatrix::DataType::QStringFactor) {
				to.update(name, custom::assigned(to.get_const_qstring_reference(name), order, from.get_const_qstring_reference(name), from_order), CustomMatrix::DataType::QStringFactor);
			}
			else {
				to.update(name, custom::assigned(QStringList(nrow, "NA"), order, from.get_const_qstring_reference(name), from_order), CustomMatrix::DataType::QStringFactor);
			}
		}
		else if (data_type == CustomMatrix::DataType::QString) {
			if (type == 2 && to.contains(name) && to.data_type_[name] == CustomMatrix::DataType::QString) {
				to.update(name, custom::assigned(to.get_const_qstring_reference(name), order, from.get_const_qstring_reference(name), from_order), CustomMatrix::DataType::QString);
			}
			else {
				to.update(name, custom::assigned(QStringList(nrow, "NA"), order, from.get_const_qstring_reference(name), from_order), CustomMatrix::DataType::QString);
			}
		}
		else if (data_type == CustomMatrix::DataType::IntegerFactor) {
			if (type == 2 && to.contains(name) && to.data_type_[name] == CustomMatrix::DataType::IntegerFactor) {
				to.update(name, custom::assigned(to.get_const_integer_reference(name), order, from.get_const_integer_reference(name), from_order), CustomMatrix::DataType::IntegerFactor);
			}
			else {
				to.update(name, custom::assigned(QVector<int>(nrow, -1), order, from.get_const_integer_reference(name), from_order), CustomMatrix::DataType::IntegerFactor);
			}
		}
	}

	return "Finished.";
};

void DataFrameItem::s_promote_to_metadata_by_metadata() {

	G_GETLOCK;

	QMap<QString, QStringList> choices;
	QMap<QString, Metadata*> metadatas;
	QStringList available_data;
	for (auto&& [name, info] : this->signal_emitter_->variable_information_) {
		auto type = info.first;
		if (type == soap::VariableType::Metadata) {
			Metadata* d = static_cast<Metadata*>(info.second);
			auto string_names = d->mat_.get_type_names(CustomMatrix::DataType::QString);
			if (!string_names.isEmpty()) {
				choices[name] = string_names;
				metadatas[name] = d;
			}
		}
	}
	if (choices.isEmpty()) {
		G_NOTICE("No suitable data found.");
		G_UNLOCK;
		return;
	}
	QStringList promote_settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose the data attached",
		{ "Data", "(if)Missing Value"},
		{ soap::InputStyle::FactorChoice, soap::InputStyle::ComboBox },
		{  { "Stop", "Fill", "Fill if No Original Value" } },
		{choices}
	);

	if (promote_settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto [data_name, level] = factor_choice_to_pair(promote_settings[0]);

	if (level.size() != 1) {
		G_UNLOCK;
		return;
	}

	auto* metadata = metadatas[data_name];
	auto& d = metadata->mat_;

	QString factor_name = level[0];

	auto factor = d.get_qstring(level[0]);

	if (!custom::is_unique(factor)) {
		G_WARN("Choosed Metadata is not unique for mapping");
		G_UNLOCK;
		return;
	}

	QString missing = promote_settings[1];
	int type{ 0 };

	if (missing == "Stop") {
		type = 0;
	}
	else if (missing == "Fill") {
		type = 1;
	}
	else {
		type = 2;
	}

	if (!this->signal_emitter_->try_lock(this->signal_emitter_->search(metadata))) {
		G_WARN("Please waiting for computation in progress.");
		G_UNLOCK;
		return;
	}

	G_LOG(metadata_insert(this->data()->mat_, d, 0, type, this->data()->mat_.rownames_, factor));

	this->signal_emitter_->unlock(this->signal_emitter_->search(metadata));
	this->signal_emitter_->x_update_interface();
	G_UNLOCK;
}

void DataFrameItem::s_promote_to_metadata_by_rowname() {
	G_GETLOCK;

	QStringList available_data;
	for (auto&& [name, info] : this->signal_emitter_->variable_information_) {
		auto type = info.first;
		if (type == soap::VariableType::Metadata) {
			available_data << name;
		}
	}
	if (available_data.isEmpty()) {
		G_NOTICE("No suitable data found.");
		G_UNLOCK;
		return;
	}
	QStringList promote_settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Choose the data attached",
		{ "Data", "Style", "(if)Missing Value" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::ComboBox, soap::InputStyle::ComboBox},
		{ available_data, { "by Row Name", "by Row Order" }, { "Stop", "Fill", "Fill if No Original Value" }}
	);

	if (promote_settings.isEmpty()) {
		G_UNLOCK;
		return;
	}

	auto metadata_name = promote_settings[0];
	QString promote_style = promote_settings[1];
	QString missing = promote_settings[2];

	int style{ 0 }, type{ 0 };
	if (promote_style == "by Row Name") {
		style = 0;
	}
	else {
		style = 1;
	}

	if (missing == "Stop") {
		type = 0;
	}
	else if(missing == "Fill") {
		type = 1;
	}
	else {
		type = 2;
	}

	Metadata* metadata = static_cast<Metadata*>(this->signal_emitter_->get_variable(metadata_name));

	if (!this->signal_emitter_->try_lock(this->signal_emitter_->search(metadata))) {
		G_WARN("Please waiting for computation in progress.");
		G_UNLOCK;
		return;
	}

	G_LOG(metadata_insert(this->data()->mat_, metadata->mat_, style, type, this->data()->mat_.rownames_, metadata->mat_.rownames_));

	this->signal_emitter_->unlock(this->signal_emitter_->search(metadata));
	this->signal_emitter_->x_update_interface();
	G_UNLOCK;
};

void DataFrameItem::s_delete_first_row() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(this->data()->mat_.rows(), true);
	filter[0] = false;
	this->data()->mat_.row_slice(filter);

	this->__s_update_interface();
}

void DataFrameItem::s_delete_first_column() {
	
	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto& mat = this->data()->mat_;
	mat.remove(mat.colnames_[0]);

	this->__s_update_interface();
};

void DataFrameItem::s_first_column_as_rowname() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto& mat = this->data()->mat_;
	mat.rownames_ = custom::make_unique(mat.get_qstring(mat.colnames_[0]));

	this->__s_update_interface();
};

void DataFrameItem::s_choose_columns_remained() {
	
	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	QStringList settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose Colunms Remained",
		{"Columns"},
		{soap::InputStyle::MultiCheckBox},
		{custom::paste(this->data()->mat_.colnames_, ":true")}
	);

	if (settings.isEmpty()) {
		return;
	}

	QStringList columns_remained = multi_check_box_to_list(settings[0]);

	if (columns_remained.isEmpty()) {
		G_WARN("No Column is selected.");
		return;
	}

	this->data()->mat_.reserve(columns_remained);

	this->__s_update_interface();

};

void DataFrameItem::s_delete_column() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_, 
		"Choose Column",
		{ "Column" },
		{ soap::InputStyle::ComboBox},
		{ this->data()->mat_.colnames_}
	);

	if (setting.isEmpty()) {
		return;
	}

	this->data()->mat_.remove(setting[0]);

	this->__s_update_interface();
};

void DataFrameItem::s_rownames_as_column() {

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

void DataFrameItem::s_column_as_rowname() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto settings = CommonDialog::get_response(
		this->signal_emitter_, 
		"Choose Column",
		{ "Column" },
		{ soap::InputStyle::ComboBox},
		{ this->data()->mat_.colnames_}
	);

	if (settings.isEmpty()) {
		return;
	}

	auto& mat = this->data()->mat_;
	mat.rownames_ = custom::make_unique(mat.get_qstring(settings[0]));

	this->__s_update_interface();
};

void DataFrameItem::col_slice(const Eigen::ArrayX<bool>& filter, bool in_place) {

	if (in_place) {
		this->data()->mat_.col_slice(filter);

		this->__s_update_interface();
	}
	else {

		auto* df = new DataFrame();

		df->mat_ = this->data()->mat_.col_sliced(filter);

		this->signal_emitter_->x_data_create_soon(df, soap::VariableType::DataFrame, "Sliced DataFrame");
	}
	
};

void DataFrameItem::row_slice(const Eigen::ArrayX<bool>& filter, bool in_place) {

	if (in_place) {
		this->data()->mat_.row_slice(filter);

		this->__s_update_interface();
	}
	else {

		auto* df = new DataFrame();

		df->mat_ = this->data()->mat_.row_sliced(filter);

		this->signal_emitter_->x_data_create_soon(df, soap::VariableType::DataFrame, "Sliced DataFrame");
	}

};

void DataFrameItem::s_filter_by_feature() {

	G_GETLOCK;
	G_UNLOCK;

	LogicHandler lh(&this->data()->mat_);

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
		return;
	}

	auto filter = lh.resolve(settings[0]);
	bool in_place = switch_to_bool(settings[1]);

	if (filter.count() == 0) {
		G_WARN("No Data Meets Requirements.");
		return;
	}
	if (filter.count() == filter.size()) {
		G_WARN("No Data is Excluded.");
		return;
	}

	this->row_slice(filter, in_place);

	this->__s_update_interface();
};

void DataFrameItem::s_first_row_as_colname() {

	G_GETLOCK;
	G_UNLOCK;

	auto& mat = this->data()->mat_;
	QStringList row1 = mat.get_row(0), column_names = mat.colnames_;
	row1 = custom::make_unique(row1);
	int size = row1.size();
	mat.colnames_ = row1;
	for (int i = 0; i < size; ++i) {
		QString column_name = column_names[i], new_name = row1[i];

		auto ptr = mat.data_[column_name];
		auto type = mat.data_type_[column_name];
		mat.data_.erase(column_name);
		mat.data_type_.erase(column_name);
		mat.data_[new_name] = ptr;
		mat.data_type_[new_name] = type;
		if (type == CustomMatrix::DataType::QStringFactor) {
			auto factor = mat.string_factors_[column_name];
			mat.string_factors_.erase(column_name);
			mat.string_factors_[new_name] = factor;
		}
		else if (type == CustomMatrix::DataType::IntegerFactor) {
			auto int_factor = mat.integer_factors_[column_name];
			mat.integer_factors_.erase(column_name);
			mat.integer_factors_[new_name] = int_factor;
		}
	}

	this->__s_update_interface();
};


void DataFrameItem::s_remove_prefix() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_suffix() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_until_first() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_from_first() {

	G_GETLOCK;

	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_until_last() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_from_last() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_all() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_first() {

	G_GETLOCK;
	G_UNLOCK;

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_add_suffix() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_add_prefix() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_remove_last() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

void DataFrameItem::s_slice_data() {

	G_GETLOCK;
	G_UNLOCK;

	this->__data_delete_soon();

	auto&& metadata = this->data()->mat_;

	auto valid_names = metadata.get_type_names(CustomMatrix::DataType::QStringFactor) << metadata.get_type_names(CustomMatrix::DataType::QString);
	if (valid_names.isEmpty()) {
		G_WARN("No Valid DataFrame.");
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

