#include "LogicHandler.h"

LogicHandler::LogicHandler(SingleCellRna* single_cell_rna) :
	source_type_(LogicHandler::SourceType::SingleCellRna),
	source_(single_cell_rna),
	logic_width_(single_cell_rna->counts()->cols())
{
	auto metadata = single_cell_rna->metadata();

	if (metadata != nullptr) {
		this->process_cm(&metadata->mat_);
	}
	auto counts = single_cell_rna->counts();

	if (counts != nullptr) {
		for (auto&& n : counts->rownames_) {
			this->integer_numeric_names_.insert(n);
		}
	}

	this->finalize();
};

LogicHandler::LogicHandler(SingleCellAtac* single_cell_atac) :
	source_type_(LogicHandler::SourceType::SingleCellAtac),
	source_(single_cell_atac),
	logic_width_(single_cell_atac->counts()->cols())
{
	auto metadata = single_cell_atac->metadata();

	if (metadata != nullptr) {
		this->process_cm(&metadata->mat_);
	}
	auto counts = single_cell_atac->counts();

	if (counts != nullptr) {
		for (auto&& n : counts->rownames_) {
			this->integer_numeric_names_.insert(n);
		}
	}

	this->finalize();
};

LogicHandler::LogicHandler(DataField* data_field) :
	source_type_(LogicHandler::SourceType::DataField),
	source_(data_field),
	logic_width_(data_field->counts()->cols())
{

	auto counts = data_field->counts();

	if (counts != nullptr) {
		for (auto&& n : counts->rownames_) {
			this->integer_numeric_names_.insert(n);
		}
	}

	this->finalize();
};

LogicHandler::LogicHandler(SingleCellMultiome* single_cell_multiome) :
	source_type_(LogicHandler::SourceType::SingleCellMultiome),
	source_(single_cell_multiome),
	logic_width_(single_cell_multiome->rna_counts()->cols())
{
	auto metadata = single_cell_multiome->metadata();

	if (metadata != nullptr) {
		this->process_cm(&metadata->mat_);
	}
	auto counts = single_cell_multiome->rna_counts();

	if (counts != nullptr) {
		for (auto&& n : counts->rownames_) {
			this->integer_numeric_names_.insert(n);
		}
	}

	this->finalize();
};

LogicHandler::LogicHandler(CustomMatrix* mat) :
	source_type_(LogicHandler::SourceType::CustomMatrix),
	source_(mat),
	logic_width_(mat->rows())
{
	this->process_cm(mat);

	this->finalize();
};


Eigen::ArrayX<bool> LogicHandler::resolve_triplet_single_cell_multiome(
	const QString& feature,
	const QString& compare,
	const QString& value) const
{
	Eigen::ArrayX<bool> ret;

	auto single_cell_multiome = static_cast<SingleCellMultiome*>(this->source_);

	auto& metadata = single_cell_multiome->metadata()->mat_;

	if (metadata.contains(feature)) {

		CustomMatrix::DataType datatype = metadata.data_type_[feature];
		if (datatype == CustomMatrix::DataType::QStringFactor || datatype == CustomMatrix::DataType::IntegerFactor) {
			if (compare == "=") {
				ret = custom::equal(metadata.get_qstring(feature), value);
			}
			else {
				ret = custom::not_equal(metadata.get_qstring(feature), value);
			}
		}
		else {
			QVector<double> vec = metadata.get_double(feature);

			if (compare == u"=") {
				ret = custom::equal(vec, value.toDouble());
			}
			else if (compare == u"≠") {
				ret = custom::not_equal(vec, value.toDouble());
			}
			else if (compare == u"≤") {
				ret = custom::less_equal(vec, value.toDouble());
			}
			else if (compare == u"<") {
				ret = custom::less_than(vec, value.toDouble());
			}
			else if (compare == u"≥") {
				ret = custom::greater_equal(vec, value.toDouble());
			}
			else {
				ret = custom::greater_than(vec, value.toDouble());
			}
		}
	}
	else if (single_cell_multiome->rna_counts()->rownames_.contains(feature))
	{
		Eigen::ArrayXd vec = single_cell_multiome->rna_counts()->get_row(feature).cast<double>();

		if (compare == u"=") {
			ret = (vec == value.toDouble());
		}
		else if (compare == u"≠") {
			ret = (vec != value.toDouble());
		}
		else if (compare == u"≤") {
			ret = (vec <= value.toDouble());
		}
		else if (compare == u"<") {
			ret = (vec < value.toDouble());
		}
		else if (compare == u"≥") {
			ret = (vec >= value.toDouble());
		}
		else {
			ret = (vec > value.toDouble());
		}
	}
	else if (single_cell_multiome->atac_counts()->rownames_.contains(feature))
	{
		Eigen::ArrayXd vec = single_cell_multiome->atac_counts()->get_row(feature).cast<double>();

		if (compare == u"=") {
			ret = (vec == value.toDouble());
		}
		else if (compare == u"≠") {
			ret = (vec != value.toDouble());
		}
		else if (compare == u"≤") {
			ret = (vec <= value.toDouble());
		}
		else if (compare == u"<") {
			ret = (vec < value.toDouble());
		}
		else if (compare == u"≥") {
			ret = (vec >= value.toDouble());
		}
		else {
			ret = (vec > value.toDouble());
		}
	}

	return ret;
};

Eigen::ArrayX<bool> LogicHandler::resolve_triplet_cm(
	const QString& feature, 
	const QString& compare, 
	const QString& value) const {

	Eigen::ArrayX<bool> ret;

	auto mat = static_cast<CustomMatrix*>(this->source_);

	if (mat->contains(feature)) {

		CustomMatrix::DataType datatype = mat->data_type_[feature];
		if (datatype == CustomMatrix::DataType::QStringFactor || datatype == CustomMatrix::DataType::IntegerFactor) {
			if (compare == "=") {
				ret = custom::equal(mat->get_qstring(feature), value);
			}
			else {
				ret = custom::not_equal(mat->get_qstring(feature), value);
			}
		}
		else {
			QVector<double> vec = mat->get_double(feature);

			if (compare == u"=") {
				ret = custom::equal(vec, value.toDouble());
			}
			else if (compare == u"≠") {
				ret = custom::not_equal(vec, value.toDouble());
			}
			else if (compare == u"≤") {
				ret = custom::less_equal(vec, value.toDouble());
			}
			else if (compare == u"<") {
				ret = custom::less_than(vec, value.toDouble());
			}
			else if (compare == u"≥") {
				ret = custom::greater_equal(vec, value.toDouble());
			}
			else {
				ret = custom::greater_than(vec, value.toDouble());
			}
		}
	}

	return ret;
};

Eigen::ArrayX<bool> LogicHandler::resolve_triplet_single_cell_rna(
	const QString& feature,
	const QString& compare,
	const QString& value) const
{
	Eigen::ArrayX<bool> ret;

	auto single_cell_rna = static_cast<SingleCellRna*>(this->source_);

	auto& metadata = single_cell_rna->metadata()->mat_;

	if (metadata.contains(feature)) {

		CustomMatrix::DataType datatype = metadata.data_type_[feature];
		if (datatype == CustomMatrix::DataType::QStringFactor || datatype == CustomMatrix::DataType::IntegerFactor) {
			if (compare == "=") {
				ret = custom::equal(metadata.get_qstring(feature), value);
			}
			else {
				ret = custom::not_equal(metadata.get_qstring(feature), value);
			}
		}
		else {
			QVector<double> vec = metadata.get_double(feature);

			if (compare == u"=") {
				ret = custom::equal(vec, value.toDouble());
			}
			else if (compare == u"≠") {
				ret = custom::not_equal(vec, value.toDouble());
			}
			else if (compare == u"≤") {
				ret = custom::less_equal(vec, value.toDouble());
			}
			else if (compare == u"<") {
				ret = custom::less_than(vec, value.toDouble());
			}
			else if (compare == u"≥") {
				ret = custom::greater_equal(vec, value.toDouble());
			}
			else {
				ret = custom::greater_than(vec, value.toDouble());
			}
		}
	}
	else if (single_cell_rna->counts()->rownames_.contains(feature))
	{
		Eigen::ArrayXd vec = single_cell_rna->counts()->get_row(feature).cast<double>();

		if (compare == u"=") {
			ret = (vec == value.toDouble());
		}
		else if (compare == u"≠") {
			ret = (vec != value.toDouble());
		}
		else if (compare == u"≤") {
			ret = (vec <= value.toDouble());
		}
		else if (compare == u"<") {
			ret = (vec < value.toDouble());
		}
		else if (compare == u"≥") {
			ret = (vec >= value.toDouble());
		}
		else {
			ret = (vec > value.toDouble());
		}
	}

	return ret;
};

Eigen::ArrayX<bool> LogicHandler::resolve_triplet_single_cell_atac(
	const QString& feature,
	const QString& compare,
	const QString& value) const
{
	Eigen::ArrayX<bool> ret;

	auto single_cell_atac = static_cast<SingleCellAtac*>(this->source_);

	auto& metadata = single_cell_atac->metadata()->mat_;

	if (metadata.contains(feature)) {

		CustomMatrix::DataType datatype = metadata.data_type_[feature];
		if (datatype == CustomMatrix::DataType::QStringFactor || datatype == CustomMatrix::DataType::IntegerFactor) {
			if (compare == "=") {
				ret = custom::equal(metadata.get_qstring(feature), value);
			}
			else {
				ret = custom::not_equal(metadata.get_qstring(feature), value);
			}
		}
		else {
			QVector<double> vec = metadata.get_double(feature);

			if (compare == u"=") {
				ret = custom::equal(vec, value.toDouble());
			}
			else if (compare == u"≠") {
				ret = custom::not_equal(vec, value.toDouble());
			}
			else if (compare == u"≤") {
				ret = custom::less_equal(vec, value.toDouble());
			}
			else if (compare == u"<") {
				ret = custom::less_than(vec, value.toDouble());
			}
			else if (compare == u"≥") {
				ret = custom::greater_equal(vec, value.toDouble());
			}
			else {
				ret = custom::greater_than(vec, value.toDouble());
			}
		}
	}
	else if (single_cell_atac->counts()->rownames_.contains(feature))
	{
		Eigen::ArrayXd vec = single_cell_atac->counts()->get_row(feature).cast<double>();

		if (compare == u"=") {
			ret = (vec == value.toDouble());
		}
		else if (compare == u"≠") {
			ret = (vec != value.toDouble());
		}
		else if (compare == u"≤") {
			ret = (vec <= value.toDouble());
		}
		else if (compare == u"<") {
			ret = (vec < value.toDouble());
		}
		else if (compare == u"≥") {
			ret = (vec >= value.toDouble());
		}
		else {
			ret = (vec > value.toDouble());
		}
	}

	return ret;
};

Eigen::ArrayX<bool> LogicHandler::resolve_triplet_data_field(
	const QString& feature,
	const QString& compare,
	const QString& value) const
{
	Eigen::ArrayX<bool> ret;

	auto data_field = static_cast<DataField*>(this->source_);

	if (data_field->counts()->rownames_.contains(feature)) {

		Eigen::ArrayXd vec = data_field->counts()->get_row(feature).cast<double>();

		if (compare == u"=") {
			ret = (vec == value.toDouble());
		}
		else if (compare == u"≠") {
			ret = (vec != value.toDouble());
		}
		else if (compare == u"≤") {
			ret = (vec <= value.toDouble());
		}
		else if (compare == u"<") {
			ret = (vec < value.toDouble());
		}
		else if (compare == u"≥") {
			ret = (vec >= value.toDouble());
		}
		else {
			ret = (vec > value.toDouble());
		}
	}

	return ret;
};

Eigen::ArrayX<bool> LogicHandler::resolve_triplet(const QString& feature, const QString& comp, const QString& val) const {

	if (this->source_type_ == SourceType::SingleCellRna) {
		return this->resolve_triplet_single_cell_rna(feature, comp, val);
	}

	if (this->source_type_ == SourceType::SingleCellAtac) {
		return this->resolve_triplet_single_cell_atac(feature, comp, val);
	}

	if (this->source_type_ == SourceType::SingleCellMultiome) {
		return this->resolve_triplet_single_cell_multiome(feature, comp, val);
	}

	if (this->source_type_ == SourceType::DataField) {
		return this->resolve_triplet_data_field(feature, comp, val);
	}

	if (this->source_type_ == SourceType::CustomMatrix) {
		return this->resolve_triplet_cm(feature, comp, val);
	}

	return Eigen::ArrayX<bool>{};
};

Eigen::ArrayX<bool> LogicHandler::resolve(const QString& s) const {

	Eigen::ArrayX<bool> filter;

	auto or_logics = s.split(SOAP_DELIMITER3);

	for (auto&& or_logic : or_logics) {
		auto logics = or_logic.split(SOAP_DELIMITER2);

		Eigen::ArrayX<bool> or_filter;

		for (auto&& logic : logics) {
			auto comps = logic.split(SOAP_DELIMITER);

			if (comps.size() == 3) {

				auto f = this->resolve_triplet(comps[0], comps[1], comps[2]);

				if (f.size() != 0) {
					if (or_filter.size() == 0) /*find finished logic for 1st time*/ {
						or_filter = f;
					}
					else {
						or_filter += f;
					}
				}
			}
		}

		if (or_filter.size() != 0) {
			if (filter.size() != 0) {
				filter *= or_filter;
			}
			else {
				filter = or_filter;
			}
		}
	}

	if (filter.size() == 0) {
		filter.resize(this->logic_width_);
		filter.setConstant(false);
	}

	return filter;
};

std::tuple<bool, LogicHandler::DataType, QStringList> LogicHandler::get_content(const QString& name) const {

	auto iter1 = this->string_factor_contents_.find(name);
	if (iter1 != this->string_factor_contents_.cend()) {
		return { true, LogicHandler::DataType::StringFactor, iter1->second };
	}

	auto iter2 = this->integer_factor_contents_.find(name);
	if (iter2 != this->integer_factor_contents_.cend()) {
		return { true, LogicHandler::DataType::IntegerFactor, iter2->second };
	}

	if (this->integer_numeric_names_.contains(name)) {
		return { true, LogicHandler::DataType::IntegerNumeric, {} };
	}

	if (this->double_numeric_names_.contains(name)) {
		return { true, LogicHandler::DataType::DoubleNumeric, {} };
	}

	return { false, LogicHandler::DataType::NoType, {} };
};

LogicHandler::DataType LogicHandler::get_type(const QString& name) const {

	if (this->string_factor_names_.contains(name)) {
		return LogicHandler::DataType::StringFactor;
	}

	if (this->integer_factor_names_.contains(name)) {
		return LogicHandler::DataType::IntegerFactor;
	}

	if (this->integer_numeric_names_.contains(name)) {
		return LogicHandler::DataType::IntegerNumeric;
	}

	if (this->double_numeric_names_.contains(name)) {
		return LogicHandler::DataType::DoubleNumeric;
	}

	return LogicHandler::DataType::NoType;
};

const QStringList& LogicHandler::data_names() {

	return this->data_names_;
};

void LogicHandler::finalize() {

	int name_size = this->string_factor_names_.size() + this->integer_factor_names_.size() \
		+ this->integer_numeric_names_.size() + this->double_numeric_names_.size();

	this->data_names_.reserve(name_size);

	for (auto&& n : this->string_factor_names_) {
		this->data_names_ << n;
	}

	for (auto&& n : this->integer_factor_names_) {
		this->data_names_ << n;
	}

	for (auto&& n : this->integer_numeric_names_) {
		this->data_names_ << n;
	}

	for (auto&& n : this->double_numeric_names_) {
		this->data_names_ << n;
	}

	this->data_names_ = custom::unique(this->data_names_);
};

void LogicHandler::process_cm(CustomMatrix* mat) {

	for (auto&& [name, type] : mat->data_type_) {
		switch (type)
		{
		case CustomMatrix::DataType::QStringFactor:
			this->string_factor_names_.insert(name);
			this->string_factor_contents_[name] = mat->string_factors_[name];
			break;
		case CustomMatrix::DataType::IntegerFactor:
			this->integer_factor_names_.insert(name);
			this->integer_factor_contents_[name] = custom::cast<QString>(mat->integer_factors_[name]);
			break;
		case CustomMatrix::DataType::QString:
			break;
		case CustomMatrix::DataType::IntegerNumeric:
			this->integer_numeric_names_.insert(name);
			break;
		case CustomMatrix::DataType::DoubleNumeric:
			this->double_numeric_names_.insert(name);
			break;
		default:
			break;
		}
	}
};