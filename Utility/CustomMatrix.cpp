#include "CustomMatrix.h"
#include "Custom.h"

CustomMatrix::CustomMatrix(int nrow)
{
	this->rownames_ = custom::cast<QString>(custom::seq_n(0, nrow));
};

CustomMatrix::CustomMatrix(const QStringList& rownames) :
	rownames_(custom::make_unique(rownames))
{};

CustomMatrix::CustomMatrix(const CustomMatrix& mat) :
	rownames_(mat.rownames_),
	colnames_(mat.colnames_)
{
	for (const auto& [name, type] : mat.data_type_) {

		if (type == DataType::IntegerNumeric) {
			this->data_[name] = static_cast<void*>(new QVector<int>(*static_cast<QVector<int> *>(mat.data_.at(name))));
		}
		else if (type == DataType::IntegerFactor) {
			this->data_[name] = static_cast<void*>(new QVector<int>(*static_cast<QVector<int> *>(mat.data_.at(name))));
			this->integer_factors_[name] = mat.integer_factors_.at(name);
		}
		else if (type == DataType::DoubleNumeric) {
			this->data_[name] = static_cast<void*>(new QVector<double>(*static_cast<QVector<double> *>(mat.data_.at(name))));
		}
		else if (type == DataType::QStringFactor) {
			this->data_[name] = static_cast<void*>(new QStringList(*static_cast<QStringList*>(mat.data_.at(name))));
			this->string_factors_[name] = mat.string_factors_.at(name);
		}
		else if (type == DataType::QString) {
			this->data_[name] = static_cast<void*>(new QStringList(*static_cast<QStringList*>(mat.data_.at(name))));
		}

		this->data_type_[name] = type;
	}
};

CustomMatrix& CustomMatrix::operator=(CustomMatrix&& mat) noexcept {

	for (auto& [name, data] : this->data_) {
		auto type = this->data_type_[name];

		if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
			delete static_cast<QVector<int>*>(data);
		}
		else if (type == DataType::DoubleNumeric) {
			delete static_cast<QVector<double>*>(data);
		}
		else if (type == DataType::QStringFactor || type == DataType::QString) {
			delete static_cast<QStringList*>(data);
		}
	}

	this->rownames_ = mat.rownames_;
	this->colnames_ = mat.colnames_;
	this->data_ = mat.data_;
	this->data_type_ = mat.data_type_;
	this->string_factors_ = mat.string_factors_;
	this->integer_factors_ = mat.integer_factors_;

	mat.data_.clear();
	mat.data_type_.clear();
	mat.string_factors_.clear();
	mat.integer_factors_.clear();
	mat.rownames_.clear();
	mat.colnames_.clear();

	return *this;
};

CustomMatrix& CustomMatrix::operator=(const CustomMatrix& mat) {

	for (auto&& [name, data] : this->data_) {
		auto type = this->data_type_[name];

		if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
			delete static_cast<QVector<int>*>(data);
		}
		else if (type == DataType::DoubleNumeric) {
			delete static_cast<QVector<double>*>(data);
		}
		else if (type == DataType::QStringFactor || type == DataType::QString) {
			delete static_cast<QStringList*>(data);
		}
	}

	this->rownames_ = mat.rownames_;
	this->colnames_ = mat.colnames_;
	this->data_.clear();
	this->data_type_.clear();
	this->string_factors_.clear();
	this->integer_factors_.clear();

	for (const auto& [name, type] : mat.data_type_) {

		if (type == DataType::IntegerNumeric) {
			this->data_[name] = static_cast<void*>(new QVector<int>(*static_cast<QVector<int> *>(mat.data_.at(name))));
		}
		else if (type == DataType::IntegerFactor) {
			this->data_[name] = static_cast<void*>(new QVector<int>(*static_cast<QVector<int> *>(mat.data_.at(name))));
			this->integer_factors_[name] = mat.integer_factors_.at(name);
		}
		else if (type == DataType::DoubleNumeric) {
			this->data_[name] = static_cast<void*>(new QVector<double>(*static_cast<QVector<double> *>(mat.data_.at(name))));
		}
		else if (type == DataType::QStringFactor) {
			this->data_[name] = static_cast<void*>(new QStringList(*static_cast<QStringList*>(mat.data_.at(name))));
			this->string_factors_[name] = mat.string_factors_.at(name);
		}
		else if (type == DataType::QString) {
			this->data_[name] = static_cast<void*>(new QStringList(*static_cast<QStringList*>(mat.data_.at(name))));
		}

		this->data_type_[name] = type;
	}

	return *this;
}

CustomMatrix::CustomMatrix(CustomMatrix&& mat) noexcept :
	rownames_(mat.rownames_),
	colnames_(mat.colnames_),
	data_(mat.data_),
	data_type_(mat.data_type_),
	string_factors_(mat.string_factors_),
	integer_factors_(mat.integer_factors_)
{

	mat.data_.clear();
	mat.data_type_.clear();
	mat.string_factors_.clear();
	mat.integer_factors_.clear();
	mat.rownames_.clear();
	mat.colnames_.clear();
};

CustomMatrix::~CustomMatrix() {

	this->clear();
}

void CustomMatrix::adjust_type() {

	for (const auto& data_name : this->colnames_) {

		if (this->data_type_[data_name] == DataType::QString) {

			const QStringList& data = this->get_const_qstring_reference(data_name);

			int type = custom::detect_type_mt(data);

			if (type == 2) { // integer
				const qsizetype size = data.size();
				QVector<int>* new_data = new QVector<int>(size);
				for (qsizetype i = 0; i < size; ++i) {
					new_data->operator[](i) = data[i].toInt();
				}

				delete (QStringList*)this->data_[data_name];

				this->data_[data_name] = new_data;

				this->data_type_[data_name] = DataType::IntegerNumeric;
			}

			if (type == 1) {

				const qsizetype size = data.size();

				QVector<double>* new_data = new QVector<double>(size);
				for (qsizetype i = 0; i < size; ++i) {
					new_data->operator[](i) = data[i].toDouble();
				}

				delete (QStringList*)this->data_[data_name];
				this->data_[data_name] = new_data;

				this->data_type_[data_name] = DataType::DoubleNumeric;
			}

			if (type == 0) {

				QStringList factors = custom::unique(data);
				int n_unique = factors.size();
				int size = data.size();

				if (n_unique < std::sqrt(size)) {
					this->data_type_[data_name] = DataType::QStringFactor;
					this->string_factors_[data_name] = factors;
				}
			}
		}
	}
}

void CustomMatrix::clear() {

	for (auto& [name, data] : this->data_) {

		auto type = this->data_type_[name];
		if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
			delete static_cast<QVector<int>*>(data);
		}
		else if (type == DataType::DoubleNumeric) {
			delete static_cast<QVector<double>*>(data);
		}
		else if (type == DataType::QStringFactor || type == DataType::QString) {
			delete static_cast<QStringList*>(data);
		}
	}

	this->data_.clear();
	this->data_type_.clear();
	this->string_factors_.clear();
	this->integer_factors_.clear();
	this->rownames_.clear();
	this->colnames_.clear();
};

void CustomMatrix::to_integer_numeric(const QString& name) {

	if (!this->contains(name))return;

	auto data_type = this->data_type_[name];

	if (data_type == DataType::IntegerNumeric) {
		return;
	}
	else if (data_type == DataType::DoubleNumeric) {
		QVector<double>* vec = static_cast<QVector<double> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](double d) {return (int)d; }), DataType::IntegerNumeric);
	}
	else if (data_type == DataType::IntegerFactor) {
		this->data_type_[name] = DataType::IntegerNumeric;
		this->integer_factors_.erase(name);
	}
	else if (data_type == DataType::QString || data_type == DataType::QStringFactor) {
		QStringList* vec = static_cast<QStringList*>(this->at(name));
		this->update(name, custom::sapply(*vec, [](const QString& str) {return str.toInt(); }), DataType::IntegerNumeric);
	}
};

void CustomMatrix::to_integer_factor(const QString& name) {

	if (!this->contains(name))return;

	auto data_type = this->data_type_[name];

	if (data_type == DataType::IntegerFactor) {
		return;
	}
	else if (data_type == DataType::DoubleNumeric) {
		QVector<double>* vec = static_cast<QVector<double> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](double d) {return (int)d; }), DataType::IntegerFactor);
	}
	else if (data_type == DataType::IntegerNumeric) {
		this->data_type_[name] = DataType::IntegerFactor;
		this->integer_factors_[name] = custom::unique(this->get_const_integer_reference(name));
	}
	else if (data_type == DataType::QString || data_type == DataType::QStringFactor) {
		QStringList* vec = static_cast<QStringList*>(this->at(name));
		this->update(name, custom::sapply(*vec, [](const QString& str) {return str.toInt(); }), DataType::IntegerFactor);
	}
};

void CustomMatrix::to_double_numeric(const QString& name) {

	if (!this->contains(name))return;

	auto data_type = this->data_type_[name];

	if (data_type == DataType::DoubleNumeric) {
		return;
	}
	else if (data_type == DataType::IntegerNumeric || data_type == DataType::IntegerFactor) {
		QVector<int>* vec = static_cast<QVector<int> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](int i) {return (double)i; }), DataType::DoubleNumeric);
	}
	else if (data_type == DataType::QString || data_type == DataType::QStringFactor) {
		QStringList* vec = static_cast<QStringList*>(this->at(name));
		this->update(name, custom::sapply(*vec, [](const QString& str) {return str.toDouble(); }), DataType::DoubleNumeric);
	}
};

void CustomMatrix::to_string(const QString& name) {

	if (!this->contains(name))return;

	auto data_type = this->data_type_[name];

	if (data_type == DataType::QString) {
		return;
	}
	else if (data_type == DataType::IntegerNumeric || data_type == DataType::IntegerFactor) {
		QVector<int>* vec = static_cast<QVector<int> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](int i) {return QString::number(i); }), DataType::QString);
	}
	else if (data_type == DataType::DoubleNumeric) {
		QVector<double>* vec = static_cast<QVector<double> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](double d) {return QString::number(d); }), DataType::QString);
	}
	else if (data_type == DataType::QStringFactor) {
		this->data_type_[name] = DataType::QString;
		this->string_factors_.erase(name);
	}
};

void CustomMatrix::to_string_factor(const QString& name) {

	if (!this->contains(name))return;

	auto data_type = this->data_type_[name];

	if (data_type == DataType::QStringFactor) {
		return;
	}
	else if (data_type == DataType::IntegerNumeric || data_type == DataType::IntegerFactor) {
		QVector<int>* vec = static_cast<QVector<int> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](int i) {return QString::number(i); }), DataType::QStringFactor);
	}
	else if (data_type == DataType::DoubleNumeric) {
		QVector<double>* vec = static_cast<QVector<double> *>(this->at(name));
		this->update(name, custom::sapply(*vec, [](double d) {return QString::number(d); }), DataType::QStringFactor);
	}
	else if (data_type == DataType::QString) {
		this->data_type_[name] = DataType::QStringFactor;
		this->string_factors_[name] = custom::unique(this->get_const_qstring_reference(name));
	}
};

void CustomMatrix::change_name(const QString& original_name, const QString& new_name) {

	if (this->contains(new_name))return;

	int index = this->colnames_.indexOf(original_name);
	if (index == -1)return;

	DataType type = this->data_type_[original_name];

	this->colnames_[index] = new_name;

	this->data_type_[new_name] = type;
	this->data_type_.erase(original_name);

	this->data_[new_name] = this->data_[original_name];
	this->data_.erase(original_name);

	if (type == DataType::IntegerFactor) {
		this->integer_factors_[new_name] = this->integer_factors_[original_name];
		this->integer_factors_.erase(original_name);
	}
	else if (type == DataType::QStringFactor) {
		this->string_factors_[new_name] = this->string_factors_[original_name];
		this->string_factors_.erase(original_name);
	}
};

void CustomMatrix::reassign_type(const QString& name, DataType type) {

	switch (type)
	{
	case CustomMatrix::DataType::NoType:
		break;
	case CustomMatrix::DataType::IntegerNumeric:
		this->to_integer_numeric(name);
		break;
	case CustomMatrix::DataType::DoubleNumeric:
		this->to_double_numeric(name);
		break;
	case CustomMatrix::DataType::QString:
		this->to_string(name);
		break;
	case CustomMatrix::DataType::QStringFactor:
		this->to_string_factor(name);
		break;
	case CustomMatrix::DataType::IntegerFactor:
		this->to_integer_factor(name);
		break;
	default:
		break;
	}
};

qsizetype CustomMatrix::rows() const {

	return this->rownames_.size();
};

qsizetype CustomMatrix::cols() const {

	return this->colnames_.size();
};

void* CustomMatrix::at(const QString& name) const {

	return this->data_.at(name);
}

bool CustomMatrix::contains(const QString& name, DataType data_type) const {

	auto iter = this->data_type_.find(name);

	if (iter != this->data_type_.cend()) {
		if (data_type == iter->second) {
			return true;
		}
	}

	return false;
};

bool CustomMatrix::contains(const QString& name) const {

	return this->data_.contains(name);
};

bool CustomMatrix::contains(const QStringList& names) const {

	for (auto&& name : names) {

		if (!this->data_.contains(name)) {
			return false;
		}
	}

	return true;
};

void CustomMatrix::update(const QString& name, const QStringList& data, DataType data_type) {

	this->remove(name);

	if (this->rownames_.isEmpty()) {
		this->set_nrow(data.size());
	}

	this->data_[name] = new QStringList(data);
	this->data_type_[name] = data_type;

	if (data_type == DataType::QStringFactor) {
		this->string_factors_[name] = custom::unique(data);
	}

	this->colnames_ << name;
};

void CustomMatrix::update(const QString& name, const QVector<int>& data, DataType data_type) {

	this->remove(name);

	if (this->rownames_.isEmpty()) {
		this->set_nrow(data.size());
	}

	this->data_[name] = new QVector<int>(data);
	this->data_type_[name] = data_type;

	if (data_type == DataType::IntegerFactor) {
		this->integer_factors_[name] = custom::unique(data);
	}

	this->colnames_ << name;
};

void CustomMatrix::update(const QString& name, const QVector<double>& data, DataType data_type) {

	this->remove(name);

	if (this->rownames_.isEmpty()) {
		this->set_nrow(data.size());
	}

	this->data_[name] = new QVector<double>(data);
	this->data_type_[name] = data_type;

	this->colnames_ << name;
};

void CustomMatrix::reserve(const QStringList& names) {

	QStringList colnames = this->colnames_;

	for (auto&& colname : colnames) {
		if (!names.contains(colname)) {
			this->remove(colname);
		}
	}
};

void CustomMatrix::remove(const QString& name) {

	if (this->contains(name)) {

		void* data_ptr = this->data_[name];

		if (this->data_type_[name] == DataType::IntegerNumeric) {
			delete static_cast<QVector<int>*>(data_ptr);
		}
		else if (this->data_type_[name] == DataType::DoubleNumeric) {
			delete static_cast<QVector<double>*>(data_ptr);
		}
		else if (this->data_type_[name] == DataType::QString) {
			delete static_cast<QStringList*>(data_ptr);
		}
		else if (this->data_type_[name] == DataType::QStringFactor) {
			delete static_cast<QStringList*>(data_ptr);
			this->string_factors_.erase(name);
		}
		else if (this->data_type_[name] == DataType::IntegerFactor) {
			delete static_cast<QVector<int>*>(data_ptr);
			this->integer_factors_.erase(name);
		}

		this->data_.erase(name);
		this->data_type_.erase(name);
		this->colnames_.removeOne(name);

	}
};

QStringList CustomMatrix::get_qstring(const QString& name) const {

	if (!this->data_.contains(name)) {
		return {};
	}

	auto type = this->data_type_.at(name);
	auto ptr = this->data_.at(name);

	if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
		return custom::cast<QString>(*static_cast<QVector<int> *>(ptr));
	}
	else if (type == DataType::DoubleNumeric) {
		return custom::cast<QString>(*static_cast<QVector<double> *>(ptr));
	}
	else if (type == DataType::QStringFactor || type == DataType::QString) {
		return *static_cast<QStringList*>(ptr);
	}
	else {
		return {};
	}
}

QStringList& CustomMatrix::get_mutable_qstring_reference(const QString& name) {

	return *static_cast<QStringList*>(this->data_[name]);
};

QVector<int>& CustomMatrix::get_mutable_integer_reference(const QString& name) {

	return *static_cast<QVector<int>*>(this->data_[name]);
};

QVector<double>& CustomMatrix::get_mutable_double_reference(const QString& name) {

	return *static_cast<QVector<double>*>(this->data_[name]);
};

const QStringList& CustomMatrix::get_const_qstring_reference(const QString& name) const {

	return *static_cast<QStringList*>(this->data_.at(name));
};

const QVector<int>& CustomMatrix::get_const_integer_reference(const QString& name) const {

	return *static_cast<QVector<int>*>(this->data_.at(name));
};

const QVector<double>& CustomMatrix::get_const_double_reference(const QString& name) const {

	return *static_cast<QVector<double>*>(this->data_.at(name));
};

QString CustomMatrix::get_qstring(const QString& name, int index) const {

	if (!this->data_.contains(name)) {
		return QString();
	}

	auto type = this->data_type_.at(name);
	auto ptr = this->data_.at(name);

	if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
		return QString::number(static_cast<QVector<int> *>(ptr)->operator[](index));
	}
	else if (type == DataType::DoubleNumeric) {
		return QString::number(static_cast<QVector<double> *>(ptr)->operator[](index));
	}
	else if (type == DataType::QStringFactor || type == DataType::QString) {
		return static_cast<QStringList*>(ptr)->operator[](index);
	}
	else {
		return QString();
	}
};

QString CustomMatrix::get_qstring(int row, int col) const {

	return this->get_qstring(this->colnames_[col], row);
};

QStringList CustomMatrix::get_row(int index) const {

	if (index < 0 || index > this->rows())return QStringList();

	QStringList ret;

	for (const auto& name : this->colnames_) {

		auto type = this->data_type_.at(name);

		if (type == DataType::QStringFactor || type == DataType::QString) {
			ret << this->get_const_qstring_reference(name).at(index);
		}
		else if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
			ret << QString::number(this->get_const_integer_reference(name).at(index));
		}
		else if (type == DataType::DoubleNumeric) {
			ret << QString::number(this->get_const_double_reference(name).at(index));
		}
	}
	return ret;
};

QStringList CustomMatrix::get_row(QString index) const {
	return this->get_row(this->rownames_.indexOf(index));
};

QVector<int> CustomMatrix::get_integer(const QString& name) const {

	if (!this->data_.contains(name)) {
		return QVector<int>();
	}

	auto type = this->data_type_.at(name);

	if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
		return *static_cast<QVector<int> *>(this->data_.at(name));
	}
	else if (type == DataType::DoubleNumeric) {
		return custom::cast<int>(*static_cast<QVector<double> *>(this->data_.at(name)));
	}
	else if (type == DataType::QStringFactor || type == DataType::QString) {
		return custom::cast<int>(*static_cast<QStringList*>(this->data_.at(name)));
	}
	else {
		return QVector<int>();
	}
}

QVector<double> CustomMatrix::get_double(const QString& name) const {

	if (!this->data_.contains(name)) {
		return QVector<double>();
	}

	auto type = this->data_type_.at(name);
	if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
		return custom::cast<double>(*static_cast<QVector<int> *>(this->data_.at(name)));
	}
	else if (type == DataType::DoubleNumeric) {
		return *static_cast<QVector<double> *>(this->data_.at(name));
	}
	else if (type == DataType::QStringFactor || type == DataType::QString) {
		return custom::cast<double>(*static_cast<QStringList*>(this->data_.at(name)));
	}
	else {
		return QVector<double>();
	}
}

CustomMatrix CustomMatrix::row_reordered(const QStringList& order) const {
	return this->row_reordered(custom::index_of(order, this->rownames_));
};

std::pair<QStringList, QList<CustomMatrix::DataType>> CustomMatrix::get_type_information() const {
	QStringList names;
	QList<DataType> types;

	for (const auto& [name, type] : this->data_type_) {
		names << name;
		types << type;
	}
	return std::make_pair(names, types);
};

QMap<QString, QStringList> CustomMatrix::get_factor_information(bool get_single_factor_data) const {
	QMap<QString, QStringList> ret;

	for (const auto& [name, type] : this->data_type_) {
		if (type == DataType::QStringFactor) {
			if (!get_single_factor_data) {
				if (this->string_factors_.at(name).size() == 1) {
					continue;
				}
			}
			ret[name] = this->string_factors_.at(name);
		}
		else if (type == DataType::IntegerFactor) {
			if (!get_single_factor_data) {
				if (this->integer_factors_.at(name).size() == 1) {
					continue;
				}
			}
			ret[name] = custom::cast<QString>(this->integer_factors_.at(name));
		}
	}

	return ret;
};

QStringList CustomMatrix::get_factor_name(bool get_single_factor_data) const {

	QStringList ret;

	for (const auto& [name, factor] : this->string_factors_) {
		if (!get_single_factor_data) {
			if (this->string_factors_.at(name).size() == 1) {
				continue;
			}
		}
		ret << name;
	}

	for (const auto& [name, factor] : this->integer_factors_) {
		if (!get_single_factor_data) {
			if (this->integer_factors_.at(name).size() == 1) {
				continue;
			}
		}
		ret << name;
	}

	return ret;
};

QStringList CustomMatrix::get_type_names(CustomMatrix::DataType type) const {

	QStringList ret;

	for (const auto& [name, dtype] : this->data_type_) {
		if (dtype == type) {
			ret << name;
		}
	}

	return ret;
};

QStringList CustomMatrix::get_colnames() const {

	return this->colnames_;
};

bool CustomMatrix::is_empty() const {

	return this->rownames_.isEmpty() || this->colnames_.isEmpty();
};

void CustomMatrix::set_rownames(const QStringList& rownames) {

	this->rownames_ = custom::make_unique(rownames);
};

void CustomMatrix::set_nrow(int nrow) {

	this->rownames_ = custom::cast<QString>(custom::seq_n(0, nrow));
};