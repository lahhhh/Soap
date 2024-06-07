#pragma once

#include "Identifier.h"

#include <QMap>

#include "CustomTypeTraits.h"
#include "CustomTemplates.h"

class CustomMatrix
{
public:

	CustomMatrix() = default;
	CustomMatrix(int nrow);
	CustomMatrix(const QStringList& rownames);
	CustomMatrix(const CustomMatrix& mat);
	CustomMatrix(CustomMatrix&& mat) noexcept;
	~CustomMatrix();

	CustomMatrix& operator=(const CustomMatrix& mat);
	CustomMatrix& operator=(CustomMatrix&& mat) noexcept;

	enum class DataType : int  { 
		NoType = 0, 
		IntegerNumeric = 1, 
		IntegerFactor = 2, 
		DoubleNumeric = 3, 
		QString = 4, 
		QStringFactor = 5
	};

	QStringList rownames_;
	QStringList colnames_;

	std::map<QString, void* > data_;
	std::map<QString, DataType> data_type_;
	std::map<QString, QStringList> string_factors_;
	std::map<QString, QList<int> > integer_factors_;	

	qsizetype rows() const;
	qsizetype cols() const;

	void* at(const QString& name) const;

	bool contains(const QString& name) const;
	bool contains(const QStringList& names) const;
	bool contains(const QString& name, DataType data_type) const;

	void update(const QString& name, const QStringList& data, DataType = DataType::QString);
	void update(const QString& name, const QVector<int>& data, DataType = DataType::IntegerNumeric);
	void update(const QString& name, const QVector<double>& data, DataType = DataType::DoubleNumeric);

	void to_integer_numeric(const QString& name);
	void to_double_numeric(const QString& name);
	void to_string(const QString& name);
	void to_string_factor(const QString& name);
	void to_integer_factor(const QString& name);

	void remove(const QString& name);
	void reserve(const QStringList& names);


	QStringList get_qstring(const QString& name) const;
	const QStringList& get_const_qstring_reference(const QString& name) const;
	QStringList& get_mutable_qstring_reference(const QString& name);

	QString get_qstring(const QString& name, int index) const;
	QString get_qstring(int row, int col) const;

	QStringList get_row(int index) const;
	QStringList get_row(QString index) const;

	QVector<int> get_integer(const QString& name) const;
	const QVector<int>& get_const_integer_reference(const QString& name) const;
	QVector<int>& get_mutable_integer_reference(const QString& name);

	QVector<double> get_double(const QString& name) const;
	const QVector<double>& get_const_double_reference(const QString& name) const;
	QVector<double>& get_mutable_double_reference(const QString& name);

	template <typename SliceType>
	requires _Cs is_slice_container<SliceType>
	void edit(const QString& feature_name, const SliceType& filter, const QString& value) {

		if (!this->colnames_.contains(feature_name)) {
			return;
		}

		auto type = this->data_type_[feature_name];
		
		if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {

			int val = value.toInt();
			
			QVector<int>* vec = static_cast<QVector<int> *>(this->at(feature_name));			
			_Cs assign(*vec, val, filter);
			
			if (type == DataType::IntegerFactor) {
				this->integer_factors_[feature_name] = _Cs unique(*vec);
			}
		}
		else if (type == DataType::DoubleNumeric) {
			double val = value.toDouble();

			QVector<double>* vec = static_cast<QVector<double> *>(this->at(feature_name));
			_Cs assign(*vec, val, filter);
		}
		else if (type == DataType::QString || type == DataType::QStringFactor) {

			QStringList* vec = static_cast<QStringList*>(this->at(feature_name));
			_Cs assign(*vec, value, filter);
			
			if (type == DataType::QStringFactor) {
				this->string_factors_[feature_name] = _Cs unique(*vec);
			}
		}
	};

	template <typename SliceType>
	requires _Cs is_slice_container<SliceType>
	void row_slice(const SliceType& slice) {
		for (auto& [name, data] : this->data_) {
			DataType type = this->data_type_[name];

			if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
				auto tmp = _Cs sliced_ptr(*static_cast<QVector<int> *>(data), slice);
				delete static_cast<QVector<int>*>(data);
				data = tmp;
				if (type == DataType::IntegerFactor) {
					this->integer_factors_[name] = _Cs unique(*tmp);
				}
			}
			else if (type == DataType::DoubleNumeric) {
				auto tmp = _Cs sliced_ptr(*static_cast<QVector<double> *>(data), slice);
				delete static_cast<QVector<double>*>(data);
				data = tmp;
			}
			else if (type == DataType::QStringFactor || type == DataType::QString) {
				auto tmp = _Cs sliced_ptr(*static_cast<QStringList*>(data), slice);
				delete static_cast<QStringList*>(data);
				data = tmp;
				if (type == DataType::QStringFactor) {
					this->string_factors_[name] = _Cs unique(*tmp);
				}
			}

		}

		this->rownames_ = _Cs sliced(this->rownames_, slice);
	}

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	void col_slice(const SliceType& slice) {

		QStringList removed = _Cs sliced(this->colnames_, _Cs flip(slice));
		
		for (const auto& name : removed) {
			this->remove(name);
		}
	}

	template <typename SliceType1, typename SliceType2>
	requires _Cs is_slice_container<SliceType1> && _Cs is_slice_container<SliceType2>
	void slice(const SliceType1& row_slice, const SliceType2& col_slice) {

		this->row_slice(row_slice);
		
		this->col_slice(col_slice);
	};

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	CustomMatrix row_sliced(const SliceType& slice) const {

		CustomMatrix ret(_Cs sliced(this->rownames_, slice));

		for (auto&& name : this->colnames_) {

			auto type = this->data_type_.at(name);

			if (type == DataType::QStringFactor || type == DataType::QString) {
				ret.update(name, _Cs sliced(*static_cast<QStringList*>(this->at(name)), slice), type);
			}
			else if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
				ret.update(name, _Cs sliced(*static_cast<QVector<int> *>(this->at(name)), slice), type);
			}
			else if (type == DataType::DoubleNumeric) {
				ret.update(name, _Cs sliced(*static_cast<QVector<double> *>(this->at(name)), slice), type);
			}
		}

		return ret;
	};

	template <typename OrderContainer>
	requires _Cs is_order_container<OrderContainer>
	CustomMatrix row_reordered(const OrderContainer& order) const {

		CustomMatrix ret(_Cs reordered(this->rownames_, order));

		for (auto&& name : this->colnames_) {

			auto type = this->data_type_.at(name);
			if (type == DataType::QStringFactor || type == DataType::QString) {
				ret.update(name, _Cs reordered(*static_cast<QStringList*>(this->at(name)), order), type);
			}
			else if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
				ret.update(name, _Cs reordered(*static_cast<QVector<int> *>(this->at(name)), order), type);
			}
			else if (type == DataType::DoubleNumeric) {
				ret.update(name, _Cs reordered(*static_cast<QVector<double> *>(this->at(name)), order), type);
			}
		}

		return ret;
	};

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	CustomMatrix col_sliced(const SliceType& slice) const {

		CustomMatrix ret(this->rownames_);
		
		QStringList remained = _Cs sliced(this->colnames_, slice);
		
		for (const auto& name : remained) {

			auto type = this->data_type_.at(name);

			if (type == DataType::QStringFactor || type == DataType::QString) {
				ret.update(name, *static_cast<QStringList*>(this->at(name)), type);
			}
			else if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
				ret.update(name, *static_cast<QVector<int> *>(this->at(name)), type);
			}
			else if (type == DataType::DoubleNumeric) {
				ret.update(name, *static_cast<QVector<double> *>(this->at(name)), type);
			}
		}

		return ret;
	};

	template <typename SliceType1, typename SliceType2>
		requires _Cs is_slice_container<SliceType1>&& _Cs is_slice_container<SliceType2>
	CustomMatrix sliced(const SliceType1& row_slice, const SliceType2& col_slice) const {
		return this->row_sliced(row_slice).col_sliced(col_slice);
	};

	template <typename OrderContainer>
		requires _Cs is_order_container<OrderContainer>
	void row_reorder(const OrderContainer& order) {

		for (auto& [name, data] : this->data_) {
			DataType type = this->data_type_[name];

			if (type == DataType::IntegerNumeric || type == DataType::IntegerFactor) {
				auto tmp = _Cs reordered_ptr(*static_cast<QVector<int> *>(data), order);
				delete static_cast<QVector<int>*>(data);
				data = tmp;
				if (type == DataType::IntegerFactor) {
					this->integer_factors_[name] = _Cs unique(*tmp);
				}
			}
			else if (type == DataType::DoubleNumeric) {
				auto tmp = _Cs reordered_ptr(*static_cast<QVector<double> *>(data), order);
				delete static_cast<QVector<double>*>(data);
				data = tmp;
			}
			else if (type == DataType::QStringFactor || type == DataType::QString) {
				auto tmp = _Cs reordered_ptr(*static_cast<QStringList*>(data), order);
				delete static_cast<QStringList*>(data);
				data = tmp;
				if (type == DataType::QStringFactor) {
					this->string_factors_[name] = _Cs unique(*tmp);
				}
			}

		}

		this->rownames_ = _Cs reordered(this->rownames_, order);
	};

	CustomMatrix row_reordered(const QStringList& order) const;

	std::pair<QStringList, QList<DataType> > get_type_information() const;
	QMap<QString, QStringList> get_factor_information(bool get_single_factor_data = true) const;
	QStringList get_colnames() const;
	QStringList get_type_names(CustomMatrix::DataType type) const;
	QStringList get_factor_name(bool get_single_factor_data = true) const;

	bool is_empty() const;

	void set_rownames(const QStringList& rownames_);

	void set_nrow(int nrow);

	void clear();

	void reassign_type(const QString& name, DataType type);

	void change_name(const QString& original_name, const QString& new_name);

	void adjust_type();

};
