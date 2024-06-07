#pragma once

#include "Identifier.h"

#include <QVector>

#include "Custom.h"

template <typename T>
class RunLengthEncoding
{
public:
	RunLengthEncoding() = default;

	RunLengthEncoding(const RunLengthEncoding& rle) :
		data_(rle.data_),
		index_(rle.index_),
		size_(rle.size_),
		encoding_size_(rle.encoding_size_)
	{};

	template<typename ContainerType>
		requires _Cs is_container<ContainerType>
	RunLengthEncoding(const ContainerType& vec) {

		auto iter = std::cbegin(vec);
		auto end = std::cend(vec);

		for (; iter != end; ++iter) {
			this->append(*iter);
		}
	}

	QVector<T> data_;

	QVector<qsizetype> index_; // index of segments

	qsizetype size_ = 0; // length 

	qsizetype encoding_size_ = 0; // length of segments

	bool is_empty() const {
		return this->size_ == 0;
	}

	void clear() {
		this->data_.clear();
		this->index_.clear();
		this->size_ = 0;
		this->encoding_size_ = 0;
	};

	QVector<T> to_qvector() const {
		QVector<T> result(this->size_);

		for (qsizetype i = 0; i < this->encoding_size_; ++i) {
			const T& value = this->data_[i];
			qsizetype begin = this->index_[i];
			qsizetype end = this->index_[i + 1];
			for (qsizetype j = begin; j < end; ++j) {
				result[j] = value;
			}
		}

		return result;
	}

	template<typename Comp>
	Eigen::ArrayX<bool> compare(const T& val, Comp op) const {

		Eigen::ArrayX<bool> result(this->size_);

		for (qsizetype i = 0; i < this->encoding_size_; ++i) {

			result.segment(this->index_[i], this->index_[i + 1] - this->index_[i]) = op(this->data_[i], val);
		}

		return result;
	}

	Eigen::ArrayX<bool> operator==(const T& val) const {

		return this->compare(val, std::equal_to{});
	}

	Eigen::ArrayX<bool> operator!=(const T& val) const {

		return this->compare(val, std::not_equal_to{});
	}

	Eigen::ArrayX<bool> operator<=(const T& val) const {

		return this->compare(val, std::less_equal{});
	}

	Eigen::ArrayX<bool> operator<(const T& val) const {

		return this->compare(val, std::less{});
	}

	Eigen::ArrayX<bool> operator>=(const T& val) const {

		return this->compare(val, std::greater_equal{});
	}

	Eigen::ArrayX<bool> operator>(const T& val) const {

		return this->compare(val, std::greater{});
	}

	template<typename ContainerType, typename ElementType = _Cs get_element_raw_type<ContainerType>>
	requires _Cs is_container<ContainerType>
	static RunLengthEncoding<ElementType> from_container(const ContainerType& vec) {

		RunLengthEncoding<ElementType> ret;

		auto iter = std::cbegin(vec);
		auto end = std::cend(vec);

		for (; iter != end; ++iter) {
			ret.append(*iter);
		}

		return ret;
	}

	// TODO : optimize
	template<typename S>
	requires _Cs is_slice_container<S>
	RunLengthEncoding<T> sliced(const S& slice) const{

		return RunLengthEncoding<T>::from_container(_Cs sliced(this->to_qvector(), slice));
	};

	// TODO : optimize
	RunLengthEncoding<T> sliced(int start, int n) const {

		return RunLengthEncoding<T>::from_container(this->to_qvector().sliced(start, n));
	};

	template <typename Order>
		requires _Cs is_order_container<Order>
	RunLengthEncoding<T> reordered(const Order& order) const {

		return RunLengthEncoding<T>::from_container(_Cs reordered(this->to_qvector(), order));
	};

	template <typename Order>
		requires _Cs is_order_container<Order>
	void reorder(const Order& order) {

		*this = this->reordered(order);
	};

	void assign(const T& val, const Eigen::ArrayX<bool>& index) {
		QVector<T> vec = this->to_qvector();
		const Eigen::Index index_size = index.size();
		for (Eigen::Index i = 0; i < index_size; ++i) {
			if (index[i]) {
				vec[i] = val;
			}
		}
		*this = RunLengthEncoding<T>::from_container(vec);
	};

	QVector<int> match(const T& val) const {

		return _Cs which(*this == val);
	};

	[[nodiscard]] QVector<T> unique() const{

		return _Cs unique(this->data_);
	}

	qsizetype size() const {

		return this->size_;
	}

	RunLengthEncoding& append(const T& val) {
		if (this->size_ == 0) {
			this->data_ << val;
			this->index_ << 0 << 1;
			this->size_ = 1;
			this->encoding_size_ = 1;
			return *this;
		}
		if (val == this->data_[this->encoding_size_ - 1]) {
			++this->index_[this->encoding_size_];
		}
		else {
			this->data_ << val;
			this->index_ << this->index_[this->encoding_size_] + 1;
			++this->encoding_size_;
		}
		++this->size_;
		return *this;
	}

	RunLengthEncoding& append(const T& val, const qsizetype n) {
		if (this->size_ == 0) {
			this->data_ << val;
			this->index_ << 0 << n;
			this->size_ = n;
			this->encoding_size_ = 1;
			return *this;
		}
		if (val == this->data_[this->encoding_size_ - 1]) {
			this->index_[this->encoding_size_] += n;
		}
		else {
			this->data_ << val;
			this->index_ << this->index_[this->encoding_size_] + n;
			++this->encoding_size_;
		}
		this->size_ += n;
		return *this;
	}

	RunLengthEncoding& append(const RunLengthEncoding<T>& rle) {
		if (rle.size() == 0) {
			return *this;
		}
		if (this->size_ == 0) {
			this->data_ = rle.data_;
			this->index_ = rle.index_;
			this->size_ = rle.size_;
			this->encoding_size_ = rle.encoding_size_;

			return *this;
		}
		if (this->last() == rle.last()) {
			if (rle.encoding_size_ == 1) {
				this->index_[this->encoding_size_] += rle.size();
				this->size_ += rle.size();
			}
			else {
				this->data_ << rle.data_.sliced(1);
				this->index_[this->encoding_size_] += rle.index_[1];
				this->index_ << _Cs add(rle.index_.sliced(2), this->size());
				this->encoding_size_ += (rle.encoding_size_ - 1);
				this->size_ += rle.size();
			}
		}
		else {
			this->data_ << rle.data_;
			this->index_ << _Cs add(rle.index_.sliced(1), this->size());
			this->encoding_size_ += rle.encoding_size_;
			this->size_ += rle.size();
		}
		return *this;
	}

	RunLengthEncoding& operator<<(const T& val) {
		return this->append(val);
	}

	RunLengthEncoding& operator<<(const RunLengthEncoding<T>& rle) {
		return this->append(rle);
	}

	T operator[](const qsizetype loc) const{

		auto upper_bound_it = std::ranges::upper_bound(this->index_, loc);

		int loc1 = std::distance(this->index_.begin(), upper_bound_it) - 1;

		return this->data_[loc1];
	}

	const T& last() const{
		return this->data_.last();
	}
};

