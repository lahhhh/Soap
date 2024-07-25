#pragma once

#include "Identifier.h"

#include "RunLengthEncoding.h"
#include "IRange.h"
#include "SeqInfo.h"
#include "CustomMatrix.h"

class GenomicRange
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(GenomicRange);

	GenomicRange(
		const RunLengthEncoding<QString>& sequence_names, 
		const IRange& ranges, 
		const RunLengthEncoding<char>& strand
	);

	DataType data_type_{ DataType::Plain };

	RunLengthEncoding<QString> sequence_names_;

	IRange ranges_;

	RunLengthEncoding<char> strand_;

	SeqInfo sequence_information_;

	CustomMatrix metadata_;

	std::size_t size() const;

	int rows() const;
	int cols() const;

	QString get_qstring(int row, int col) const;

	QStringList colnames() const;

	QString get_range_name(int row);
	QStringList get_range_names();

	QVector<int> get_sequence_end() const;

	void finalize();

	template <typename OrderType>
	requires _Cs is_order_container<OrderType>
	void row_reorder(const OrderType& order) {

		const qsizetype order_size = std::size(order), original_size = this->sequence_names_.size();

		this->sequence_names_.reorder(order);
		this->ranges_.start_ = _Cs reordered(this->ranges_.start_, order);
		this->ranges_.width_ = _Cs reordered(this->ranges_.width_, order);
		this->strand_.reorder(order);

		if (this->metadata_.rows() == original_size && this->metadata_.cols() > 0) {
			this->metadata_.row_reorder(order);
		}
		else {
			this->metadata_.clear();
		}
	};

	template <typename SliceType>
	requires _Cs is_slice_container<SliceType>
	void row_slice(const SliceType& slice) {

		const int original_size = this->sequence_names_.size();

		this->sequence_names_ = this->sequence_names_.sliced(slice);
		this->ranges_.start_ = _Cs sliced(this->ranges_.start_, slice);
		this->ranges_.width_ = _Cs sliced(this->ranges_.width_, slice);
		this->strand_ = this->strand_.sliced(slice);
		
		if (this->metadata_.cols() > 0 && this->metadata_.rows() == original_size) {
			this->metadata_.row_slice(slice);
		}
		else {
			this->metadata_.clear();
		}
	};

	template <typename SliceType>
		requires _Cs is_slice_container<SliceType>
	[[nodiscard]] GenomicRange row_sliced(const SliceType& slice) const{

		const int original_size = this->sequence_names_.size();

		GenomicRange ret;

		ret.sequence_names_ = this->sequence_names_.sliced(slice);
		ret.ranges_.start_ = _Cs sliced(this->ranges_.start_, slice);
		ret.ranges_.width_ = _Cs sliced(this->ranges_.width_, slice);
		ret.strand_ = this->strand_.sliced(slice);

		if (this->metadata_.cols() > 0 && this->metadata_.rows() == original_size) {
			ret.metadata_ = this->metadata_.row_sliced(slice);
		}

		return ret;
	};

	template <typename OrderType>
		requires _Cs is_order_container<OrderType>
	GenomicRange row_reordered(const OrderType& order) const {

		GenomicRange new_genomic_range(*this);

		new_genomic_range.row_reorder(order);

		return new_genomic_range;
	};

	void row_slice(int start, int length);

	void append(const GenomicRange& genomic_range);
	void append(const QString& sequence_name, const int start, const int width, const char strand);
	void append2(const QString& sequence_name, const int start, const int end, const char strand);
	void append(const std::tuple<QString, int, int, char>& segment);

	void extend(int upstream, int downstream, bool from_midpoint = false);
	GenomicRange extended(int upstream, int downstream, bool from_midpoint = false) const;

	// it returns width, not end
	std::tuple<QString, int, int, char> operator[](qsizetype index) const;

	std::tuple<QString, int, int, char> at(qsizetype index) const;

	void merge_nearby_range();
	
	void subtract(const GenomicRange& rhs);
	int count_overlap(const GenomicRange& rhs) const;

	QStringList find_gene(const QString& region) const;
	QStringList find_chromosome_gene(const QString& chromosome) const;

	bool is_empty() const;

	void clear();

	void show(int showsize = 100000);

	G_SET_IDENTIFIER("GenomicRange");

};

