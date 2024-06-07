#include "GenomicRange.h"

#include <QDebug>

#include "Custom.h"
#include "GenomeUtility.h"

GenomicRange::GenomicRange(
	const RunLengthEncoding<QString>& sequence_names, 
	const IRange& ranges, 
	const RunLengthEncoding<char>& strand
):
	sequence_names_(sequence_names),
	ranges_(ranges),
	strand_(strand),
	metadata_(sequence_names.size())
{};

void GenomicRange::clear() {
	this->sequence_names_.clear();
	this->ranges_.clear();
	this->strand_.clear();
	this->sequence_information_.clear();
	this->metadata_.clear();
};

QStringList GenomicRange::find_gene(const QString& region) const {

	if (!this->metadata_.contains(METADATA_GENOMIC_RANGE_GENE_NAME, CustomMatrix::DataType::QString)) {
		return {};
	}

	auto [seq_name, start, end, success] = _Cs string_to_peak(region);

	if (!success) {
		return {};
	}
	
	auto sequence_filter = this->sequence_names_ == seq_name;
	if (sequence_filter.count() == 0) {
		return {};
	}

	auto g_start = _Cs sliced(this->ranges_.start_, sequence_filter);
	auto g_end = _Cs sliced(_Cs add(this->ranges_.start_, this->ranges_.width_), sequence_filter);

	QStringList gene_name = _Cs sliced(this->metadata_.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME), sequence_filter);
	QStringList res;
	const qsizetype sequence_size = g_start.size();
	for (qsizetype i = 0; i < sequence_size; ++i) {

		if (start < g_end[i] && end > g_start[i]) {
			res << gene_name[i];
		}

	}

	return _Cs unique(res);
};

bool GenomicRange::is_empty() const {
	return this->sequence_names_.is_empty();
};

void GenomicRange::row_slice(int start, int n) {

	const int original_size = this->sequence_names_.size();

	this->sequence_names_ = this->sequence_names_.sliced(start, n);
	this->ranges_.start_ = this->ranges_.start_.sliced(start, n);
	this->ranges_.width_ = this->ranges_.width_.sliced(start, n);
	this->strand_ = this->strand_.sliced(start, n);

	if (this->metadata_.cols() > 0 && this->metadata_.rows() == original_size) {

		Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(original_size, false);
		filter.segment(start, n) = true;

		this->metadata_.row_slice(filter);
	}
	else {
		this->metadata_.clear();
	}
};

void GenomicRange::extend(int upstream, int downstream, bool from_midpoint) {

	Eigen::ArrayX<bool> direction = this->strand_ == '+' || this->strand_ == '*';

	if (from_midpoint) {
		auto midpoints = _Cs add(this->ranges_.start_, _Cs integer_divide(this->ranges_.width_, 2));
		this->ranges_.start_ = _Cs minus(midpoints, _Cs ifelse<QVector<int>>(direction, upstream, downstream));
		auto new_end = _Cs add(midpoints, _Cs ifelse<QVector<int>>(direction, downstream, upstream));
		this->ranges_.width_ = _Cs minus(new_end, this->ranges_.start_);
	}
	else {
		auto end = this->get_sequence_end();
		this->ranges_.start_ = _Cs minus(this->ranges_.start_, _Cs ifelse<QVector<int>>(direction, upstream, downstream));
		end = _Cs add(end, _Cs ifelse<QVector<int>>(direction, downstream, upstream));
		this->ranges_.width_ = _Cs minus(end, this->ranges_.start_);
	}

	this->finalize();
};

GenomicRange GenomicRange::extended(int upstream, int downstream, bool from_midpoint) const {
	
	GenomicRange res(*this);

	res.extend(upstream, downstream, from_midpoint);

	return res;
};

std::tuple<QString, int, int, char> 
GenomicRange::operator[](qsizetype index) const{
	return std::make_tuple(
		this->sequence_names_[index], 
		this->ranges_.start_[index], 
		this->ranges_.width_[index],
		this->strand_[index]
	);
};

std::tuple<QString, int, int, char> 
GenomicRange::at(qsizetype index) const {

	return std::make_tuple(
		this->sequence_names_[index],
		this->ranges_.start_[index],
		this->ranges_.width_[index] + this->ranges_.start_[index],
		this->strand_[index]
	);
};

int GenomicRange::count_overlap(const GenomicRange& rhs) const {
	auto seq_names = this->sequence_names_.unique();
	auto seq_names2 = rhs.sequence_names_.unique();
	seq_names = _Cs intersect(seq_names, seq_names2);

	int c{ 0 };

	for (auto&& seq_name : seq_names) {

		auto this_filter = this->sequence_names_ == seq_name;
		auto this_index = _Cs which(this_filter);
		auto this_start = _Cs sliced(this->ranges_.start_, this_filter);
		auto this_width = _Cs sliced(this->ranges_.width_, this_filter);
		auto this_end = _Cs add(this_start, this_width);

		auto sub_filter = rhs.sequence_names_ == seq_name;
		auto sub_start = _Cs sliced(rhs.ranges_.start_, sub_filter);
		auto sub_width = _Cs sliced(rhs.ranges_.width_, sub_filter);
		auto sub_end = _Cs add(sub_start, sub_width);

		int n_sub = sub_start.size();
		int n_this = this_index.size();

		int this_ind{ 0 };
		int sub_ind{ 0 };
		while (this_ind < n_this && sub_ind < n_sub) {

			int s1 = this_start[this_ind];
			int e1 = this_end[this_ind];

			int s2 = sub_start[sub_ind];
			int e2 = sub_end[sub_ind];

			if (e1 <= s2) {
				++this_ind;
				continue;
			}

			if (e2 <= s1) {
				++sub_ind;
				continue;
			}

			++c;
			++this_ind;
		}
	}

	return c;
};

void GenomicRange::subtract(const GenomicRange& rhs) {

	auto seq_names = this->sequence_names_.unique(); 
	auto sub_seq_names = rhs.sequence_names_.unique();
	seq_names = _Cs intersect(seq_names, sub_seq_names);

	const int n_range = this->sequence_names_.size();
	Eigen::ArrayX<bool> filter = Eigen::ArrayX<bool>::Constant(n_range, true);

	for (auto&& seq_name : seq_names) {

		auto this_filter = this->sequence_names_ == seq_name;
		auto this_index = _Cs which(this_filter);
		auto this_start = _Cs sliced(this->ranges_.start_, this_filter);
		auto this_width = _Cs sliced(this->ranges_.width_, this_filter);
		auto this_end = _Cs add(this_start, this_width);

		auto sub_filter = rhs.sequence_names_ == seq_name;
		auto sub_start = _Cs sliced(rhs.ranges_.start_, sub_filter);
		auto sub_width = _Cs sliced(rhs.ranges_.width_, sub_filter);
		auto sub_end = _Cs add(sub_start, sub_width);

		int n_sub = sub_start.size();
		int n_this = this_index.size();

		int this_ind{ 0 };
		int sub_ind{ 0 };
		while (this_ind < n_this && sub_ind < n_sub) {

			int s1 = this_start[this_ind];
			int e1 = this_end[this_ind];

			int s2 = sub_start[sub_ind];
			int e2 = sub_end[sub_ind];

			if (e1 <= s2) {
				++this_ind;
				continue;
			}

			if (e2 <= s1) {
				++sub_ind;
				continue;
			}

			filter[this_index[this_ind]] = false;
			++this_ind;
		}
	}

	this->row_slice(filter);
};

void GenomicRange::merge_nearby_range() {

	this->metadata_.clear();

	QVector<int> new_start, new_end;
	RunLengthEncoding<QString> new_sequence_names;
	RunLengthEncoding<char> new_strand;

	auto seq_names = this->sequence_names_.unique();

	for (auto&& seq_name : seq_names) {
		auto filter = this->sequence_names_ == seq_name;
		auto start = _Cs sliced(this->ranges_.start_, filter);
		auto width = _Cs sliced(this->ranges_.width_, filter);
		auto end = _Cs add(start, width);

		int n_range = start.size(), new_n_range{ 0 };
		int now_start = start[0], now_end = end[0];
		for (int i = 1; i < n_range; ++i) {

			if (start[i] <= now_end){
				if (end[i] > now_end) {
					now_end = end[i];
				}
			}
			else {
				new_start << now_start;
				new_end << now_end;
				++new_n_range;

				now_start = start[i];
				now_end = end[i];
			}
		}

		new_start << now_start;
		new_end << now_end;
		++new_n_range;

		new_sequence_names.append(seq_name, new_n_range);
		new_strand.append('*', new_n_range);
	}

	this->sequence_names_ = new_sequence_names;
	this->ranges_.start_ = new_start;
	this->ranges_.width_ = _Cs minus(new_end, new_start);
	this->strand_ = new_strand;
};

void GenomicRange::append2(
	const QString& sequence_name,
	const int start,
	const int end,
	const char strand
) {
	this->sequence_names_ << sequence_name;
	this->ranges_.append(start, end - start);
	this->strand_.append(strand);
};

void GenomicRange::append(
	const QString& sequence_name,
	const int start, 
	const int width, 
	const char strand
) {
	this->sequence_names_ << sequence_name;
	this->ranges_.append(start, width);
	this->strand_.append(strand);
};

void GenomicRange::append(const std::tuple<QString, int, int, char>& segment) {
	this->sequence_names_ << std::get<0>(segment);
	this->ranges_.append(std::get<1>(segment), std::get<2>(segment));
	this->strand_.append(std::get<3>(segment));
};

void GenomicRange::append(const GenomicRange& genomic_range) {
	this->sequence_names_ << genomic_range.sequence_names_;
	this->ranges_.append(genomic_range.ranges_);
	this->strand_ << genomic_range.strand_;
	this->sequence_information_.append(genomic_range.sequence_information_);
	this->metadata_.clear();
	this->finalize();
};

QVector<int> GenomicRange::get_sequence_end() const {
	return _Cs add(this->ranges_.start_, this->ranges_.width_);
};

std::size_t GenomicRange::size() const {

	return this->sequence_names_.size();
};

QString GenomicRange::get_range_name(int row) {
	int start_position = this->ranges_.start_[row], width = this->ranges_.width_[row];
	return this->sequence_names_[row] + ":" + QString::number(start_position) + "-" + QString::number(start_position + width);
};

QStringList GenomicRange::get_range_names() {
	QVector<int>& start_position = this->ranges_.start_;
	QVector<int> end_position = _Cs add(start_position, this->ranges_.width_);

	qsizetype size = start_position.size();

	QStringList range_names(size);
	for (qsizetype i = 0; i < size; ++i) {
		range_names[i] = this->sequence_names_[i] + ":" + QString::number(start_position[i]) + "-" + QString::number(end_position[i]);
	}
	return range_names;
}

QString GenomicRange::get_qstring(int row, int col) const {

	switch (col)
	{
	case 0:
		return this->sequence_names_[row];
	case 1:
		return QString::number(this->ranges_.start_[row]);
	case 2:
		return QString::number(this->ranges_.start_[row] + this->ranges_.width_[row]);
	case 3:
		return QString(QChar(this->strand_[row]));
	default:
		return this->metadata_.get_qstring(row, col - 4);
		break;
	}

};

int GenomicRange::rows() const {
	return this->sequence_names_.size();
};

int GenomicRange::cols() const {
	return this->metadata_.cols() + 4;;
};

QStringList GenomicRange::colnames() const {
	return QStringList() << "Seq Names" << "Range Start" << "Range End" << "Strand" << this->metadata_.colnames_;
};

void GenomicRange::show(int showsize) {
	const std::size_t size = this->size();
	for (std::size_t i = 0; i < size && i < showsize; ++i) {
		qDebug() << this->sequence_names_[i] << this->ranges_.start_[i] << this->ranges_.width_[i] << this->strand_[i];
	}
};

void GenomicRange::finalize() {

	QStringList unique_seqs = this->sequence_names_.unique();

	QVector<int> final_order;
	
	for (const auto& seq : unique_seqs) {

		auto order = this->sequence_names_.match(seq);

		QVector<int> start = _Cs reordered(this->ranges_.start_, order);
		QVector<int> width = _Cs reordered(this->ranges_.width_, order);
		
		auto filter = _Cs greater_than(width, 0);
		if (filter.count() < filter.size()) {
			order = _Cs sliced(order, filter);
			start = _Cs sliced(start, filter);
			// width = _Cs sliced(width, filter);
		}
		
		auto seq_order = _Cs order(start);
		order = _Cs reordered(order, seq_order);
		final_order << order;
	}

	this->row_reorder(final_order);
};