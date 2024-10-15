#pragma once

#include <QStringList>

#include "Custom.h"
#include "GenomeUtility.h"

struct ChromIndex {

	std::vector<int> start;

	std::vector<int> end;

	std::vector<int> id;

	void append(int start, int end, int id) {
		this->start.push_back(start);
		this->end.push_back(end);
		this->id.push_back(id);
	}

	void append(int start, int end) {
		this->start.push_back(start);
		this->end.push_back(end);
	}
};

enum class IndexMode : int { ID, Position, Strand, IdStrand };


/*
	Note : do not handle overlapped ranges
*/
template <IndexMode Mode>
class GenomeIndex
{
public:

	GenomeIndex() = default;
	GenomeIndex(const GenomeIndex&) = default;
	GenomeIndex(GenomeIndex&&) = default;
	GenomeIndex& operator=(const GenomeIndex&) = default;
	GenomeIndex& operator=(GenomeIndex&&) = default;
	~GenomeIndex() = default;

	GenomeIndex(const GenomicRange& genomic_range) {

		this->set(genomic_range);
	}

	GenomeIndex(const QStringList& peak_names) requires (Mode == IndexMode::Position) {

		this->set(peak_names);
	}

	void set_error_code(int n) {

		this->error_code_ = n;
	}

	bool success(int n) const {

		return this->error_code_ != n;
	}

	std::pair<int, int>	find_offset(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::Position);

	std::tuple<int, int, int, int>	find_offset(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::Strand);

	std::tuple<int, int, int, int, int, int, int, int>	find_offset(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::IdStrand);

	bool is_overlap(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::Position);

	QVector<int> find_overlap_ranges(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::ID);

	std::pair<int, int>	find_location(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::ID);

	// Note : we count insertion here, so both ends are included
	int	count_insertion(const QString& sequence_name, int start, int end) const
		requires (Mode == IndexMode::Position);

	void clear() {

		this->major_index_.clear();

		this->minor_index_.clear();
	}

	void set(const GenomicRange& genomic_range);

	void set(const QStringList& names) requires (Mode == IndexMode::Position);

	void reorganize();

	void append_index(const QString& sequence_name, int start, int end) requires (Mode == IndexMode::Position);

	void append_index(const QString& sequence_name, int start, int end, int id) requires (Mode == IndexMode::ID);

	void append_index(const QString& sequence_name, int start, int end, char strand) requires (Mode == IndexMode::Strand);

private:

	std::unordered_map<QString, ChromIndex> major_index_;

	std::unordered_map<QString, ChromIndex> minor_index_;

	int error_code_{ -1 };
};

template <IndexMode Mode>
void GenomeIndex<Mode>::append_index(const QString& sequence_name, int start, int end, int id)
	requires (Mode == IndexMode::ID)
{
	this->major_index_[sequence_name].append(start, end, id);
}

template <IndexMode Mode>
void GenomeIndex<Mode>::append_index(const QString& sequence_name, int start, int end, char strand)
	requires (Mode == IndexMode::Strand)
{
	if (strand == '-') {

		this->minor_index_[sequence_name].append(start, end);
	}
	else {

		this->major_index_[sequence_name].append(start, end);
	}
}

template <IndexMode Mode>
void GenomeIndex<Mode>::append_index(const QString& sequence_name, int start, int end)
	requires (Mode == IndexMode::Position)
{
	this->major_index_[sequence_name].append(start, end);
}

template <IndexMode Mode>
void GenomeIndex<Mode>::reorganize() {

	for (auto&& [seq_name, seq_index] : this->major_index_) {

		auto&& [start, end, id] = seq_index;

		auto order = custom::order(start);

		start = custom::reordered(start, order);
		end = custom::reordered(end, order);

		if (!id.empty()) {
			id = custom::reordered(id, order);
		}
	}

	for (auto&& [seq_name, seq_index] : this->minor_index_) {

		auto&& [start, end, id] = seq_index;

		auto order = custom::order(start);

		start = custom::reordered(start, order);
		end = custom::reordered(end, order);
		if (!id.empty()) {
			id = custom::reordered(id, order);
		}

	}
}

template <IndexMode Mode>
std::pair<int, int> GenomeIndex<Mode>::find_location(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::ID)
{
	auto ret = std::make_pair(this->error_code_, this->error_code_);

	auto iter = this->major_index_.find(sequence_name);
	if (iter == this->major_index_.cend()) {
		return ret;
	}

	auto&& [chr_start, chr_end, chr_id] = iter->second;

	if (chr_start.empty()) {
		return ret;
	}

	auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
	int loc = std::distance(chr_start.begin(), upper_bound_it);

	if (loc != 0) {
		int prev_loc = loc - 1;
		int prev_end = chr_end[prev_loc];
		int prev_id = chr_id[prev_loc];

		if (prev_end > start) {
			ret.first = prev_id;
		}
		if (prev_end > end) {
			ret.second = prev_id;
			return ret;
		}
	}

	int size = chr_start.size();
	for (int i = loc; i < size; ++i) {

		int current_start = chr_start[i];

		if (current_start > end) {
			return ret;
		}

		if (current_start <= end && chr_end[i] > end) {
			ret.second = chr_id[i];
			return ret;
		}
	}

	return ret;
}

template <IndexMode Mode>
std::tuple<int, int, int, int, int, int, int, int> GenomeIndex<Mode>::find_offset(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::IdStrand)
{
	auto ret = std::make_tuple(
		this->error_code_,
		this->error_code_,
		this->error_code_,
		this->error_code_,
		this->error_code_,
		this->error_code_,
		this->error_code_,
		this->error_code_
	);

	auto iter = this->major_index_.find(sequence_name);
	if (iter != this->major_index_.cend()) {

		auto&& [chr_start, chr_end, chr_id] = iter->second;

		if (chr_start.empty()) {
			goto next_strand_2;
		}

		auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
		int loc = std::distance(chr_start.begin(), upper_bound_it);

		if (loc != 0) {
			int prev_loc = loc - 1;
			int prev_start = chr_start[prev_loc];
			int prev_end = chr_end[prev_loc];
			int prev_id = chr_id[prev_loc];

			if (prev_end > start) {
				std::get<0>(ret) = prev_id;
				std::get<1>(ret) = start - prev_start;
			}
			if (prev_end > end) {
				std::get<2>(ret) = prev_id;
				std::get<3>(ret) = end - prev_start;
				goto next_strand_2;
			}
		}

		int size = chr_start.size();
		for (int i = loc; i < size; ++i) {

			int current_start = chr_start[i];

			if (current_start > end) {
				goto next_strand_2;
			}

			if (current_start <= end && chr_end[i] > end) {
				std::get<2>(ret) = chr_id[i];
				std::get<3>(ret) = end - current_start;
				goto next_strand_2;
			}
		}
	}

next_strand_2:
	auto iter2 = this->minor_index_.find(sequence_name);
	if (iter2 != this->minor_index_.cend()) {

		auto&& [chr_start, chr_end, chr_id] = iter2->second;

		if (chr_start.empty()) {
			return ret;
		}

		auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
		int loc = std::distance(chr_start.begin(), upper_bound_it);

		if (loc != 0) {
			int prev_loc = loc - 1;
			int prev_start = chr_start[prev_loc];
			int prev_end = chr_end[prev_loc];
			int prev_id = chr_id[prev_loc];

			if (prev_end > start) {
				std::get<4>(ret) = prev_id;
				std::get<5>(ret) = prev_end - start - 1;
			}
			if (prev_end > end) {
				std::get<6>(ret) = prev_id;
				std::get<7>(ret) = prev_end - end - 1;
				return ret;
			}
		}

		int size = chr_start.size();
		for (int i = loc; i < size; ++i) {

			int current_start = chr_start[i];

			if (current_start > end) {
				return ret;
			}

			if (current_start <= end && chr_end[i] > end) {
				std::get<6>(ret) = chr_id[i];
				std::get<7>(ret) = chr_end[i] - end - 1;
				return ret;
			}
		}

		return ret;
	}

	return ret;
}

template <IndexMode Mode>
std::tuple<int, int, int, int> GenomeIndex<Mode>::find_offset(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::Strand)
{
	auto ret = std::make_tuple(
		this->error_code_,
		this->error_code_,
		this->error_code_,
		this->error_code_
	);

	auto iter = this->major_index_.find(sequence_name);
	if (iter != this->major_index_.cend()) {

		auto&& [chr_start, chr_end, chr_id] = iter->second;

		if (chr_start.empty()) {
			goto next_strand;
		}

		auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
		int loc = std::distance(chr_start.begin(), upper_bound_it);

		if (loc != 0) {
			int prev_loc = loc - 1;
			int prev_start = chr_start[prev_loc];
			int prev_end = chr_end[prev_loc];

			if (prev_end > start) {
				std::get<0>(ret) = start - prev_start;
			}
			if (prev_end > end) {
				std::get<1>(ret) = end - prev_start;
				goto next_strand;
			}
		}

		int size = chr_start.size();
		for (int i = loc; i < size; ++i) {

			int current_start = chr_start[i];

			if (current_start > end) {
				goto next_strand;
			}

			if (current_start <= end && chr_end[i] > end) {
				std::get<1>(ret) = end - current_start;
				goto next_strand;
			}
		}
	}

next_strand:
	auto iter2 = this->minor_index_.find(sequence_name);
	if (iter2 != this->minor_index_.cend()) {

		auto&& [chr_start, chr_end, chr_id] = iter2->second;

		if (chr_start.empty()) {
			return ret;
		}

		auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
		int loc = std::distance(chr_start.begin(), upper_bound_it);

		if (loc != 0) {
			int prev_loc = loc - 1;
			int prev_start = chr_start[prev_loc];
			int prev_end = chr_end[prev_loc];

			if (prev_end > start) {
				std::get<2>(ret) = prev_end - start - 1;
			}
			if (prev_end > end) {
				std::get<3>(ret) = prev_end - end - 1;
				return ret;
			}
		}

		int size = chr_start.size();
		for (int i = loc; i < size; ++i) {

			int current_start = chr_start[i];

			if (current_start > end) {
				return ret;
			}

			if (current_start <= end && chr_end[i] > end) {
				std::get<3>(ret) = chr_end[i] - end - 1;
				return ret;
			}
		}

		return ret;
	}

	return ret;
}

/*
	Note : Although end is not contained in segments, we count insertion sites here. So no --end;
*/
template <IndexMode Mode>
std::pair<int, int> GenomeIndex<Mode>::find_offset(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::Position)
{
	auto ret = std::make_pair(this->error_code_, this->error_code_);

	auto iter = this->major_index_.find(sequence_name);
	if (iter != this->major_index_.cend()) {

		auto&& [chr_start, chr_end, chr_id] = iter->second;

		if (chr_start.empty()) {
			return ret;
		}

		auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
		int loc = std::distance(chr_start.begin(), upper_bound_it);

		if (loc != 0) {
			int prev_loc = loc - 1;
			int prev_start = chr_start[prev_loc];
			int prev_end = chr_end[prev_loc];

			if (prev_end > start) {
				std::get<0>(ret) = start - prev_start;
			}
			if (prev_end > end) {
				std::get<1>(ret) = end - prev_start;
				return ret;
			}
		}

		int size = chr_start.size();
		for (int i = loc; i < size; ++i) {

			int current_start = chr_start[i];

			if (current_start > end) {
				return ret;
			}

			if (current_start <= end && chr_end[i] > end) {
				std::get<1>(ret) = end - current_start;
				return ret;
			}
		}
	}

	return ret;
}

template <IndexMode Mode>
bool GenomeIndex<Mode>::is_overlap(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::Position)
{

	auto iter = this->major_index_.find(sequence_name);
	if (iter == this->major_index_.cend()) {
		return false;
	}

	auto&& [chr_start, chr_end, chr_id] = iter->second;

	if (chr_start.empty()) {
		return false;
	}

	auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
	int loc = std::distance(chr_start.begin(), upper_bound_it);

	if (loc != 0) {
		int prev_loc = loc - 1;
		int prev_end = chr_end[prev_loc];

		if (prev_end > start) {
			return true;
		}
	}

	if (loc < chr_start.size()) {
		if (chr_start[loc] < end) {
			return true;
		}
		else {
			return false;
		}
	}

	return false;
}

template <IndexMode Mode>
QVector<int> GenomeIndex<Mode>::find_overlap_ranges(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::ID)
{

	QVector<int> ret;

	auto iter = this->major_index_.find(sequence_name);
	if (iter == this->major_index_.cend()) {
		return ret;
	}

	auto&& [chr_start, chr_end, chr_id] = iter->second;

	if (chr_start.empty()) {
		return ret;
	}

	auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
	int loc = std::distance(chr_start.begin(), upper_bound_it);

	if (loc != 0) {
		int prev_loc = loc - 1;
		int prev_end = chr_end[prev_loc];

		if (prev_end > start) {
			ret << chr_id[prev_loc];
		}
	}

	int size = chr_start.size();
	for (int i = loc; i < size; ++i) {
		int current_start = chr_start[i];

		if (end <= current_start) {
			return ret;
		}

		ret << chr_id[i];
	}

	return ret;
}

template <IndexMode Mode>
void GenomeIndex<Mode>::set(const GenomicRange& genomic_range) {

	this->clear();

	int n_range = genomic_range.rows();

	if constexpr (Mode == IndexMode::Position) {
		for (int i = 0; i < n_range; ++i) {
			QString seq = genomic_range.sequence_names_[i];
			int start_position = genomic_range.ranges_.start_[i];
			int end_position = genomic_range.ranges_.width_[i] + start_position;
			this->major_index_[seq].append(start_position, end_position);
		}
	}

	if constexpr (Mode == IndexMode::ID) {
		for (int i = 0; i < n_range; ++i) {
			QString seq = genomic_range.sequence_names_[i];
			int start_position = genomic_range.ranges_.start_[i];
			int end_position = genomic_range.ranges_.width_[i] + start_position;
			this->major_index_[seq].append(start_position, end_position, i);
		}
	}

	if constexpr (Mode == IndexMode::Strand) {
		for (int i = 0; i < n_range; ++i) {
			QString seq = genomic_range.sequence_names_[i];
			int start_position = genomic_range.ranges_.start_[i];
			int end_position = genomic_range.ranges_.width_[i] + start_position;

			if (genomic_range.strand_[i] == '-') {
				this->minor_index_[seq].append(start_position, end_position);
			}
			else {
				this->major_index_[seq].append(start_position, end_position);
			}
		}
	}

	if constexpr (Mode == IndexMode::IdStrand) {
		for (int i = 0; i < n_range; ++i) {
			QString seq = genomic_range.sequence_names_[i];
			int start_position = genomic_range.ranges_.start_[i];
			int end_position = genomic_range.ranges_.width_[i] + start_position;

			if (genomic_range.strand_[i] == '-') {
				this->minor_index_[seq].append(start_position, end_position, i);
			}
			else {
				this->major_index_[seq].append(start_position, end_position, i);
			}
		}
	}
}

template <IndexMode Mode>
void GenomeIndex<Mode>::set(const QStringList& names)
	requires (Mode == IndexMode::Position)
{

	this->clear();

	for (const auto& peak_name : names) {
		auto [sequence_name, start, end, success] = custom::string_to_peak(peak_name);
		if (success) {
			this->major_index_[sequence_name].append(start, end);
		}
	}

	this->reorganize();
}

template <IndexMode Mode>
int	GenomeIndex<Mode>::count_insertion(const QString& sequence_name, int start, int end) const
	requires (Mode == IndexMode::Position)
{

	auto iter = this->major_index_.find(sequence_name);
	if (iter == this->major_index_.cend()) {
		return 0;
	}

	int ret{ 0 };

	auto&& [chr_start, chr_end, chr_id] = iter->second;

	if (chr_start.empty()) {
		return ret;
	}

	auto upper_bound_it = std::ranges::upper_bound(chr_start, start);
	int loc = std::distance(chr_start.begin(), upper_bound_it);

	if (loc != 0) {
		int prev_loc = loc - 1;
		int prev_end = chr_end[prev_loc];

		if (prev_end > start) {
			++ret;
		}

		if (prev_end > end) {
			++ret;
			return ret;
		}
	}

	int size = chr_start.size();
	for (int i = loc; i < size; ++i) {
		int current_start = chr_start[i];

		if (end < current_start) {
			return ret;
		}

		if (current_start <= end && end < chr_end[i]) {
			++ret;
			return ret;
		}
	}

	return ret;
}