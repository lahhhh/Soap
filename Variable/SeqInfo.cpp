#include "SeqInfo.h"

void SeqInfo::append(const SeqInfo& seq_info) {
	this->sequence_names_ << seq_info.sequence_names_;
	this->sequence_lengths_ << seq_info.sequence_lengths_;
	this->is_circular_ << seq_info.is_circular_;
	this->genome_ << seq_info.genome_;
};