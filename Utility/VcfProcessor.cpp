#include "VcfProcessor.h"

VcfProcessor::VcfProcessor(const std::string& vcf_file_name, int cache_size) :
	file_name_(vcf_file_name), cache_size_(cache_size)
{
	if (this->cache_size_ <= 0) {
		this->cache_size_ = 10000;
	}

	this->cache_.resize(this->cache_size_);
};

bool VcfProcessor::check_header(const StringSplitView& view) {

	if (view.size() < 10) {
		return false;
	}

	if (view[0] != "#CHROM") {
		return false;
	}
	if (view[1] != "POS") {
		return false;
	}
	if (view[2] != "ID") {
		return false;
	}
	if (view[3] != "REF") {
		return false;
	}
	if (view[4] != "ALT") {
		return false;
	}
	if (view[5] != "QUAL") {
		return false;
	}
	if (view[6] != "FILTER") {
		return false;
	}
	if (view[7] != "INFO") {
		return false;
	}
	if (view[8] != "FORMAT") {
		return false;
	}

	return true;
};