#include "VcfProcessor.h"

VcfProcessor::VcfProcessor(const std::string& bam_file_name, int cache_size) :
	file_name_(bam_file_name), cache_size_(cache_size)
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

bool VcfProcessor::generate_snp_matrix(const std::string& output_name, int min_mutation_depth) {

	this->in_.open(this->file_name_, std::ios::in);

	if (!this->in_.is_open()) {
		return false;
	}

	std::string line;
	while (std::getline(this->in_, line)) {
		if (!line.starts_with("##")) {
			break;
		}
	}

	if (this->in_.eof()) {
		return false;
	}

	auto sp = custom::split_view(line, '\t');

	if (!this->check_header(sp)) {
		return false;
	}

	std::string header{ "SNP" };

	int n_column = sp.size();
	int n_sample = n_column - 9;

	for (int i = 0; i < n_sample; ++i) {
		header += '\t';
		header += sp[i + 9];
	}
	header += '\n';

	this->out_.open(output_name, std::ios::out);

	if (!this->out_.is_open()) {
		return false;
	}

	this->out_ << header;

	int cache_storage{ 0 };
	bool error{ false };

	while (true) {

		while (cache_storage < this->cache_size_
			&& std::getline(this->in_, this->cache_[cache_storage])) {

			++cache_storage;
		}

	#pragma omp parallel for
		for (int i = 0; i < cache_storage; ++i) {
			auto sv = custom::split_view(this->cache_[i], '\t');

			if (sv.size() != n_column) {
				error = true;
				break;
			}

			std::string snp;
			snp.append(sv[0].data(), sv[0].size());
			snp += '-';
			snp.append(sv[1].data(), sv[1].size());
			snp += '-';
			snp.append(sv[3].data(), sv[3].size());
			snp += '-';
			snp.append(sv[4].data(), sv[4].size());

			std::string content(n_sample * 2 + 1, '\t');
			for (int j = 9; j < n_column; ++j) {
				if (sv[j].starts_with("0/0") || 
					sv[j].starts_with("0|0") ) {
					content[2 * (j - 9) + 1] = '0';
				}
				else if (
					sv[j].starts_with("./.") ||
					sv[j].starts_with(".|.")) {
					content[2 * (j - 9) + 1] = '2';
				}
				else{

					auto sample_info = custom::split(sv[j], ':');

					if (sample_info.size() < 2) {
						error = true;
						break;
					}
					else {

						auto&& genotype = sample_info[0];

						if (genotype.size() < 3) {
							error = true;
							break;
						}
					
						auto allelic_depth = custom::split(sample_info[1], ',');
						int n_alle = allelic_depth.size();
						if (n_alle < 2) {
							error = true;
							break;
						}

						bool is_mut{ false };
						for (int k = 1; k < n_alle; ++k) {
							int depth = std::stoi(allelic_depth[k]);
							if (depth >= min_mutation_depth) {

								content[2 * (j - 9) + 1] = '1';
								is_mut = true;
								break;
							}
						}
						if (!is_mut) {
							content[2 * (j - 9) + 1] = '0';
						}
					}

				}				
			}

			content[n_sample * 2] = '\n';

			this->cache_[i] = snp + content;
		}

		if (error) {
			break;
		}

		for (int i = 0; i < cache_storage; ++i) {
			this->out_ << this->cache_[i];
		}

		if (cache_storage < this->cache_size_) {
			break;
		}

		cache_storage = 0;
	}

	return !error;
};