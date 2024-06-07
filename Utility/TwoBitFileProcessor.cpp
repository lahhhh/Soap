#include "TwoBitFileProcessor.h"

#include <QDebug>

template<int N>
char get_base(uint8_t byte) {
	constexpr int offset =  8 - N * 2;
	uint8_t base = (byte >> offset) & 0x3;
	switch (base)
	{
	case 0:
		return 'T';
	case 1:
		return 'C';
	case 2:
		return 'A';
	case 3:
		return 'G';
	default:
		return 'N';
	}
}

char get_base(uint8_t byte, int N) {
	const int offset =  8 - N * 2;
	uint8_t base = (byte >> offset) & 0x3;
	switch (base)
	{
	case 0:
		return 'T';
	case 1:
		return 'C';
	case 2:
		return 'A';
	case 3:
		return 'G';
	default:
		return 'N';
	}
}

// 1 <= N <= 16
char get_base(uint32_t block, int N) {
	int byte_loc = (N - 1) / 4;
	int byte_offset = 8 - ((N - 1) % 4 + 1) * 2;
	uint8_t base = (block >> (byte_loc * 8 + byte_offset)) & 0x3;
	switch (base)
	{
	case 0:
		return 'T';
	case 1:
		return 'C';
	case 2:
		return 'A';
	case 3:
		return 'G';
	default:
		return 'N';
	}
}

TwoBitFileProcessor::TwoBitFileProcessor(const QString& file_name): 
	file_name_(file_name),
	file_(file_name)
{};

QStringList TwoBitFileProcessor::get_sequence(const GenomicRange& genomic_range) {

	if (!this->loaded_) {
		if (!this->load_sequence()) {
			return {};
		}
	}

	const qsizetype size = genomic_range.size();
	QStringList ret(size);
	for (qsizetype i = 0; i < size; ++i) {
		auto [sequence_name, start, end, strand] = genomic_range.at(i);
		ret[i] = this->get_sequence(sequence_name, start, end);
	}
	return ret;
};

std::vector<std::string> TwoBitFileProcessor::get_std_sequence(const GenomicRange& genomic_range) {

	if (!this->loaded_) {
		if (!this->load_sequence()) {
			return {};
		}
	}

	const qsizetype size = genomic_range.size();
	std::vector<std::string> ret(size);
	for (qsizetype i = 0; i < size; ++i) {
		auto [sequence_name, start, end, strand] = genomic_range.at(i);
		ret[i] = this->get_std_sequence(sequence_name, start, end);
	}
	return ret;
};

std::string TwoBitFileProcessor::get_std_whole_sequence(const QString& sequence_name) {

	if (!this->loaded_) {
		if (!this->load_sequence()) {
			return {};
		}
	}

	if (!this->sequence_dna_size_.contains(sequence_name)) {
		return {};
	}
	uint32_t start = 1;
	uint32_t end = this->sequence_dna_size_[sequence_name];
	return this->get_std_sequence(sequence_name, start, end + 1);
};

std::string TwoBitFileProcessor::get_std_sequence(const QString& sequence_name, uint32_t start, uint32_t end) {

	if (!this->loaded_) {
		if (!this->load_sequence()) {
			return {};
		}
	}

	--end;

	if (!this->sequence_dna_size_.contains(sequence_name)) {
		return {};
	}
	uint32_t dna_size = this->sequence_dna_size_[sequence_name];
	if (start < 1 || end > dna_size) {
		return {};
	}
	if (start > end) {
		return {};
	}

	uint32_t region_size = end - start + 1;
	uint32_t offset = this->sequence_offset2_[sequence_name];
	uint32_t start_block = (start - 1) / 16;
	uint32_t local_offset = offset + start_block * 4;
	uint32_t end_block = (end - 1) / 16;

	this->file_.seek(local_offset);

	uint32_t first_base_offset = start % 16;
	if (first_base_offset == 0) {
		first_base_offset = 16;
	}
	uint32_t last_base_offset = end % 16;
	if (last_base_offset == 0) {
		last_base_offset = 16;
	}

	uint32_t blocks_to_read = end_block - start_block + 1;

	QVector<uint32_t> four_byte_dna(blocks_to_read);

	QDataStream ds(&this->file_);
	for (uint32_t i = 0; i < blocks_to_read; ++i) {
		ds >> four_byte_dna[i];
	}

	if (this->byte_swapped_) {
		for (uint32_t i = 0; i < blocks_to_read; ++i) {
			four_byte_dna[i] = _byteswap_ulong(four_byte_dna[i]);
		}
	}

	std::string ret;
	ret.resize(region_size);

	// first_base_offset : last_base_offset closed interval
	if (blocks_to_read == 1) {
		for (uint32_t i = 0; i < region_size; ++i) {
			ret[i] = get_base(four_byte_dna[0], i + first_base_offset);
		}

	}
	else {
		// first block
		int first_block_base_num = 17 - first_base_offset;
		for (uint32_t i = 0; i < first_block_base_num; ++i) {
			ret[i] = get_base(four_byte_dna[0], i + first_base_offset);
		}

		uint8_t* bytes = (uint8_t*)(four_byte_dna.data() + 1);
		for (uint32_t i = 0; i < 4 * (blocks_to_read - 2); ++i) {
			uint8_t byte = bytes[i];
			ret[i * 4 + first_block_base_num] = get_base<1>(byte);
			ret[i * 4 + first_block_base_num + 1] = get_base<2>(byte);
			ret[i * 4 + first_block_base_num + 2] = get_base<3>(byte);
			ret[i * 4 + first_block_base_num + 3] = get_base<4>(byte);
		}

		//last block
		for (uint32_t i = 1; i < last_base_offset + 1; ++i) {
			ret[i + 16 * (blocks_to_read - 2) + first_block_base_num - 1] = get_base(four_byte_dna[blocks_to_read - 1], i);
		}
	}

	// check N block
	const auto& n_starts = this->sequence_n_block_starts_[sequence_name];
	const auto& n_sizes = this->sequence_n_block_sizes_[sequence_name];
	qsizetype size = n_starts.size();
	for (qsizetype i = 0; i < size; ++i) {
		uint32_t n_start = n_starts[i], n_size = n_sizes[i], n_end = n_start + n_size;
		if (n_end < start) {
			continue;
		}
		if (n_start > end) {
			break;
		}
		if (n_start < start) {
			n_start = start;
		}
		if (n_end > end) {
			n_end = end;
		}
		for (qsizetype i = n_start - start; i < n_end - start; ++i) {
			ret[i] = 'N';
		}
	}

	// check mask block
	const auto& mask_starts = this->sequence_mask_block_starts_[sequence_name];
	const auto& mask_sizes = this->sequence_mask_block_sizes_[sequence_name];
	size = mask_starts.size();
	for (qsizetype i = 0; i < size; ++i) {
		uint32_t mask_start = mask_starts[i], mask_size = mask_sizes[i], mask_end = mask_start + mask_size;
		if (mask_end < start) {
			continue;
		}
		if (mask_start > end) {
			break;
		}
		if (mask_start < start) {
			mask_start = start;
		}
		if (mask_end > end) {
			mask_end = end;
		}
		for (qsizetype i = mask_start - start; i < mask_end - start; ++i) {
			ret[i] = std::tolower(ret[i]);
		}
	}
	return ret;
};

QString TwoBitFileProcessor::get_sequence(const QString& sequence_name, uint32_t start, uint32_t end) {

	if (!this->loaded_) {
		if (!this->load_sequence()) {
			return {};
		}
	}

	--end;

	if (!this->sequence_dna_size_.contains(sequence_name)) {
		return {};
	}
	uint32_t dna_size = this->sequence_dna_size_[sequence_name];
	if (start < 1 || end > dna_size) {
		return {};
	}
	if (start > end) {
		return {};
	}
	uint32_t region_size = end - start + 1;
	uint32_t offset = this->sequence_offset2_[sequence_name];
	uint32_t start_block = (start - 1) / 16;
	uint32_t local_offset = offset + start_block * 4;
	uint32_t end_block = (end - 1) / 16;

	this->file_.seek(local_offset);

	uint32_t first_base_offset = start % 16;
	if (first_base_offset == 0) {
		first_base_offset = 16;
	}
	uint32_t last_base_offset = end % 16;
	if (last_base_offset == 0) {
		last_base_offset = 16;
	}

	uint32_t blocks_to_read = end_block - start_block + 1;

	QVector<uint32_t> four_byte_dna(blocks_to_read);

	QDataStream ds(&this->file_);
	for (uint32_t i = 0; i < blocks_to_read; ++i) {
		ds >> four_byte_dna[i];
	}

	if (this->byte_swapped_) {
		for (uint32_t i = 0; i < blocks_to_read; ++i) {
			four_byte_dna[i] = _byteswap_ulong(four_byte_dna[i]);
		}
	}


	QString ret;
	ret.resize(region_size);

	// first_base_offset : last_base_offset closed interval
	if (blocks_to_read == 1) {
		for (uint32_t i = 0; i < region_size; ++i) {
			ret[i] = get_base(four_byte_dna[0], i + first_base_offset);
		}

	}
	else {
		// first block
		int first_block_base_num = 17 - first_base_offset;
		for (uint32_t i = 0; i < first_block_base_num; ++i) {
			ret[i] = get_base(four_byte_dna[0], i + first_base_offset);
		}

		uint8_t* bytes = (uint8_t*)(four_byte_dna.data() + 1);
		for (uint32_t i = 0; i < 4 * (blocks_to_read - 2); ++i) {
			uint8_t byte = bytes[i];
			ret[i * 4 + first_block_base_num] = get_base<1>(byte);
			ret[i * 4 + first_block_base_num + 1] = get_base<2>(byte);
			ret[i * 4 + first_block_base_num + 2] = get_base<3>(byte);
			ret[i * 4 + first_block_base_num + 3] = get_base<4>(byte);
		}

		//last block
		for (uint32_t i = 1; i < last_base_offset + 1; ++i) {
			ret[i + 16 * (blocks_to_read - 2) + first_block_base_num - 1] = get_base(four_byte_dna[blocks_to_read - 1], i);
		}
	}

	// check N block
	const auto& n_starts = this->sequence_n_block_starts_[sequence_name];
	const auto& n_sizes = this->sequence_n_block_sizes_[sequence_name];
	qsizetype size = n_starts.size();
	for (qsizetype i = 0; i < size; ++i) {
		uint32_t n_start = n_starts[i], n_size = n_sizes[i], n_end = n_start + n_size;
		if (n_end < start) {
			continue;
		}
		if (n_start > end) {
			break;
		}
		if (n_start < start) {
			n_start = start;
		}
		if (n_end > end) {
			n_end = end;
		}
		for (qsizetype i = n_start - start; i < n_end - start; ++i) {
			ret[i] = 'N';
		}
	}

	// check mask block
	const auto& mask_starts = this->sequence_mask_block_starts_[sequence_name];
	const auto& mask_sizes = this->sequence_mask_block_sizes_[sequence_name];
	size = mask_starts.size();
	for (qsizetype i = 0; i < size; ++i) {
		uint32_t mask_start = mask_starts[i], mask_size = mask_sizes[i], mask_end = mask_start + mask_size;
		if (mask_end < start) {
			continue;
		}
		if (mask_start > end) {
			break;
		}
		if (mask_start < start) {
			mask_start = start;
		}
		if (mask_end > end) {
			mask_end = end;
		}
		for (qsizetype i = mask_start - start; i < mask_end - start; ++i) {
			ret[i] = ret[i].toLower();
		}
	}

	return ret;
};

void TwoBitFileProcessor::release() {
	if (this->file_.isOpen()) {
		this->file_.close();
	}
	this->sequence_count_ = 0;
	this->file_size_ = 0;
	this->sequence_offset_.clear();
	this->sequence_offset2_.clear();
	this->sequence_dna_size_.clear();
	this->sequence_packed_size_.clear();
	this->sequence_bytes_.clear();
	this->sequence_n_block_starts_.clear();
	this->sequence_n_block_sizes_.clear();
	this->sequence_mask_block_starts_.clear();
	this->sequence_mask_block_sizes_.clear();
	this->loaded_ = false;
};

void TwoBitFileProcessor::set_sequence_file(const QString& sequence_file_name) {

	if (this->loaded_) {
		this->release();
	}
	
	this->file_.setFileName(sequence_file_name);
	this->file_name_ = sequence_file_name;
}

bool TwoBitFileProcessor::load_sequence() {

	if (!this->file_.open(QIODevice::ReadOnly)) {
		return false;
	}
	
	this->file_size_ = this->file_.size();
	QDataStream ds(&this->file_);
	int signature, version, sequence_count, reserved;
	ds >> signature >> version >> sequence_count >> reserved;
	if (signature != 0x1A412743) {
		this->byte_swapped_ = true;
		this->sequence_count_ = _byteswap_ulong(sequence_count);
	}
	else {
		this->byte_swapped_ = false;
		this->sequence_count_ = sequence_count;
	}
	this->load_index(ds);
	this->initialize_sequence();
	this->loaded_ = true;

	return true;
}

void TwoBitFileProcessor::load_index(QDataStream& ds) {
	qsizetype remaining = this->sequence_count_;
	while (remaining > 0) {
		uint8_t name_size;
		ds >> name_size;
		QVector<char> name(name_size + 1);
		for (uint8_t i = 0; i < name_size; ++i) {
			ds >> name[i];
		}
		name[name_size] = '\0';
		QString sequence_name = _Cs standardize_chromosome_name(QString::fromUtf8(name.data()));
		uint32_t offset;
		ds >> offset;
		if (this->byte_swapped_) {
			offset = _byteswap_ulong(offset);
		}
		this->sequence_offset_[sequence_name] = offset;
		--remaining;
	}
};

void TwoBitFileProcessor::initialize_sequence() {
	for (const auto& sequence_name : this->sequence_offset_.keys()) {
		
		uint32_t offset = this->sequence_offset_[sequence_name];
		this->file_.seek(offset);
		QDataStream ds(&this->file_);
		uint32_t dna_size, n_block_count;
		ds >> dna_size >> n_block_count;
		if (this->byte_swapped_) {
			dna_size = _byteswap_ulong(dna_size);
			n_block_count = _byteswap_ulong(n_block_count);
		}
		this->sequence_dna_size_[sequence_name] = dna_size;
		this->sequence_bytes_[sequence_name] = (dna_size + 3) / 4;
		this->sequence_packed_size_[sequence_name] = (dna_size + 15) / 16;


		auto& starts = this->sequence_n_block_starts_[sequence_name];
		starts.resize(n_block_count);
		for (uint32_t i = 0; i < n_block_count; ++i) {
			uint32_t start;
			ds >> start;
			if (this->byte_swapped_) {
				start = _byteswap_ulong(start);
			}
			starts[i] = start;
		}


		auto& sizes = this->sequence_n_block_sizes_[sequence_name];
		sizes.resize(n_block_count);
		for (uint32_t i = 0; i < n_block_count; ++i) {
			uint32_t size;
			ds >> size;
			if (this->byte_swapped_) {
				size = _byteswap_ulong(size);
			}
			sizes[i] = size;
		}

		uint32_t mask_raw_count;
		ds >> mask_raw_count;
		if (this->byte_swapped_) {
			mask_raw_count = _byteswap_ulong(mask_raw_count);
		}

		auto& mask_starts = this->sequence_mask_block_starts_[sequence_name];
		mask_starts.resize(mask_raw_count);
		for (uint32_t i = 0; i < mask_raw_count; ++i) {
			uint32_t start;
			ds >> start;
			if (this->byte_swapped_) {
				start = _byteswap_ulong(start);
			}
			mask_starts[i] = start;
		}

		auto& mask_sizes = this->sequence_mask_block_sizes_[sequence_name];
		mask_sizes.resize(mask_raw_count);
		for (uint32_t i = 0; i < mask_raw_count; ++i) {
			uint32_t size;
			ds >> size;
			if (this->byte_swapped_) {
				size = _byteswap_ulong(size);
			}
			mask_sizes[i] = size;
		}

		this->sequence_offset2_[sequence_name] = this->file_.pos() + 4;
	}
};

