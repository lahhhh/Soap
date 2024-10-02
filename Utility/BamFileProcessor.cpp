#include "BamFileProcessor.h"

#define STORAGE_CAPACITY 800000
#define BUFFER_CAPACITY 65536

#define NEED_BYTES(n_bytes) \
	if(!this->need_bytes(n_bytes)){ return false; }

#define READ_VAR(var_name) \
	this->read(&var_name, sizeof(var_name));

BamFileProcessor::BamFileProcessor(const std::string& bam_file_name) :
	bam_file_name_(bam_file_name)
{
	this->uncompressed_ = new char[STORAGE_CAPACITY];

	this->buffer_ = new char[BUFFER_CAPACITY];

	this->stream_.zalloc = Z_NULL;
	this->stream_.zfree = Z_NULL;
	this->stream_.opaque = Z_NULL;

};

BamFileProcessor::~BamFileProcessor() {
	delete[] this->uncompressed_;
	delete[] this->buffer_;
};

bool BamFileProcessor::read_block() {

	auto buffer_left_size = this->uncompressed_end_ - this->reading_location_;

	// no enough space 
	if (STORAGE_CAPACITY - buffer_left_size < 65536) {
		return false;
	}

	memmove(this->uncompressed_, this->uncompressed_ + this->reading_location_, buffer_left_size);

	this->reading_location_ = 0;

	this->uncompressed_end_ = buffer_left_size;

	long long start_pos = this->in_.tellg();

	long long left_size = this->file_size_ - start_pos;

	if (left_size < 26) {
		return false;
	}

	this->in_.seekg(16, std::ios::cur);

	uint16_t block_size{ 0 };

	this->in_.read((char*)&block_size, 2);

	if (left_size < block_size) {
		return false;
	}

	this->in_.seekg(block_size - 21, std::ios::cur);

	uint32_t input_size{ 0 };

	this->in_.read((char*)&input_size, 4);

	if (input_size == 0) {
		if (this->in_.eof()) {
			return false;
		}
		else {
			return this->read_block();
		}
	}

	this->in_.seekg(- block_size - 1, std::ios::cur);

	this->in_.read(this->buffer_, block_size + 1);

	this->stream_.next_in = (z_const Bytef*)this->buffer_;
	this->stream_.avail_in = block_size + 1;

	this->stream_.next_out = (z_const Bytef*)(this->uncompressed_ + buffer_left_size);
	this->stream_.avail_out = input_size;

	inflateInit2(&this->stream_, 16 + MAX_WBITS);

	inflate(&this->stream_, Z_NO_FLUSH);

	if (this->stream_.avail_in == 0) {
		
		this->uncompressed_end_ += input_size;

		inflateEnd(&this->stream_);

		return true;
	}
	else {
		inflateEnd(&this->stream_);

		return false; // file reading will terminate.
	}
};

bool BamFileProcessor::open_file() {

	this->in_.open(this->bam_file_name_, std::ios::in | std::ios::binary);

	if (this->in_.is_open()) {

		this->in_.seekg(0, std::ios::end);

		this->file_size_ = this->in_.tellg();

		this->in_.seekg(0, std::ios::beg);

		return true;
	}
	else {
		return false;
	}
};

void BamFileProcessor::advance(std::size_t length) {

	this->reading_location_ += length;
}

std::string BamFileProcessor::ref_name() {

	int32_t ref_id = this->get_ref_id();

	if (ref_id >= 0 && ref_id < this->ref_names_.size()) {
		return this->ref_names_[ref_id];
	}
	else {
		return std::string("*");
	}
}

uint16_t BamFileProcessor::get_n_cigar_op() {

	uint16_t n_cigar_op{ 0 };

	std::memcpy(&n_cigar_op, this->alignment_ + 12, 2);

	return n_cigar_op;
}

char* BamFileProcessor::now() {

	return this->uncompressed_ + this->reading_location_;
}

std::string BamFileProcessor::get_base_quality() {

	auto seq_length = this->get_seq_length();

	char* base_quality_start = this->alignment_ + 32 + this->get_read_name_length() + 4 * this->get_n_cigar_op() + (seq_length + 1) / 2;

	std::string qual(base_quality_start, seq_length);

	for (uint32_t i = 0; i < seq_length; ++i) {
		qual[i] += 33;
	}

	return qual;
}

std::string BamFileProcessor::get_tag_string() {
	char* tag_start = this->alignment_ + 32 + this->get_read_name_length() + 4 * this->get_n_cigar_op() + (3 * this->get_seq_length() + 1) / 2;

	char* end = this->alignment_ + this->alignment_size_;

	std::string tag_string;

	char tag_name[3] = "\0\0";

	char value_type{ '\0' };

	while (tag_start < end) {

		if (value_type != '\0') {
			tag_string += '\t';
		}

		std::memcpy(tag_name, tag_start, 2);

		tag_start += 2;

		tag_string += tag_name;
		tag_string += ':';

		value_type = *tag_start;

		tag_string += value_type;
		tag_string += ':';


		++tag_start;

		switch (value_type)
		{
		case 'A':
			tag_string += *tag_start;
			++tag_start;
			break;

		case 'c':
			tag_string += std::to_string(*(int8_t*)tag_start);
			++tag_start;
			break;

		case 'C':
			tag_string += std::to_string(*(uint8_t*)tag_start);
			++tag_start;
			break;

		case 's':
			tag_string += std::to_string(*(int16_t*)tag_start);
			tag_start += 2;
			break;

		case 'S':
			tag_string += std::to_string(*(uint16_t*)tag_start);
			tag_start += 2;
			break;

		case 'i':
			tag_string += std::to_string(*(int32_t*)tag_start);
			tag_start += 4;
			break;

		case 'I':
			tag_string += std::to_string(*(uint32_t*)tag_start);
			tag_start += 4;
			break;

		case 'f':
			tag_string += std::to_string(*(float*)tag_start);
			tag_start += 4;
			break;

		case 'Z':
		{
			auto string_size_before = tag_string.size();
			tag_string += tag_start;
			auto string_size_after = tag_string.size();
			tag_start += (string_size_after - string_size_before + 1);
			break;
		}
		case 'H':
		{
			auto string_size_before = tag_string.size();
			tag_string += tag_start;
			auto string_size_after = tag_string.size();
			tag_start += (string_size_after - string_size_before + 1);
			break;
		}
		case 'B':
		{
			char sub_type = *tag_start;

			++tag_start;

			uint32_t count = *(uint32_t*)tag_start;

			tag_start += 4;

			switch (sub_type)
			{
			case 'c':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(int8_t*)tag_start);
					++tag_start;
				}
				break;

			case 'C':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(uint8_t*)tag_start);
					++tag_start;
				}
				break;

			case 's':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(int16_t*)tag_start);
					tag_start += 2;
				}
				break;

			case 'S':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(uint16_t*)tag_start);
					tag_start += 2;
				}
				break;

			case 'i':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(int32_t*)tag_start);
					tag_start += 4;
				}
				break;

			case 'I':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(uint32_t*)tag_start);
					tag_start += 4;
				}
				break;

			case 'f':
				for (uint32_t i = 0; i < count; ++i) {
					tag_string += std::to_string(*(float*)tag_start);
					tag_start += 4;
				}
				break;

			default:
				break;
			}

			break;
		}
		default:
			break;
		}

	}
	return tag_string;
}

uint8_t BamFileProcessor::get_mapq() {

	uint8_t mapq{ 0 };

	std::memcpy(&mapq, this->alignment_ + 9, 1);

	return mapq;
}

std::string BamFileProcessor::get_cigar_string() {
	const char* cigar_start = this->alignment_ + 32 + this->get_read_name_length();

	uint16_t n_cigar_op = this->get_n_cigar_op();

	std::string cigar_string;

	for (uint16_t i = 0; i < n_cigar_op; ++i) {

		uint32_t op = *(uint32_t*)(cigar_start + 4 * i);

		uint32_t op_length = op >> 4;

		cigar_string += std::to_string(op_length);

		uint32_t op_type = op & 0x15;

		char type_char{ '\0' };

		switch (op_type)
		{
		case 0:
			type_char = 'M';
			break;
		case 1:
			type_char = 'I';
			break;
		case 2:
			type_char = 'D';
			break;
		case 3:
			type_char = 'N';
			break;
		case 4:
			type_char = 'S';
			break;
		case 5:
			type_char = 'H';
			break;
		case 6:
			type_char = 'P';
			break;
		case 7:
			type_char = '=';
			break;
		case 8:
			type_char = 'X';
			break;
		default:
			break;
		}

		cigar_string += type_char;
	}

	return cigar_string;
}

std::string BamFileProcessor::get_seq() {

	const char* seq_start = this->alignment_ + 32 + this->get_read_name_length() + 4 * this->get_n_cigar_op();

	uint32_t seq_length = this->get_seq_length();

	std::string seq(seq_length, '\0');

	for (uint32_t i = 0; i < seq_length; ++i) {

		uint8_t type{ 0 };

		if (i % 2 == 0) {

			type = *((uint8_t*)(seq_start + i / 2));

			type >>= 4;
		}
		else {
			type = *((uint8_t*)(seq_start + i / 2));


			type &= 0xf;
		}

		switch (type)
		{
		case 0:
			seq[i] = '=';
			break;
		case 1:
			seq[i] = 'A';
			break;
		case 2:
			seq[i] = 'C';
			break;
		case 3:
			seq[i] = 'M';
			break;
		case 4:
			seq[i] = 'G';
			break;
		case 5:
			seq[i] = 'R';
			break;
		case 6:
			seq[i] = 'S';
			break;
		case 7:
			seq[i] = 'V';
			break;
		case 8:
			seq[i] = 'T';
			break;
		case 9:
			seq[i] = 'W';
			break;
		case 10:
			seq[i] = 'Y';
			break;
		case 11:
			seq[i] = 'H';
			break;
		case 12:
			seq[i] = 'K';
			break;
		case 13:
			seq[i] = 'D';
			break;
		case 14:
			seq[i] = 'B';
			break;
		case 15:
			seq[i] = 'N';
			break;
		default:
			break;
		}
	}

	return seq;
}

std::string BamFileProcessor::get_read_name() {

	return std::string(this->alignment_ + 32);
}

int32_t BamFileProcessor::get_template_length() {

	int32_t template_length{ 0 };

	std::memcpy(&template_length, this->alignment_ + 28, 4);

	return template_length;
}

int32_t BamFileProcessor::get_next_POS() {

	return this->get_next_pos() + 1;
}

int32_t BamFileProcessor::get_next_pos() {
	int32_t next_pos{ 0 };

	std::memcpy(&next_pos, this->alignment_ + 24, 4);

	return next_pos;
}

std::string BamFileProcessor::get_next_ref() {

	int32_t next_ref_id = this->get_ref_id();

	if (next_ref_id >= 0 && next_ref_id < this->ref_names_.size()) {
		return this->ref_names_[next_ref_id];
	}
	else {
		return std::string("*");
	}
}

int32_t BamFileProcessor::get_next_ref_id() {
	int32_t next_ref_id{ 0 };

	std::memcpy(&next_ref_id, this->alignment_, 4);

	return next_ref_id;
}

uint32_t BamFileProcessor::get_seq_length() {

	uint32_t seq_length{ 0 };

	std::memcpy(&seq_length, this->alignment_ + 16, 4);

	return seq_length;
}

bool BamFileProcessor::get_strand_direction() {

	uint16_t flag{ 0 };

	std::memcpy(&flag, this->alignment_ + 14, 2);

	return !(flag & 0x10);
}

bool BamFileProcessor::is_unmapped() {

	uint16_t flag{ 0 };

	std::memcpy(&flag, this->alignment_ + 14, 2);

	return flag & 0x4;
}

uint16_t BamFileProcessor::get_flag() {

	uint16_t flag{ 0 };

	std::memcpy(&flag, this->alignment_ + 14, 2);

	return flag;
}

uint8_t BamFileProcessor::get_read_name_length() {

	uint8_t read_name_length{ 0 };

	std::memcpy(&read_name_length, this->alignment_ + 8, 1);

	return read_name_length;
}

int32_t BamFileProcessor::get_POS() {

	return this->get_pos() + 1;
}

int32_t BamFileProcessor::get_pos() {

	int32_t pos{ 0 };

	std::memcpy(&pos, this->alignment_ + 4, 4);

	return pos;
}

int32_t BamFileProcessor::get_ref_id() {

	int32_t ref_id{ 0 };

	std::memcpy(&ref_id, this->alignment_, 4);

	return ref_id;
}

bool BamFileProcessor::need_bytes(std::size_t n_bytes)
{
	while (this->uncompressed_end_ - this->reading_location_ < n_bytes) {
		if (!this->read_block()) {
			return false;
		}
	}

	return true;
}

void BamFileProcessor::read(void* dest, std::size_t length) {
	std::memcpy(dest, this->uncompressed_ + this->reading_location_, length);

	this->reading_location_ += length;
}

bool BamFileProcessor::process_header() {

	if (!this->open_file()) {
		return false;
	}

	if (!this->skip_header()) {
		return false;
	}

	if (!this->record_reference()) {
		return false;
	}

	return true;
}

bool BamFileProcessor::skip_header() {

	NEED_BYTES(8);

	this->advance(4);

	uint32_t text_length{ 0 };

	READ_VAR(text_length);

	NEED_BYTES(text_length);

	this->advance(text_length);

	return true;

};

bool BamFileProcessor::record_reference() {

	uint32_t n_ref{ 0 };

	NEED_BYTES(4);

	READ_VAR(n_ref);

	this->ref_names_.reserve(n_ref);
	this->ref_length_.reserve(n_ref);

	uint32_t name_length{ 0 }, ref_length{ 0 };

	for (uint32_t i = 0; i < n_ref; ++i) {
		NEED_BYTES(4);

		READ_VAR(name_length);
				
		NEED_BYTES(4 + name_length);

		this->ref_names_.emplace_back(this->now(), name_length - 1);// null terminated

		this->advance(name_length);

		READ_VAR(ref_length);

		this->ref_length_.emplace_back(ref_length);
	}

	return true;

};

std::string BamFileProcessor::get_read() {

	/******* read name *******/

	std::string read = std::string(this->alignment_ + 32);

	read += '\t';

	/******* flag *******/

	uint16_t flag = this->get_flag();

	read += std::to_string(flag);

	read += '\t';

	/******* reference name *******/

	read += this->ref_name();

	read += '\t';

	/******* pos *******/

	int32_t pos = this->get_POS();

	read += std::to_string(pos);

	read += '\t';

	uint8_t read_name_length{ 0 };

	std::memcpy(&read_name_length, this->alignment_ + 8, 1);

	/******* mapq *******/

	uint8_t mapq = this->get_mapq();

	read += std::to_string(mapq);

	read += '\t';

	/******* cigar *******/

	char* cigar_start = this->alignment_ + 32 + read_name_length;

	uint16_t n_cigar_op = this->get_n_cigar_op();

	for (uint16_t i = 0; i < n_cigar_op; ++i) {

		uint32_t op = *(uint32_t*)(cigar_start + 4 * i);

		uint32_t op_length = op >> 4;

		read += std::to_string(op_length);

		uint32_t op_type = op & 0xf;

		char type_char{ '\0' };

		switch (op_type)
		{
		case 0:
			type_char = 'M';
			break;
		case 1:
			type_char = 'I';
			break;
		case 2:
			type_char = 'D';
			break;
		case 3:
			type_char = 'N';
			break;
		case 4:
			type_char = 'S';
			break;
		case 5:
			type_char = 'H';
			break;
		case 6:
			type_char = 'P';
			break;
		case 7:
			type_char = '=';
			break;
		case 8:
			type_char = 'X';
			break;
		default:
			break;
		}

		read += type_char;
	}

	read += '\t';

	/******* rnext *******/

	read += this->get_next_ref();

	read += '\t';

	/******* pnext *******/

	int32_t pnext = this->get_next_POS();

	read += std::to_string(pnext);

	read += '\t';

	/******* tlen *******/

	int32_t tlen = this->get_template_length();

	read += std::to_string(tlen);

	read += '\t';
	/******* seq *******/	

	read += this->get_seq();

	read += '\t';

	/******* base quality *******/

	read += this->get_base_quality();

	/******* tag *******/

	char* tag_start = this->alignment_ + 32 + read_name_length + 4 * n_cigar_op + (3 * this->get_seq_length() + 1) / 2;

	char* end = this->alignment_ + this->alignment_size_;

	char tag_name[3] = "\0\0";

	char value_type{ '\0' };

	while (tag_start < end) {

		read += '\t';

		std::memcpy(tag_name, tag_start, 2);

		tag_start += 2;

		read += tag_name;
		read += ':';

		value_type = *tag_start;

		read += value_type;
		read += ':';


		++tag_start;

		switch (value_type)
		{
		case 'A':
			read += *tag_start;
			++tag_start;
			break;

		case 'c':
			read += std::to_string(*(int8_t*)tag_start);
			++tag_start;
			break;

		case 'C':
			read += std::to_string(*(uint8_t*)tag_start);
			++tag_start;
			break;

		case 's':
			read += std::to_string(*(int16_t*)tag_start);
			tag_start += 2;
			break;

		case 'S':
			read += std::to_string(*(uint16_t*)tag_start);
			tag_start += 2;
			break;

		case 'i':
			read += std::to_string(*(int32_t*)tag_start);
			tag_start += 4;
			break;

		case 'I':
			read += std::to_string(*(uint32_t*)tag_start);
			tag_start += 4;
			break;

		case 'f':
			read += std::to_string(*(float*)tag_start);
			tag_start += 4;
			break;

		case 'Z':
		{
			auto string_size_before = read.size();
			read += tag_start;
			auto string_size_after = read.size();
			tag_start += (string_size_after - string_size_before + 1);
			break;
		}
		case 'H':
		{
			auto string_size_before = read.size();
			read += tag_start;
			auto string_size_after = read.size();
			tag_start += (string_size_after - string_size_before + 1);
			break;
		}
		case 'B':
		{
			char sub_type = *tag_start;

			++tag_start;

			uint32_t count = *(uint32_t*)tag_start;

			tag_start += 4;

			switch (sub_type)
			{
			case 'c':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(int8_t*)tag_start);
					++tag_start;
				}
				break;

			case 'C':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(uint8_t*)tag_start);
					++tag_start;
				}
				break;

			case 's':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(int16_t*)tag_start);
					tag_start += 2;
				}
				break;

			case 'S':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(uint16_t*)tag_start);
					tag_start += 2;
				}
				break;

			case 'i':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(int32_t*)tag_start);
					tag_start += 4;
				}
				break;

			case 'I':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(uint32_t*)tag_start);
					tag_start += 4;
				}
				break;

			case 'f':
				for (uint32_t i = 0; i < count; ++i) {
					read += std::to_string(*(float*)tag_start);
					tag_start += 4;
				}
				break;

			default:
				break;
			}

			break;
		}
		default:
			break;
		}

	}

	return read;
};

bool BamFileProcessor::next() {

	this->advance(this->alignment_size_);

	NEED_BYTES(4);

	READ_VAR(this->alignment_size_);

	NEED_BYTES(this->alignment_size_);

	this->alignment_ = this->now();

	return true;
};