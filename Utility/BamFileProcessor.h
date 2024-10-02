#pragma once

#include "Identifier.h"

#include <iostream>
#include <fstream>

#include <zlib.h>
#include <vector>
#include <string>


/*
* 
*	Example:
* 
	BamFileProcessor bfp(file_name);

	bfp.process_header();

	while (bfp.next()){
		...
	}
*/
class BamFileProcessor
{
public:

	explicit BamFileProcessor(const std::string& bam_file_name);

	~BamFileProcessor();

	z_stream stream_;

	std::ifstream in_;

	std::string bam_file_name_;

	std::string message_;

	long long file_size_ = 0;

	uint32_t alignment_size_ = 0;

	std::ptrdiff_t uncompressed_end_ = 0;

	std::ptrdiff_t reading_location_ = 0;

	std::vector<std::string> ref_names_;

	std::vector<uint32_t> ref_length_;

	char* uncompressed_ = nullptr;

	char* buffer_ = nullptr;

	char* alignment_ = nullptr;

	bool read_block();

	bool open_file();

	bool skip_header();

	bool record_reference();

	bool next();

	std::string get_read();

	bool process_header();

	// no check
	void read(void* dest, std::size_t length);

	// no check 
	void advance(std::size_t length);

	// no check
	bool need_bytes(std::size_t n_bytes);

	// no check
	std::string ref_name();

	// no check
	int32_t get_ref_id();

	// no check
	int32_t get_pos();

	// no check
	int32_t get_POS();

	uint8_t get_read_name_length();

	uint16_t get_flag();

	bool is_unmapped();

	bool get_strand_direction();

	uint32_t get_seq_length();

	int32_t get_next_ref_id();

	std::string get_next_ref();

	int32_t get_next_pos();

	int32_t get_next_POS();

	int32_t get_template_length();

	std::string get_read_name();

	std::string get_seq();

	std::string get_cigar_string();

	uint8_t get_mapq();

	uint16_t get_n_cigar_op();

	std::string get_tag_string();

	std::string get_base_quality();

	char* now();
};

