#pragma once

#include "Identifier.h"

#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include "custom.h"

class VcfProcessor
{
public:
	VcfProcessor(const std::string& bam_file_name, int cache_size = 10000);

	bool generate_snp_matrix(const std::string& output_name, int min_mutation_depth = 10);

private:

	bool check_header(const StringSplitView& view);

	std::string file_name_;

	std::ifstream in_;
	std::ofstream out_;

	std::vector<std::string> cache_;
	int cache_size_{ 10000 };

};

