#include "VelocytoWorker.h"

#include <QFile>
#include <QTextStream>

#include <fstream>

#include "BamFileProcessor.h"

bool VelocytoWorker::work() {

	this->create_index();

	soap::Species species;

	if (this->mode_ == WorkMode::SingleCellRna) {
		species = this->single_cell_rna_->species_;
	}
	else if (this->mode_ == WorkMode::SingleCellMultiome) {
		species = this->single_cell_multiome_->species_;
	}

	if (species == soap::Species::Human) {

		std::string hg38_path = (FILE_HUMAN_HG38_GTF).toStdString();
		if (!this->parse_gtf(hg38_path.data())) {
			return false;
		}
	}
	else {
		std::string mm10_path = (FILE_MOUSE_MM10_GTF).toStdString();
		if (!this->parse_gtf(mm10_path.data())) {
			return false;
		}
	}

	this->filter_transcript_model();

	if (species == soap::Species::Human) {
		std::string hg38_path = (FILE_HUMAN_HG38_MASK_GTF).toStdString();
		if (!this->parse_mask(hg38_path.data())) {
			return false;
		}
	}
	else {
		std::string mm10_path = (FILE_MOUSE_MM10_MASK_GTF).toStdString();
		if (!this->parse_mask(mm10_path.data())) {
			return false;
		}
	}

	this->filter_mask_intervals();

	/* skip mark up introns */

	int n_gene = this->gene_names_.size();
	int n_cell = this->cell_names_.size();

	this->spliced_counts_.resize(n_gene, n_cell);
	this->spliced_counts_.setZero();

	this->unspliced_counts_.resize(n_gene, n_cell);
	this->unspliced_counts_.setZero();

	if (!this->count_reads()) {
		return false;
	}

	this->umi_set_.clear();

	this->generate_count_matrix();

	return true;
};

void VelocytoWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_velocyto_ready(this->res_.release());

	G_TASK_END;

};

void VelocytoWorker::create_index() {

	if (this->mode_ == WorkMode::SingleCellRna) {

		int n_cell{ 0 };

		auto& metadata = this->single_cell_rna_->metadata()->mat_;
		if (metadata.data_type_.contains(METADATA_BARCODES) && metadata.data_type_.at(METADATA_BARCODES) == CustomMatrix::DataType::QString) {
			const auto& barcodes = metadata.get_const_qstring_reference(METADATA_BARCODES);
			n_cell = barcodes.size();
			for (int i = 0; i < n_cell; ++i) {
				this->barcode_index_[barcodes[i].toStdString()] = i;
			}
		}
		else {
			const auto& barcodes = this->single_cell_rna_->counts()->colnames_;
			n_cell = barcodes.size();
			for (int i = 0; i < n_cell; ++i) {
				this->barcode_index_[barcodes[i].toStdString()] = i;
			}
			G_TASK_NOTICE("Did not find barcodes metadata, use cell names instead");
		}

		this->umi_set_.resize(n_cell);

		this->cell_names_ = this->single_cell_rna_->counts()->colnames_;

		this->gene_names_ = this->single_cell_rna_->counts()->rownames_;

		int n_gene = this->gene_names_.size();

		for (int i = 0; i < n_gene; ++i) {
			this->gene_name_index_[this->gene_names_[i].toStdString()] = i;
		}

	}
	else if (this->mode_ == WorkMode::SingleCellMultiome) {

		int n_cell{ 0 };

		auto& metadata = this->single_cell_multiome_->metadata()->mat_;
		if (metadata.data_type_.contains(METADATA_BARCODES) && metadata.data_type_.at(METADATA_BARCODES) == CustomMatrix::DataType::QString) {
			const auto& barcodes = metadata.get_const_qstring_reference(METADATA_BARCODES);
			n_cell = barcodes.size();
			for (int i = 0; i < n_cell; ++i) {
				this->barcode_index_[barcodes[i].toStdString()] = i;
			}
		}
		else {
			const auto& barcodes = this->single_cell_multiome_->rna_counts()->colnames_;
			n_cell = barcodes.size();
			for (int i = 0; i < n_cell; ++i) {
				this->barcode_index_[barcodes[i].toStdString()] = i;
			}
			G_TASK_NOTICE("Did not find barcodes metadata, use cell names instead");
		}

		this->umi_set_.resize(n_cell);

		this->cell_names_ = this->single_cell_multiome_->rna_counts()->colnames_;

		this->gene_names_ = this->single_cell_multiome_->rna_counts()->rownames_;

		int n_gene = this->gene_names_.size();

		for (int i = 0; i < n_gene; ++i) {
			this->gene_name_index_[this->gene_names_[i].toStdString()] = i;
		}

	}
};

void VelocytoWorker::generate_count_matrix() {

	this->res_.reset(new VelocytoBase());

	int n_gene = this->gene_names_.size(), n_cell = this->cell_names_.size();

	auto& spliced = SUBMODULES(*this->res_, SparseInt)[VARIABLE_RNA_SPLICED];

	spliced.data_type_ = SparseInt::DataType::Spliced;

	spliced.colnames_ = this->cell_names_;
	spliced.rownames_ = this->gene_names_;

	G_TASK_LOG(QString::number(this->spliced_counts_.sum()) + " spliced reads found.");

	spliced.mat_ = this->spliced_counts_.sparseView();
	this->spliced_counts_.resize(0, 0);

	auto& unspliced = SUBMODULES(*this->res_, SparseInt)[VARIABLE_RNA_UNSPLICED];

	unspliced.data_type_ = SparseInt::DataType::Unspliced;

	unspliced.colnames_ = this->cell_names_;
	unspliced.rownames_ = this->gene_names_;

	G_TASK_LOG(QString::number(this->unspliced_counts_.sum()) + " unspliced reads found.");
	
	unspliced.mat_ = this->unspliced_counts_.sparseView();
	this->unspliced_counts_.resize(0, 0);
};

// return value - 0 : no map; 1 : unspliced; 2 : spliced; 3 : multi match
// TODO : rewrite with std::upper_bound
int VelocytoWorker::count(
	const std::string& chromosome_name, 
	const std::vector<int>& segments, 
	std::string* ptr_gene_name)
{
	const int start = segments.front(), end = segments.back();

	auto iter = this->transcript_models_.find(chromosome_name);

	if (iter == this->transcript_models_.end()) {
		return 0;
	}

	/**********  start align section  ************/

	std::string gene_name;

	const auto& models = iter->second;
	const auto model_size = models.size();

	int strand_start = models[0].start_;
	int strand_last_start = models[model_size - 1].start_;
	int strand_end = models[model_size - 1].end_;

	int count_res{ 0 };

	if (end <= strand_start) {
		return 0;
	}
	else if (start < strand_start) {
		*ptr_gene_name = models.front().gene_name_;

		auto count_res = models.front().match(segments);

		return count_res;
	}
	else if (start >= strand_last_start) {
		if (start >= strand_end) {
			return 0;
		}
		else {

			int current_model_loc = model_size - 1;
			int current_model_end = models[current_model_loc].end_;
			int n_search{ 0 };

			// peek backward
			while (current_model_end > start || n_search < 100) {

				++n_search;

				auto match_res = models[current_model_loc].match(segments);
				switch (match_res)
				{
				case 0:
					break;
				case 1:
					if (gene_name.empty()) {
						gene_name = models[current_model_loc].gene_name_;
					}
					else if (gene_name != models[current_model_loc].gene_name_) {
						return 3;
					}
					*ptr_gene_name = models[current_model_loc].gene_name_;
					count_res = 1;
					break;
				case 2:
					*ptr_gene_name = models[current_model_loc].gene_name_;
					return 2;
					break;
				default:
					break;
				}

				--current_model_loc;

				if (current_model_loc < 0) {
					break;
				}

				current_model_end = models[current_model_loc].end_;
			}

			return count_res;
		}
	}


	if (model_size == 1) {
		auto match_res = models[0].match(segments);
		switch (match_res)
		{
		case 0:
			return 0;
			break;
		case 1:
			*ptr_gene_name = models[0].gene_name_;
			return 1;
			break;
		case 2:
			*ptr_gene_name = models[0].gene_name_;
			return 2;
			break;
		default:
			break;
		}
	}

	int low = 0;
	int high = model_size - 1;
	int middle;

	const auto * const data = models.data();

	while (low <= high) {

		middle = (high + low) / 2;

		const auto* local_ptr = data + middle;

		int local_start = local_ptr->start_, next_start = (local_ptr + 1)->start_;

		if (start < local_start) {
			high = middle - 1;
		}
		else if (start >= next_start) {
			low = middle + 1;
		}
		else {

			int current_model_loc = middle;
			int current_model_end = models[current_model_loc].end_;
			int n_search{ 0 };

			// peek around for best match
			// peek backward
			while (current_model_end > start || n_search < 100) {

				auto&& current_model = models[current_model_loc];

				++n_search;

				auto match_res = current_model.match(segments);
				switch (match_res)
				{
				case 0:
					break;
				case 1:
					if (gene_name.empty()) {
						gene_name = models[current_model_loc].gene_name_;
					}
					else if (gene_name != models[current_model_loc].gene_name_) {
						return 3;
					}
					*ptr_gene_name = current_model.gene_name_;
					count_res = 1;
					break;
				case 2:
					*ptr_gene_name = current_model.gene_name_;
					return 2;
					break;
				default:
					break;
				}

				--current_model_loc;

				if (current_model_loc < 0) {
					break;
				}

				current_model_end = models[current_model_loc].end_;
			}

			current_model_loc = middle + 1;
			int current_model_start = models[current_model_loc].start_;

			n_search = 0;

			// peek forward
			while (current_model_start < end || n_search < 100) {

				auto&& current_model = models[current_model_loc];

				++n_search;

				auto match_res = current_model.match(segments);
				switch (match_res)
				{
				case 0:
					break;
				case 1:
					if (gene_name.empty()) {
						gene_name = models[current_model_loc].gene_name_;
					}
					else if (gene_name != models[current_model_loc].gene_name_) {
						return 3;
					}
					*ptr_gene_name = current_model.gene_name_;
					count_res = 1;
					break;
				case 2:
					*ptr_gene_name = current_model.gene_name_;
					return 2;
					break;
				default:
					break;
				}

				++current_model_loc;

				if (current_model_loc >= model_size) {
					break;
				}

				current_model_start = models[current_model_loc].start_;
			}

			return count_res;
		}
	}

};

bool VelocytoWorker::interval_masked(const std::string& chromosome_name, int32_t start, int32_t end) {

	auto iter = this->mask_intervals_.find(chromosome_name);

	if (iter == this->mask_intervals_.cend()) {
		return false;
	}

	const auto& intervals = iter->second;

	const auto interval_size = intervals.size() / 2;

	if (interval_size == 0)return false;

	int interval_start = intervals[0];
	int interval_last_start = intervals[2 * interval_size - 2];
	int interval_end = intervals[2 * interval_size - 1];

	if (end <= interval_start) {
		return false;
	}
	else if (start < interval_start) {
		return false;
	}
	else if (start >= interval_last_start) {
		if (end <= interval_end) {
			return true;
		}
		else{
			return false;
		}
	}

	if (interval_size == 1) {
		if (start >= interval_start && end <= interval_end) {
			return true;
		}
		else {
			return false;
		}
	}

	int low = 0;
	int high = interval_size - 1;
	int middle;

	const int* data = intervals.data();

	while (low <= high) {

		middle = (high + low) / 2;

		const int* local_ptr = data + 2 * middle;

		int local_start = *local_ptr, next_start = *(local_ptr + 2);

		if (start < local_start) {
			high = middle - 1;
		}
		else if (start >= next_start) {
			low = middle + 1;
		}
		else {

			int local_end = *(local_ptr + 1);

			if (start <= local_end) {
				if (end <= local_end) {
					return true;
				}
				else {
					return false;
				}
			}
			else {
				return false;
			}
		}
	}
};

bool VelocytoWorker::get_umi_and_barcode(const char* read, int alignment_size, int* umi, std::string* barcode) {

	

	bool umi_found{ false }, barcode_found{ false };

	uint8_t read_name_length{ 0 };

	std::memcpy(&read_name_length, read + 8, 1);

	uint16_t n_cigar_op{ 0 };

	std::memcpy(&n_cigar_op, read + 12, 2);

	uint32_t seq_length{ 0 };

	std::memcpy(&seq_length, read + 16, 4);

	const char* tag_start = read + 32 + read_name_length + 4 * n_cigar_op + (3 * seq_length + 1) / 2;

	const char* end = read + alignment_size;

	while (tag_start < end) {

		if (umi_found && barcode_found) {
			return true;
		}

		if (*tag_start == this->umi_tag_[0] && *(tag_start + 1) == this->umi_tag_[1]) {
			tag_start += 2;

			if (*tag_start != 'Z') {
				return false;
			}

			++tag_start;

			int umi_number{ 0 };

			while (*tag_start != '\0') {
				umi_number <<= 2;
				
				switch (*tag_start)
				{
				case 'A':
					break;
				case 'C':
					umi_number += 1;
					break;
				case 'G':
					umi_number += 2;
					break;
				case 'T':
					umi_number += 3;
					break;
				default:
					break;
				}

				++tag_start;
			}

			*umi = umi_number;

			umi_found = true;

			++tag_start;

			continue;
		}
		else if (*tag_start == this->barcode_tag_[0] && *(tag_start + 1) == this->barcode_tag_[1]) {
			tag_start += 2;

			if (*tag_start != 'Z') {
				return false;
			}

			++tag_start;

			barcode->assign(tag_start);

			while (*tag_start != '\0') {
				++tag_start;
			}

			++tag_start;

			barcode_found = true;

			continue;

		}

		tag_start += 2;

		char value_type = *tag_start;

		++tag_start;

		switch (value_type)
		{
		case 'A':
			++tag_start;
			break;

		case 'c':
			++tag_start;
			break;

		case 'C':
			++tag_start;
			break;

		case 's':
			tag_start += 2;
			break;

		case 'S':
			tag_start += 2;
			break;

		case 'i':
			tag_start += 4;
			break;

		case 'I':
			tag_start += 4;
			break;

		case 'f':
			tag_start += 4;
			break;

		case 'Z':
		{
			while (*tag_start != '\0') {
				++tag_start;
			}
			++tag_start;
			break;
		}
		case 'H':
		{
			while (*tag_start != '\0') {
				++tag_start;
			}
			++tag_start;
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
				tag_start += count;
				break;

			case 'C':
				tag_start += count;
				break;

			case 's':
				tag_start += count * 2;
				break;

			case 'S':
				tag_start += count * 2;
				break;

			case 'i':
				tag_start += count * 4;
				break;

			case 'I':
				tag_start += count * 4;
				break;

			case 'f':
				tag_start += count * 4;
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

	return (umi_found && barcode_found);
};


int VelocytoWorker::map_read(const char* read, int alignment_size, const std::string& ref_name) {

	int32_t pos{ 0 };

	std::memcpy(&pos, read + 4, 4);

	++pos;// 1-based

	uint8_t read_name_length{ 0 };

	std::memcpy(&read_name_length, read + 8, 1);

	uint16_t n_cigar_op{ 0 };

	std::memcpy(&n_cigar_op, read + 12, 2);

	const char* cigar_start = read + 32 + read_name_length;

	std::string barcode;

	int umi{ 0 };

	if (!this->get_umi_and_barcode(read, alignment_size, &umi, &barcode)) {
		return 1;
	}

	static std::vector<int> matched_segment;

	matched_segment.clear();

	bool masked{ true }, test_mask{ true };

	for (uint16_t i = 0; i < n_cigar_op; ++i) {

		uint32_t op = *(uint32_t*)(cigar_start + 4 * i);

		uint32_t op_length = op >> 4;

		if (op_length > 3000000) {
			return 2;
		}

		uint32_t op_type = op & 0xf;

		switch (op_type)
		{
		case 0:
			/*M*/
		{
			if (test_mask) {
				if (!this->interval_masked(ref_name, pos, pos + op_length)) {
					masked = false;
					test_mask = false;
				}
			}
			matched_segment.push_back(pos);
			matched_segment.push_back(pos + op_length);
			pos += op_length;

			break;
		}
		case 2:
			/*D*/
			pos += op_length;
			break;
		case 3:
			/*N*/
			pos += op_length;
			break;
		case 7:
			/*=*/
		{
			if (test_mask) {
				if (!this->interval_masked(ref_name, pos, pos + op_length)) {
					masked = false;
					test_mask = false;
				}
			}
			matched_segment.push_back(pos);
			matched_segment.push_back(pos + op_length);
			pos += op_length;

			break;
		}
		case 8:
			/*X*/
			pos += op_length;
			break;
		default:
			break;
		}

	}

	if (matched_segment.empty()) {
		return 2;// 0
	}

	if (masked) {
		return 3;
	}

	auto iter = this->barcode_index_.find(barcode);

	if (iter == this->barcode_index_.end()) {
		return 4;
	}

	int barcode_index = iter->second;

	std::string gene_name;
	auto count_res = this->count(ref_name, matched_segment, &gene_name);

	if (count_res == 0) {
		return 5;
	}

	if (count_res == 3) {
		return 8;
	}

	auto iter2 = this->gene_name_index_.find(gene_name);

	if (iter2 == this->gene_name_index_.end()) {
		return 6;
	}

	int gene_index = iter2->second;

	auto& umi_set = this->umi_set_[barcode_index];

	if (!umi_set.insert(umi).second) {
		return 7;
	}

	if (count_res == 2) {
		++ this->spliced_counts_(gene_index, barcode_index);
	}
	else {
		++ this->unspliced_counts_(gene_index, barcode_index);
	}

	return 0;
};


bool VelocytoWorker::count_reads() {

	std::size_t read_count{ 0 };

	for (const auto& bam_file : this->bam_files_) {
		BamFileProcessor bfp(bam_file.toStdString());

		bfp.process_header();

		while (bfp.next()) {

			++read_count;

			if (read_count % 10000000 == 0) {
				G_TASK_LOG(QString::number(read_count) + " reads has been processed.");
			}

			if (bfp.get_mapq() == 0) {
				continue;
			}

			auto ref_name = bfp.ref_name();

			if (ref_name.size() > 3 && strnicmp(ref_name.data(), "chr", 3) == 0) {
				ref_name = ref_name.substr(3);
			}

			this->map_read(bfp.alignment_, bfp.alignment_size_, ref_name);
		}
				
	}

	return true;

};

bool VelocytoWorker::filter_mask_intervals() {

	for (auto iter = this->mask_intervals_.begin(); iter != this->mask_intervals_.end(); ++iter) {
		
		auto& chr = iter->second;
		const auto size = chr.size();
		const auto n_interval = size / 2;

		std::vector<int> interval_starts(size / 2);

		for (std::size_t i = 0; i < n_interval; ++i) {
			interval_starts[i] = chr[i * 2];
		}

		const auto order = custom::order(interval_starts);		

		std::vector<int> new_chr(size);

		for (std::size_t i = 0; i < n_interval; ++i) {
			new_chr[i * 2] = chr[order[i] * 2];
			new_chr[i * 2 + 1] = chr[order[i] * 2 + 1];
		}

		chr.clear();

		int current_start = new_chr[0];
		int current_end = new_chr[1];

		constexpr int tolerance = 5;

		for (std::size_t i = 2; i < size; i += 2) {
			if (new_chr[i] <= current_end + tolerance) {
				current_end = new_chr[i + 1];
			}
			else {
				chr.push_back(current_start);
				chr.push_back(current_end);
				current_start = new_chr[i];
				current_end = new_chr[i + 1];
			}
		}

		chr.push_back(current_start);
		chr.push_back(current_end);
	}
	return true;
};

bool VelocytoWorker::filter_transcript_model() {

	QList<std::string> target_gene_names;

	if (this->mode_ == WorkMode::SingleCellRna) {
		target_gene_names = custom::sapply(this->single_cell_rna_->counts()->rownames_, [](const QString& rowname) {return rowname.toStdString(); });
	}
	else {
		target_gene_names = custom::sapply(this->single_cell_multiome_->rna_counts()->rownames_, [](const QString& rowname) {return rowname.toStdString(); });
	}

	for (auto iter = this->transcript_models_.begin(); iter != this->transcript_models_.end();) {

		auto chromosome_gene_names = custom::sapply(iter->second, [](const TranscriptModel& tm) {return tm.gene_name_; });

		auto filter = custom::in(chromosome_gene_names, target_gene_names);

		if (filter.count() == 0) {
			iter = this->transcript_models_.erase(iter);
		}
		else {
			iter->second = custom::sliced(iter->second, filter);
			std::sort(iter->second.begin(), iter->second.end());
			++iter;
		}
	}

	return !this->transcript_models_.empty();
};

bool VelocytoWorker::parse_mask(const char* file_name) {
	FILE* file = fopen(file_name, "r");

	if (file == NULL) {
		G_TASK_WARN("mask file open failed.");
		return false;
	}

	char* buffer = new char[1024];// mask line length should be ~= 100 bytes	

	if (fgets(buffer, 1024, file) == NULL) {
		delete[]buffer;
		G_TASK_WARN("mask file is empty.");
		return false;
	}

	while (buffer[0] == '#') {
		if (fgets(buffer, 1024, file) == NULL) {
			delete[]buffer;
			G_TASK_WARN("mask file is broken.");
			return false;
		}
	}

	const char* content = buffer;

	uint16_t tab_count{ 0 }, line_index{ 0 };

	while (*content != '\0') {
		if (*content == '\t') {
			++tab_count;
		}

		++content;
	}

	if (tab_count != 8) {
		delete[]buffer;
		G_TASK_WARN("the format of mask file is wrong.");
		return false;
	}

	std::string current_chromosome;

	char current_strand{ '\0' };

	int16_t index[8];

	int current_start{ 0 }, current_end{ 0 };

	do {

		tab_count = 0, line_index = 0;

		const char* line_data = buffer;

		while (*line_data != '\0') {
			if (*line_data == '\t') {
				index[tab_count++] = line_index + 1;
			}
			++line_index;
			++line_data;
		}

		current_start = custom::atoi_specialized(buffer + index[2]);
		current_end = custom::atoi_specialized(buffer + index[3]);

		if (strnicmp(buffer, "chr", 3) == 0) {
				auto& chromosome = this->mask_intervals_[std::string(buffer + 3, index[0] - 4)];

				chromosome.push_back(current_start);
				chromosome.push_back(current_end);
		}
		else {
			auto& chromosome = this->mask_intervals_[std::string(buffer, index[0] - 1)];

			chromosome.push_back(current_start);
			chromosome.push_back(current_end);
		}

	} while (fgets(buffer, 1024, file) != NULL);

	fclose(file);

	delete[] buffer;
	return true;
}

bool VelocytoWorker::parse_gtf(const char* file_name) {

	FILE* file = fopen(file_name, "r");

	if (file == NULL) {
		G_TASK_WARN("gtf file open failed.");
		return false;
	}

	char* buffer = new char[1024];// gtf line length should be ~= 100 bytes	

	if (fgets(buffer, 1024, file) == NULL) {
		delete[]buffer;
		G_TASK_WARN("gtf file is empty.");
		return false;
	}

	while (buffer[0] == '#') {
		if (fgets(buffer, 1024, file) == NULL) {
			delete[]buffer;
			G_TASK_WARN("gtf file is broken.");
			return false;
		}
	}

	const char* content = buffer;

	uint16_t tab_count{ 0 }, line_index{ 0 };

	while (*content != '\0') {
		if (*content == '\t') {
			++tab_count;
		}

		++content;
	}

	if (tab_count != 8) {
		delete[]buffer;
		G_TASK_WARN("the format of gtf file is wrong.");
		return false;
	}


	TranscriptModel current_model;

	std::string current_chromosome;

	char current_strand{ '\0' };

	int16_t index[8];

	do {

		tab_count = 0, line_index = 0;

		const char* line_data = buffer;

		while (*line_data != '\0') {
			if (*line_data == '\t') {
				index[tab_count++] = line_index + 1;
			}
			++line_index;
			++line_data;
		}

		if (strncmp(buffer + index[1], "transcript\t", 11) == 0) {
			if (!current_model.empty()) {
				current_model.finalize();
				this->transcript_models_[current_chromosome].emplace_back(current_model);
			}

			for (uint16_t i = index[7]; i < line_index - 13; ++i) {
				if (strncmp(buffer + i, "gene_name", 9) == 0) {
					uint16_t gene_name_start = i + 11;
					uint16_t j = gene_name_start + 1;

					while (buffer[j] != '\"') {
						++j;
					}

					current_model.reset(std::string(buffer + gene_name_start, j - gene_name_start));

					break;
				}
			}

			if (strnicmp(buffer, "chr", 3) == 0) {
				current_chromosome = std::string(buffer + 3, index[0] - 4);
			}
			else {
				current_chromosome = std::string(buffer, index[0] - 1);
			}

			current_strand = buffer[index[5]];
		}
		else if (strncmp(buffer + index[1], "exon\t", 5) == 0) {
			current_model.append_exon(custom::atoi_specialized(buffer + index[2]), custom::atoi_specialized(buffer + index[3]));
		}

	} while (fgets(buffer, 1024, file) != NULL);

	if (!current_model.empty()) {
		current_model.finalize();
		this->transcript_models_[current_chromosome].emplace_back(current_model);
	}

	fclose(file);

	delete[]buffer;

	return true;
};