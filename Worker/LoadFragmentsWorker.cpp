#include "LoadFragmentsWorker.h"

#include <zlib.h>

LoadFragmentsWorker::LoadFragmentsWorker(
	const QStringList& barcodes,
	const QStringList& fragments_files
) :
	barcodes_(barcodes),
	fragments_files_(fragments_files)
{};

bool LoadFragmentsWorker::work() {

	this->create_index();

	if (!this->load_fragments()) {
		G_TASK_WARN("Fragments Loading failed.");
		return false;
	}

	return true;
};

void LoadFragmentsWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_fragments_ready(this->fragments_.release());

	G_TASK_LOG("Fragments loading finished.");
	G_TASK_END;
}

void LoadFragmentsWorker::create_index() {

	const qsizetype size = this->barcodes_.size();

	this->barcodes_index_.reserve(size);
	for (qsizetype i = 0; i < size; ++i) {
		this->barcodes_index_[this->barcodes_[i]] = i;
	}

	this->n_cell_ = size;
};

bool LoadFragmentsWorker::load_fragments() {

	this->fragments_.reset(new Fragments());

	int buffer_length = 256;
	char buffer[256];

	std::size_t process = 0;

	for (const auto& fragments_file : this->fragments_files_) {
		gzFile fragments = gzopen_w((const wchar_t*)fragments_file.utf16(), "rb");

		if (fragments == NULL) {
			G_TASK_WARN("File : " + fragments_file + " is broken.");
			return false;
		}

		bool not_end_of_file = false;
		while (not_end_of_file = gzgets(fragments, buffer, buffer_length) != 0)
		{
			if (buffer[0] != '#')break;
		}
		if (!not_end_of_file) {
			gzclose(fragments);
			G_TASK_WARN("File : " + fragments_file + " is broken.");
			return false;
		}
		QString barcode;
		do {
			const char* c = buffer, * d, * seq_end;
			while (*c != '\t') {
				++c;
			}
			seq_end = c;

			int start = custom::atoi_specialized(&c);
			int end = custom::atoi_specialized(&c);

			if (start >= end) {
				G_TASK_WARN("Invalid Fragments File: " + QString::fromUtf8(buffer));
				gzclose(fragments);
				return false;
			}

			++c;
			d = c;
			while (*c != '\t') {
				++c;
			}
			++c;
			barcode = QString::fromUtf8(d, c - d - 1);

			auto iter = this->barcodes_index_.find(barcode);
			if (iter != this->barcodes_index_.end()) {

				int index = iter->second;
				QString seq_name = custom::standardize_chromosome_name(QString::fromUtf8(buffer, seq_end - buffer));
				auto& sequence_matrix = this->fragments_->data_[seq_name];

				if (sequence_matrix.empty()) {
					sequence_matrix.resize(this->n_cell_);
				}

				auto& barcode_fragments = sequence_matrix[index];
				barcode_fragments.first.emplace_back(start);
				barcode_fragments.second.emplace_back(end);
			}

			if (process % 50000000 == 0) {
				G_TASK_LOG(QString::number(process) + " sequences has been processed.");
			}

			++process;
		} while (not_end_of_file = gzgets(fragments, buffer, buffer_length) != 0);

		gzclose(fragments);
	}

	this->fragments_->finalize();

	return true;
};