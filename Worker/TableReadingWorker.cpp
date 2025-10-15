#include "TableReadingWorker.h"

#include "FileIO.h"

bool TableReadingWorker::work() {

	this->res_.reset(read_table(this->file_name_, fast_, skip_symbol_));

	if (this->res_ == nullptr) {
		G_TASK_WARN("File reading failed.");
		return false;
	}

	this->res_->adjust_type();

	return true;
};

void TableReadingWorker::run() {

	G_TASK_LOG("Start table reading...");

	if (!this->work()) {
		G_TASK_END;
	}

	G_TASK_LOG("Table reading finished.");

	emit x_data_frame_ready(this->res_.release());

	G_TASK_END;
};
