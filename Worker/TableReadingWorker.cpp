#include "TableReadingWorker.h"

#include "FileIO.h"

bool TableReadingWorker::work() {

	this->res_.reset(read_table(this->file_name_, fast_));

	if (this->res_ == nullptr) {
		G_TASK_WARN("File reading failed.");
		return false;
	}

	this->res_->adjust_type();

	return true;
};

void TableReadingWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_data_frame_ready(this->res_.release());

	G_TASK_END;
};
