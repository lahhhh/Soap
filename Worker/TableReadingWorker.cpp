#include "TableReadingWorker.h"

#include "FileIO.h"

void TableReadingWorker::run() {
	CustomMatrix* res = read_table(this->file_name_, fast_);

	if (res == nullptr) {
		G_TASK_WARN("File reading failed.");
	}
	else {
		res->adjust_type();
		emit x_data_frame_ready(res);
	}

	G_TASK_END;
};
