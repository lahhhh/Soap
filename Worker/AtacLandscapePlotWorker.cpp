#include "AtacLandscapePlotWorker.h"

#include "Custom.h"

bool AtacLandscapePlotWorker::work() {

	if (!this->create_index1()) {
		return false;
	}

	if (!this->create_index2()) {
		return false;
	}

	if (!this->calculate_matrix()) {
		return false;
	}

	return true;
};

void AtacLandscapePlotWorker::run() {

	G_TASK_LOG("Start calculating ATAC landscape...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_plot_ready(ATAC_LANDSCAPE_PLOT_ELEMENTS{ this->plot_mat_, this->cell_loc_, this->chr_loc_, this->scale_ });

	G_TASK_END;
};

bool AtacLandscapePlotWorker::calculate_matrix() {

	bool normalize_by_depth = this->normalize_method_ == "Sequencing Depth";
	bool normalize_by_chromosome = this->normalize_method_ == "Choosed Chromosomes";

	QVector<double> depth;

	for (auto&& [chr_name, chr_frag] : *this->fragments_) {

		if (normalize_by_depth) {
			if (depth.isEmpty()) {
				depth = custom::sapply(chr_frag, [](auto&& f) {return 2.0 * f.first.size(); });
			}
			else {
				depth = custom::add(depth, custom::sapply(chr_frag, [](auto&& f) {return 2.0 * f.first.size(); }));
			}
		}

		if (!this->chr_index_.contains(chr_name)) {
			continue;
		}

		int chr_size = this->chr_sizes_[chr_name];
		int chr_start = this->chr_index_[chr_name];

		int n_cell = chr_frag.size();

		for (int i = 0; i < n_cell; ++i) {
			int cell_ind = this->cell_index_[i];

			if (cell_ind == -1) {
				continue;
			}

			auto&& start = chr_frag[i].first;
			auto&& end = chr_frag[i].second;

			for (auto loc : start) {
				if (loc > chr_size || loc < 1) {
					G_TASK_WARN("Illegal fragments position:" + chr_name + ":" + QString::number(loc));
					return false;
				}

				++this->plot_mat_(cell_ind, chr_start + (loc - 1) / this->resolution_);

			}

			for (auto loc : end) {
				if (loc > chr_size || loc < 1) {
					G_TASK_WARN("Illegal fragments position:" + chr_name + ":" + QString::number(loc));
					return false;
				}

				++this->plot_mat_(cell_ind, chr_start + (loc - 1) / this->resolution_);

			}
		}
	}

	if (normalize_by_depth) {

		int n_cell = depth.size();
		int n_cell_use = this->plot_mat_.rows();

		Eigen::ArrayXd depth_use = Eigen::ArrayXd::Zero(n_cell_use);
		for (int i = 0; i < n_cell; ++i) {
			int ind = this->cell_index_[i];

			if (ind != -1) {
				depth_use[ind] = depth[i];
			}
		}

		double mean_depth = depth_use.mean();
		depth_use = mean_depth / depth_use;

		this->plot_mat_.array().colwise() *= depth_use;
	}

	if (normalize_by_chromosome) {
		Eigen::ArrayXd multiplier = this->plot_mat_.rowwise().sum();
		double mean_frag = multiplier.mean();
		multiplier = mean_frag / multiplier;

		this->plot_mat_.array().colwise() *= multiplier;
	}

	if (this->log_transform_) {
		this->plot_mat_ = (this->plot_mat_.array() + 1.0).log();
	}

	if (this->scale_) {
		custom::scale_in_place(this->plot_mat_);
	}
};

bool AtacLandscapePlotWorker::create_index2() {
	QFile file(FILE_CICERO_HG38_CHROMSIZE);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		G_TASK_WARN("Can not open resource");
		return false;
	}

	QTextStream in(&file);

	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList content = line.split('\t');
		QString chr_name = custom::standardize_chromosome_name(content[0]);
		int chr_size = content[1].toInt();

		if (chr_size <= 1) {
			G_TASK_WARN("Broken Resource File.");
			return false;
		}

		this->chr_sizes_[chr_name] = chr_size;
	}

	int ind_now{ 0 };

	for (auto&& n : this->chromosomes_) {
		if (!this->chr_sizes_.contains(n)) {
			G_TASK_WARN("Error when loading Chromosome Size.");
			return false;
		}

		int size = this->chr_sizes_[n];

		int n_block = std::ceil((double)size / this->resolution_);

		this->chr_index_[n] = ind_now;

		this->chr_loc_ << std::make_tuple(n, ind_now, n_block);

		ind_now += n_block;
	}

	this->ncol_ = ind_now;

	this->plot_mat_ = Eigen::MatrixXd::Zero(this->nrow_, this->ncol_);

	if (ind_now == 0) {
		G_TASK_WARN("Error when loading chromosome size.");
		return false;
	}
	else {
		return true;
	}

}

bool AtacLandscapePlotWorker::create_index1() {

	int ind_now{ 0 };
	int n_cell = this->fragments_->cell_names_.size();
	this->cell_index_.resize(n_cell, -1);

	for (auto&& level : this->levels_) {

		auto ind = custom::match(this->factors_, level);

		if (ind.isEmpty()) {
			continue;
		}

		int ind_size = ind.size();

		this->cell_loc_ << std::make_tuple(level, ind_now, ind_size);

		for (int i = 0; i < ind_size; ++i) {
			this->cell_index_[ind[i]] = ind_now++;
		}
	}

	this->nrow_ = ind_now;

	if (ind_now == 0) {

		G_TASK_WARN("Find no cell data.");

		return false;
	}
	else {
		return true;
	}
};