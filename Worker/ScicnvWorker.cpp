#include "ScicnvWorker.h"

#include "FileIO.h"

#include "Custom.h"

bool ScicnvWorker::filter_data() {
	Eigen::ArrayX<bool> filter = custom::row_count<double, true, false>(this->mat_.mat_, 0).cast<double>() > (2. / this->n_cell_);
	this->mat_.row_slice(filter);

	return this->mat_.rows() > 0;
};

bool ScicnvWorker::sort_expression_by_chromosome() {

	std::unique_ptr<CustomMatrix> gene_location;

	if (this->species_ == soap::Species::Human) {
		gene_location.reset(read_sv(FILE_CNV_GENE_LOCATION_HUMAN, '\t'));
	}
	else if (this->species_ == soap::Species::Mouse) {
		gene_location.reset(read_sv(FILE_CNV_GENE_LOCATION_MOUSE, '\t'));
	}
	else {
		G_TASK_WARN("InferCNV now only support human and mouse data.");
		return false;
	}

	if (gene_location == nullptr) {
		G_TASK_WARN("Resource File Loading Failed.");
		return false;
	}

	QStringList gene_list = gene_location->get_qstring("Column 1");
	QStringList chromosome_list = gene_location->get_qstring("Column 2");
	std::ranges::for_each(chromosome_list, [](auto&& t) {t = custom::standardize_chromosome_name(t); });

	QStringList unique_chromosome = custom::unique(chromosome_list);	

	if (this->species_ == soap::Species::Human) {
		unique_chromosome = custom::sliced(soap::HumanChromosomeOrder, custom::in(soap::HumanChromosomeOrder, unique_chromosome));
	}
	else if (this->species_ == soap::Species::Mouse) {
		unique_chromosome = custom::sliced(soap::MouseChromosomeOrder, custom::in(soap::MouseChromosomeOrder, unique_chromosome));
	}

	QVector<int> gene_index;

	for (const auto& chr : unique_chromosome) {
		QStringList chr_gene_list = custom::sliced(gene_list, custom::equal(chromosome_list, chr));
		auto chr_gene_index = custom::valid_index_of(chr_gene_list, this->mat_.rownames_);
		int n_gene = chr_gene_index.size();

		if (n_gene > 20) {

			this->chromosome_location_.emplace_back(chr, gene_index.size(), n_gene);

			gene_index << chr_gene_index;
		}
	}

	if (gene_index.isEmpty()) {
		return false;
	}

	this->mat_.row_reorder(gene_index);

	this->resolution_ = (int)(this->mat_.rows() / (this->sharpness_ * 50));
	this->n_gene_ = this->mat_.rows();

	return true;
};

void ScicnvWorker::compute_reference() {

	if (this->reference_.isEmpty()) {
		this->reference_expression_ = custom::row_mean(this->mat_.mat_);
	}
	else {
		this->reference_expression_ = custom::col_sliced_row_mean(this->mat_.mat_, custom::in(this->metadata_, this->reference_));
	}
};

void ScicnvWorker::compute_reference_moving_average(bool sharp) {

	this->reference_moving_average_.resize(this->n_gene_);

	int half_width = ceil(this->resolution_ / 2.);
	if (sharp) {
		half_width = ceil(half_width / (5 * this->sharpness_));
	}

	for (const auto& [chr, start, n] : this->chromosome_location_) {
		int end = start + n;
		for (int i = start; i < end; ++i) {
			int gene_start = i - half_width > start ? i - half_width : start;
			int gene_end = i + half_width < end ? i + half_width : end;
			this->reference_moving_average_[i] = this->reference_expression_.segment(gene_start, gene_end - gene_start).mean();
		}
	}
};

void ScicnvWorker::compute_moving_average(bool sharp) {

	this->moving_average_.resize(this->n_gene_, this->n_cell_);

	int half_width = ceil(this->resolution_ / 2.);
	if (sharp) {
		half_width = ceil(half_width / (5 * this->sharpness_));
	}

	const int n_chromosome = this->chromosome_location_.size();

#pragma omp parallel for
	for(int i = 0; i < n_chromosome; ++i){
		const auto &[chr, start, n] = this->chromosome_location_[i];

		const int end = start + n;

		Eigen::ArrayXd col_sum;
		for (int j = start; j < end; ++j) {

			int gene_start = j - half_width > start ? j - half_width : start;
			int gene_end = j + half_width < end ? j + half_width : end;

			if (j == start) {
				col_sum = custom::row_reorder_and_column_sum(this->mat_.mat_, custom::seq_n(gene_start, gene_end - gene_start));
			}
			else {
				if (j - half_width > start) {
					col_sum -= this->mat_.mat_.row(gene_start - 1);
				}

				if (j + half_width <= end) {
					col_sum += this->mat_.mat_.row(gene_end - 1);
				}
			}

			this->moving_average_.row(j) = col_sum / (gene_end - gene_start);
		}
	}
};

Eigen::ArrayXd ScicnvWorker::digitalize(const Eigen::ArrayXd& arr) {
	int size = arr.size();
	Eigen::ArrayXd ret(size);
	for (int i = 0; i < size; ++i) {
		double val = arr[i];
		if (val > 0) {
			ret[i] = 1;
		}
		else if (val < 0) {
			ret[i] = -1;
		}
		else {
			ret[i] = 0;
		}
	}
	return ret;
};

void ScicnvWorker::compute_W() {

	this->W_.resize(this->n_gene_, this->n_cell_);

#pragma omp parallel for
	for(int i = 0; i < this->n_cell_; ++i){
		Eigen::ArrayXd w = (this->moving_average_.col(i) - this->reference_moving_average_) / (this->moving_average_.col(i) + this->reference_moving_average_ + 0.00001);
		w -= custom::median(w);
		w = w * w * w * w * 3.3222 + w * w * w * 5.6399 + w * w * 4.2189 + 3.8956 * w;
		this->W_.col(i) = w;
	}
};

void ScicnvWorker::compute_U() {
	this->U_.resize(this->n_gene_, this->n_cell_);

#pragma omp parallel for
	for (int i = 0; i < this->n_cell_; ++i) {
		Eigen::ArrayXd w = (this->moving_average_.col(i) - this->reference_moving_average_) / (this->moving_average_.col(i) + this->reference_moving_average_ + 0.00001);
		w = custom::cumsum(digitalize(w - custom::median(w)));
		this->U_.col(i) = w;
	};
};

void ScicnvWorker::compute_V() {
	this->V_.resize(this->n_gene_, this->n_cell_);

#pragma omp parallel for
	for (int i = 0; i < this->n_cell_; ++i) {
		Eigen::ArrayXd w = this->mat_.mat_.col(i);
		w = (w - this->reference_expression_) / (w + this->reference_expression_ + 0.00001);
		double d = ((w > 0).count() - (w < 0).count()) / (double)this->n_gene_;
		w = custom::cumsum(digitalize(w - custom::median(w)) - d);
		this->V_.col(i) = w;
	};

	double scale = (this->U_.maxCoeff() - this->U_.minCoeff()) / (this->V_.maxCoeff() - this->V_.minCoeff());
	this->V_ *= scale;
	this->U_ = 0.66 * this->U_ + 0.34 * this->V_;
};

Eigen::ArrayXd ScicnvWorker::slope(const Eigen::ArrayXXd& arr) {
	int nrow = arr.rows(), ncol = arr.cols();
	Eigen::ArrayXd col = Eigen::ArrayXd::LinSpaced(nrow, 1, nrow);
	Eigen::ArrayXXd x(nrow, ncol);
	
	x.colwise() = col;

	Eigen::ArrayXd 
		x_sum_y_sum = (x * arr).colwise().sum(), 
		x_sum = x.colwise().sum(), 
		y_sum = arr.colwise().sum(), 
		x_square_sum = x.matrix().colwise().squaredNorm(), 
		x_sum_square = x_sum * x_sum, 
		xy_sum_square = x_sum * y_sum;

	return (xy_sum_square - nrow * x_sum_y_sum) / (x_sum_square - nrow * x_square_sum);
};

void ScicnvWorker::compute_gradient() {

	int half_width = ceil(this->resolution_ / 2.);

	const int n_chromosome = this->chromosome_location_.size();

#pragma omp parallel for
	for (int i = 0; i < n_chromosome; ++i) {

		const auto& [chr, start, n] = this->chromosome_location_[i];
		const int end = start + n;

		for (int j = start; j < end; ++j) {

			int gene_start = j - half_width > start ? j - half_width : start;
			int gene_end = j + half_width < end ? j + half_width : end;

			Eigen::ArrayXd gradient = slope(this->U_(Eigen::seqN(gene_start, gene_end - gene_start), Eigen::all));
			this->V_.row(j) = gradient;
		}
	};
};

void ScicnvWorker::compute_cnv() {
	this->W_ *= this->V_;

	for (int j = 0; j < this->n_cell_; ++j) {
		for (int i = 0; i < this->n_gene_; ++i) {
			double val = this->W_(i, j);
			if (val < 0 || std::abs(val) < this->threshold_) {
				this->W_(i, j) = 0;
			}
		}
	}
	
	this->W_ = this->W_.sqrt();

	for (int j = 0; j < this->n_cell_; ++j) {
		for (int i = 0; i < this->n_gene_; ++i) {
			if (this->V_(i, j) < 0) {
				this->W_(i, j) = -this->W_(i, j);
			}
		}
	}
	
};

void ScicnvWorker::run() {

	if (!this->filter_data()) {
		G_TASK_END;
	}

	if (!this->sort_expression_by_chromosome()) {
		G_TASK_END;
	}

	this->compute_reference();

	this->compute_reference_moving_average(false);

	this->compute_moving_average(false);

	this->compute_W();

	this->compute_reference_moving_average(true);

	this->compute_moving_average(true);

	this->compute_U();

	this->compute_V();

	this->compute_gradient();

	this->compute_cnv();

	this->moving_average_.resize(0, 0);
	this->V_.resize(0, 0);
	this->U_.resize(0, 0);
	this->mat_.clear();

	CNV* cnv = new CNV(this->W_, this->metadata_location_, this->chromosome_location_, CNV::DataType::SciCnv);
	emit x_cnv_ready(cnv);
	G_TASK_END;
}
