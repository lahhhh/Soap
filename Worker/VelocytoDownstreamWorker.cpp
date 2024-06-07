#include "VelocytoDownstreamWorker.h"

#include "lm.h"

void VelocytoDownstreamWorker::run() {

	if (this->gene_relative_velocity_estimates()) {
		emit x_estimate_ready(this->estimate_);
	}

	G_TASK_END;
}

void VelocytoDownstreamWorker::get_conv_mat() {

	Eigen::MatrixXd emat = this->velocyto_base_->get_spliced()->mat_.toDense().cast<double>();

	Eigen::MatrixXd nmat = this->velocyto_base_->get_unspliced()->mat_.toDense().cast<double>();

	Eigen::ArrayXd emat_size = emat.colwise().sum();
	Eigen::ArrayXd nmat_size = nmat.colwise().sum();

	Eigen::ArrayXd emat_cs = emat_size / this->scale_factor_;
	Eigen::ArrayXd nmat_cs = nmat_size / this->scale_factor_;

	this->emat_cs_ = emat_cs;	

	Eigen::MatrixXd emat_norm = emat.array().rowwise() / emat_cs.transpose();

	Eigen::MatrixXd emat_log_norm = log(emat_norm.array() + this->pseudo_count_);
	emat_norm.resize(0, 0);

	constexpr int knn_maxl = 100;
	auto knn = _Cs balanced_knn_mt(emat_log_norm, this->n_cell_, knn_maxl * this->n_cell_);
	emat_log_norm.resize(0, 0);

	const int n_gene = emat.rows(), n_cell = emat.cols();

	Eigen::MatrixXi cell_knn(n_cell, this->n_cell_ + 1);

	cell_knn << knn, Eigen::ArrayXi::LinSpaced(n_cell, 0, n_cell - 1);
	knn.resize(0, 0);

	Eigen::MatrixXd conv_emat(n_gene, n_cell);
	Eigen::MatrixXd conv_nmat(n_gene, n_cell);
	Eigen::ArrayXd conv_emat_cs(n_cell);
	Eigen::ArrayXd conv_nmat_cs(n_cell);

#pragma omp parallel for
	for (int i = 0; i < n_cell; ++i) {
		Eigen::ArrayXi cell_index = cell_knn.row(i);
		cell_index = _Cs unique(cell_index);

		Eigen::ArrayXd col_eexp = Eigen::ArrayXd::Zero(n_gene);
		Eigen::ArrayXd col_nexp = Eigen::ArrayXd::Zero(n_gene);

		double emat_cell_cs{ 0.0 }, nmat_cell_cs{ 0.0 };

		for (auto index : cell_index) {
			col_eexp += emat.col(index).array();
			col_nexp += nmat.col(index).array();
			emat_cell_cs += emat_cs(index);
			nmat_cell_cs += nmat_cs(index);
		}

		conv_emat.col(i) = col_eexp;
		conv_nmat.col(i) = col_nexp;
		conv_emat_cs(i) = emat_cell_cs;
		conv_nmat_cs(i) = nmat_cell_cs;
	};

	cell_knn.resize(0, 0);

	Eigen::MatrixXd& conv_emat_norm = this->conv_emat_norm_1;
	Eigen::MatrixXd& conv_nmat_norm = this->conv_nmat_norm_1;

	conv_emat_norm = conv_emat.array().rowwise() / conv_emat_cs.transpose();
	conv_nmat_norm = conv_nmat.array().rowwise() / conv_nmat_cs.transpose();

	conv_emat.resize(0, 0);
	conv_nmat.resize(0, 0);	

	if (this->n_gene_ > 1) {

		this->conv_emat_norm_2 = conv_emat_norm;

		this->conv_nmat_norm_2 = conv_nmat_norm;


		knn = _Cs balanced_knn_mt(log(conv_emat_norm.array() + this->pseudo_count_).transpose(), this->n_gene_, this->n_gene_ * 1.2e3);
		Eigen::MatrixXd gene_knn = _Cs create_matrix_from_knn_index(knn).cast<double>().toDense();

		for (int i = 0; i < n_gene; ++i) {
			gene_knn(i, i) = 1;
		}

		Eigen::ArrayXd gt = conv_emat_norm.rowwise().sum();

		for (int i = 0; i < n_gene; ++i) {
			gene_knn.col(i).array() *= (_Cs median(_Cs sliced(gt, gene_knn.col(i).array() > 0)) / gt).cwiseMin(1.0).array();
		}

		gene_knn.transposeInPlace();

		conv_emat_norm = gene_knn * conv_emat_norm;
		conv_nmat_norm = gene_knn * conv_nmat_norm;
	}

};

bool VelocytoDownstreamWorker::fit_linear() {
	Eigen::MatrixXd& ko = this->ko_;

	Eigen::MatrixXd& conv_emat_norm = this->conv_emat_norm_1;
	Eigen::MatrixXd& conv_nmat_norm = this->conv_nmat_norm_1;

	const int n_gene = conv_emat_norm.rows();

	ko.resize(n_gene, 3);

#pragma omp parallel for
	for (int i = 0; i < n_gene; ++i) {
		Eigen::ArrayXd e_array = conv_emat_norm.row(i);
		Eigen::ArrayXd n_array = conv_nmat_norm.row(i);

		ko(i, 2) = _Cs correlation_spearman(e_array, n_array);

		auto sorted_e = _Cs sorted(e_array);
		// fit.quantile = 0.02
		double p2 = _Cs linear_percentile_no_sort(sorted_e, 2);
		double p98 = _Cs linear_percentile_no_sort(sorted_e, 98);

		auto filter = (e_array <= p2) + (e_array >= p98);
		e_array = _Cs sliced(e_array, filter);
		n_array = _Cs sliced(n_array, filter);

		auto fit_res = lm(n_array, e_array);
		ko(i, 0) = fit_res.coefficients[0];
		ko(i, 1) = fit_res.coefficients[1];
	};

	Eigen::ArrayX<bool> nan_filter = Eigen::ArrayX<bool>::Constant(n_gene, true);

	int n_gene_remain = n_gene;

	for (int i = 0; i < n_gene; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (std::isnan(ko(i, j))) {
				nan_filter[i] = false;
				break;
			}
		}
	}

	int pass = nan_filter.count();

	if (pass == 0) {
		G_TASK_WARN("No gene passed NaN filter");
		return false;
	}
	else if (pass < n_gene_remain) {
		G_TASK_LOG(QString::number(n_gene_remain - pass) + " gene filtered out of " + QString::number(n_gene_remain) + " genes due to nan filter.");

		n_gene_remain = pass;
		ko = _Cs row_sliced(ko, nan_filter);
	}


	constexpr double min_nmat_emat_correlation = 0.05;
	Eigen::ArrayX<bool> correlation_filter = ko.col(2).array() > min_nmat_emat_correlation;

	pass = correlation_filter.count();

	if (pass == 0) {
		G_TASK_WARN("No gene passed correlation filter");
		return false;
	}
	else if (pass < n_gene_remain) {
		G_TASK_LOG(QString::number(n_gene_remain - pass) + " gene filtered out of " + QString::number(n_gene_remain) + " genes due to correlation filter.");
		n_gene_remain = pass;
		ko = _Cs row_sliced(ko, correlation_filter);
	}


	constexpr double min_nmat_emat_slope = 0.05;
	Eigen::ArrayX<bool> slope_filter = ko.col(1).array() > min_nmat_emat_slope;

	pass = slope_filter.count();

	if (pass == 0) {
		G_TASK_WARN("No gene passed slope filter");
		return false;
	}
	else if (pass < n_gene_remain) {
		G_TASK_LOG(QString::number(n_gene_remain - pass) + " gene filtered out of " + QString::number(n_gene_remain) + " genes due to slope filter.");
		n_gene_remain = pass;
		ko = _Cs row_sliced(ko, slope_filter);
	}

	auto gene_index = _Cs seq_n(0, n_gene);
	gene_index = _Cs sliced(gene_index, nan_filter);
	gene_index = _Cs sliced(gene_index, correlation_filter);
	gene_index = _Cs sliced(gene_index, slope_filter);

	this->ko_index_ = gene_index;
};

void VelocytoDownstreamWorker::calculate_velocity_shift() {

	Eigen::MatrixXd& ko = this->ko_;
	Eigen::MatrixXd& conv_emat_norm_1 = this->conv_emat_norm_1;
	Eigen::MatrixXd& conv_nmat_norm_1 = this->conv_nmat_norm_1;

	if (this->n_gene_ == 1) {

		Eigen::ArrayXd egt = ( - ko.col(1).array()).exp();
		Eigen::MatrixXd y = (conv_nmat_norm_1(this->ko_index_, Eigen::all).array().colwise() - this->ko_.col(0).array()).cwiseMax(0.0);
		this->estimate_->deltaE_ = conv_emat_norm_1(this->ko_index_, Eigen::all).array().colwise() * egt + (y.array().colwise() * (1 - egt)).colwise() / this->ko_.col(1).array() - conv_emat_norm_1(this->ko_index_, Eigen::all).array();
		this->gamma_index_ = this->ko_index_;
		this->estimate_->gene_names_ = _Cs reordered(this->velocyto_base_->get_spliced()->rownames_, this->gamma_index_);
		return;
	}
		

	Eigen::MatrixXd n_prediction = ((conv_emat_norm_1(this->ko_index_, Eigen::all).array().colwise() * ko.col(1).array()).colwise() + ko.col(0).array()).cwiseMax(0.0);

	conv_emat_norm_1.resize(0, 0);

	Eigen::MatrixXd mval = log2(conv_nmat_norm_1(this->ko_index_, Eigen::all).array() + this->pseudo_count_) - log2(n_prediction.array() + this->pseudo_count_);

	conv_nmat_norm_1.resize(0, 0);

	Eigen::MatrixXd& conv_emat_norm = this->conv_emat_norm_2;
	Eigen::MatrixXd& conv_nmat_norm = this->conv_nmat_norm_2;

	Eigen::MatrixXd am = conv_nmat_norm(this->ko_index_, Eigen::all).array().colwise() - ko.col(0).array();

	am = am.array().cwiseMax(0.0);

	Eigen::MatrixXd fm = log2(am.array()) - mval.array() - log2(conv_emat_norm(this->ko_index_, Eigen::all).array());

	const int n_gene_remain = this->ko_index_.size(), n_cell = conv_emat_norm.cols();
	Eigen::ArrayXd row_sum_wm = Eigen::ArrayXd::Zero(n_gene_remain);

	for (int j = 0; j < n_cell; ++j) {
		for (int i = 0; i < n_gene_remain; ++i) {
			bool finite = std::isfinite(fm(i, j));
			if (finite) {
				++row_sum_wm[i];
			}
			else {
				fm(i, j) = 0;
			}
		}
	}

	Eigen::ArrayX<bool> gamma_filter = row_sum_wm > 0;

	auto gene_index2 = _Cs sliced(this->ko_index_, gamma_filter);
	this->gamma_index_ = gene_index2;
	this->estimate_->gene_names_ = _Cs reordered(this->velocyto_base_->get_spliced()->rownames_, this->gamma_index_);
	auto gene_index3 = _Cs which(gamma_filter);

	Eigen::ArrayXd row_sum_fm = fm.array().rowwise().sum();

	Eigen::ArrayXd gammaA = pow(2, _Cs sliced(row_sum_fm, gamma_filter) / _Cs sliced(row_sum_wm, gamma_filter));

	Eigen::ArrayXXd r = pow(2.0, mval(gene_index3, Eigen::all).array());

	// get projected delta from log2ratio delta = 1.0
	Eigen::MatrixXd deltaE = (conv_emat_norm(gene_index2, Eigen::all).array() + 1e-4) * \
		((1 - r) * exp(-gammaA) + r) - conv_emat_norm(gene_index2, Eigen::all).array();

	conv_emat_norm.resize(0, 0);
	conv_nmat_norm.resize(0, 0);

	this->estimate_->deltaE_ = deltaE;
};

void VelocytoDownstreamWorker::calculate_extrapolated_cell_state() {

	this->estimate_->current_ = this->velocyto_base_->get_spliced()->mat_.toDense().cast<double>();

	this->estimate_->current_ = this->estimate_->current_(this->gamma_index_, Eigen::all).eval();

	auto& emat = this->estimate_->current_;

	double delta = 1.0;
	this->estimate_->projected_ = emat.array() + this->estimate_->deltaE_.array().rowwise() * this->emat_cs_.transpose() * delta;
	this->estimate_->projected_ = this->estimate_->projected_.array().cwiseMax(0.0);


	Eigen::ArrayXd new_cell_size = this->emat_cs_ + (this->estimate_->projected_.array() - emat.array()).colwise().sum() / this->scale_factor_;
	this->estimate_->projected_.array().rowwise() / new_cell_size.transpose();

	this->estimate_->current_.array().rowwise() /= this->emat_cs_.transpose();
};

bool VelocytoDownstreamWorker::gene_relative_velocity_estimates() {

	this->estimate_ = new VelocityEstimate();

	this->get_conv_mat();

	if (!this->fit_linear()) {
		delete this->estimate_;
		return false;
	}
	this->calculate_velocity_shift();

	this->calculate_extrapolated_cell_state();

	return true;
};
