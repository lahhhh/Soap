#include "PcaWorker.h"
#include "TruncatedSvd.h"

#include "Custom.h"

/*
Modified from R package Seurat : FindVariableFeatures.
*/
static Eigen::ArrayX<bool> find_variable_features(
	const Eigen::SparseMatrix<double>& mat, 
	int n_variable_feature) 
{
	
	const int ncol = mat.cols();
	const int nrow = mat.rows();
	double clipmax = std::sqrt(ncol);
	Eigen::ArrayXd row_mean = custom::row_mean(mat);
	Eigen::ArrayXd row_var = custom::row_var(mat);

	Eigen::ArrayX<bool> not_const = row_var > 0;
	const int not_const_row = not_const.count();
	if (not_const_row < n_variable_feature) {
		return {};
	}
	auto not_const_index = custom::which(not_const);

	Eigen::ArrayXd not_const_row_means = custom::sliced(row_mean, not_const);
	Eigen::ArrayXd not_const_row_var = custom::sliced(row_var, not_const);

	Eigen::ArrayXd fitted = custom::loess_mt(log10(not_const_row_var), log10(not_const_row_means), 2, 0.3, 50);

	Eigen::ArrayXd expected_var = Eigen::ArrayXd::Zero(nrow);
	expected_var(not_const_index) = pow(10, fitted);

	Eigen::ArrayXd expected_sd = sqrt(expected_var);
	Eigen::ArrayXd standardized_row_var = Eigen::ArrayXd::Zero(nrow);
	Eigen::ArrayXi n_zero = Eigen::ArrayXi::Constant(nrow, ncol);

	for (int i = 0; i < ncol; ++i) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
		{
			int row = it.row();
			--n_zero[row];
			double s = std::min(clipmax, std::abs((it.value() - row_mean[row]) / expected_sd[row]));
			standardized_row_var[row] += s * s;
		}
	}

	for (int i = 0; i < nrow; ++i) {
		if (expected_sd[i] == 0)continue;
		standardized_row_var[i] += std::pow(row_mean[i] / expected_sd[i], 2) * n_zero[i];
	}
	standardized_row_var /= (ncol - 1);
	Eigen::ArrayXi standardized_row_var_sorted_index = custom::order(standardized_row_var, true);
	Eigen::ArrayX<bool> res = Eigen::ArrayX<bool>::Constant(nrow, false);

	for (int i = 0; i < n_variable_feature && i < standardized_row_var_sorted_index.size(); ++i) {
		int index = standardized_row_var_sorted_index[i];
		res[index] = true;
	}
	return res;
}


/*
Modified from R package Seurat : FindVariableFeatures.
*/
static Eigen::ArrayX<bool> find_variable_features(
	const Eigen::MatrixXd& mat,
	int n_variable_feature)
{

	const int ncol = mat.cols();
	const int nrow = mat.rows();
	double clipmax = std::sqrt(ncol);
	Eigen::ArrayXd row_mean = mat.rowwise().mean();
	Eigen::ArrayXd row_var = (mat.array().colwise() - row_mean).square().rowwise().sum() / (ncol - 1);

	Eigen::ArrayX<bool> not_const = row_var > 0;
	const int not_const_row = not_const.count();
	if (not_const_row < n_variable_feature) {
		return {};
	}
	auto not_const_index = custom::which(not_const);

	Eigen::ArrayXd not_const_row_means = custom::sliced(row_mean, not_const);
	Eigen::ArrayXd not_const_row_var = custom::sliced(row_var, not_const);

	Eigen::ArrayXd fitted = custom::loess_mt(log10(not_const_row_var), log10(not_const_row_means), 2, 0.3, 50);

	Eigen::ArrayXd expected_var = Eigen::ArrayXd::Zero(nrow);
	expected_var(not_const_index) = pow(10, fitted);

	Eigen::ArrayXd expected_sd = sqrt(expected_var);
	Eigen::ArrayXd standardized_row_var = Eigen::ArrayXd::Zero(nrow);
	Eigen::ArrayXi n_zero = Eigen::ArrayXi::Constant(nrow, ncol);

	for (int i = 0; i < nrow; ++i) {

		if (expected_sd[i] == 0.0) {
			continue;
		}

		standardized_row_var[i] = ((mat.array().row(i) - row_mean[i]).abs().cwiseMin(clipmax) / expected_sd[i]).square().sum() / (ncol - 1);
	}

	Eigen::ArrayXi standardized_row_var_sorted_index = custom::order(standardized_row_var, true);
	Eigen::ArrayX<bool> res = Eigen::ArrayX<bool>::Constant(nrow, false);

	for (int i = 0; i < n_variable_feature && i < standardized_row_var_sorted_index.size(); ++i) {
		int index = standardized_row_var_sorted_index[i];
		res[index] = true;
	}
	return res;
}

Eigen::MatrixXd pca_infercnv(
	const Eigen::MatrixXd& mat,
	int n_variable_feature,
	int nu,
	int random_state,
	double tol,
	int max_iter)
{

	auto vars = find_variable_features(mat, n_variable_feature);
	if (vars.size() == 0) {
		return {};
	}

	Eigen::MatrixXd m = custom::row_sliced(mat, vars).transpose();
	custom::scale_in_place(m);
	auto [U, S, V] = tsvd(&m, nu, random_state, tol, max_iter);
	Eigen::MatrixXd emb = U * S.asDiagonal();
	return emb;
};

bool PcaWorker::work() {

	Eigen::SparseMatrix<double> mat = custom::normalize(*this->mat_, 10000.0);

	int n_feature = this->feature_proportion_ <= 1.0 ? mat.rows() * this->feature_proportion_ : this->feature_proportion_;

	auto vars = find_variable_features(mat, n_feature);
	if (vars.size() == 0) {
		G_TASK_WARN("No enough variable features meeting requirement.");
		return false;
	}

	mat = custom::row_sliced(mat, vars);
	Eigen::MatrixXd scaled_matrix = custom::row_scale_mt(mat);
	mat.resize(0, 0);
	scaled_matrix.transposeInPlace();
	auto [U, S, V] = tsvd(&scaled_matrix, this->n_mat_u_, this->random_state_, this->tol_, this->maximum_iteration_);
	Eigen::ArrayXd sdev = S.array() / sqrt(scaled_matrix.rows() - 1.);
	scaled_matrix.resize(0, 0);

	this->res_ = U * S.asDiagonal();
	this->sdev_ = custom::cast<QVector>(sdev);

	return true;
};

void PcaWorker::run() {

	G_TASK_LOG("Start PCA...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_pca_ready(this->res_, this->sdev_);

	G_TASK_LOG("PCA finished.");

	G_TASK_END;
}
