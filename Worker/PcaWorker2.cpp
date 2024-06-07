#include "PcaWorker2.h"

#include "Custom.h"

/*
Modified from R package Seurat : FindVariableFeatures.
*/
Eigen::ArrayX<bool> find_variable_features(
	const Eigen::MatrixXd& mat,
	int n_variable_feature)
{

	const int ncol = mat.cols();
	const int nrow = mat.rows();
	double clipmax = std::sqrt(ncol);
	Eigen::ArrayXd row_mean = mat.rowwise().mean();
	Eigen::ArrayXd row_var = _Cs row_var_mt(mat);

	Eigen::ArrayX<bool> not_const = row_var > 0;
	const int not_const_row = not_const.count();
	if (not_const_row < n_variable_feature) {
		return {};
	}
	auto not_const_index = _Cs which(not_const);

	Eigen::ArrayXd not_const_row_means = _Cs sliced(row_mean, not_const);
	Eigen::ArrayXd not_const_row_var = _Cs sliced(row_var, not_const);

	Eigen::ArrayXd fitted = _Cs loess_mt(log10(not_const_row_var), log10(not_const_row_means), 2, 0.3, 50);

	Eigen::ArrayXd expected_var = Eigen::ArrayXd::Zero(nrow);
	expected_var(not_const_index) = pow(10, fitted);

	Eigen::ArrayXd expected_sd = sqrt(expected_var);
	Eigen::ArrayXd standardized_row_var = Eigen::ArrayXd::Zero(nrow);
	Eigen::ArrayXi n_zero = Eigen::ArrayXi::Constant(nrow, ncol);

	for (int i = 0; i < nrow; ++i) {
		if (expected_sd[i] == 0)continue;

		for (int j = 0; j < ncol; ++j) {
			double s = std::min(clipmax, std::abs((mat(i, j) - row_mean[i]) / expected_sd[i]));
			standardized_row_var[i] += s * s;
		}
	}

	standardized_row_var /= (ncol - 1);
	Eigen::ArrayXi standardized_row_var_sorted_index = _Cs order(standardized_row_var, true);
	Eigen::ArrayX<bool> res = Eigen::ArrayX<bool>::Constant(nrow, false);

	for (int i = 0; i < n_variable_feature && i < standardized_row_var_sorted_index.size(); ++i) {
		int index = standardized_row_var_sorted_index[i];
		res[index] = true;
	}
	return res;
}

void PcaWorker2::run() {

	if (this->use_variable_features_) {

		auto features = find_variable_features(this->mat_, this->n_variable_feature_);
		if (features.size() == 0) {
			G_TASK_WARN("No enough variable features meeting requirement.");
			G_TASK_END;
		}

		this->mat_ = _Cs row_sliced(this->mat_, features);
	}

	this->mat_.transposeInPlace();

	_Cs scale_in_place(this->mat_);

    auto vars = _Cs col_var_mt(this->mat_);
    double total_var = vars.sum();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(this->mat_, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd singular_values = svd.singularValues();
    Eigen::MatrixXd V = svd.matrixV();
	Eigen::MatrixXd emb = this->mat_ * V;

    Eigen::ArrayXd sdev = singular_values.array() / sqrt(this->mat_.rows() - 1.0);

    Eigen::ArrayXd variance_proportion = sdev.square() / total_var;

    emit x_pca_ready(emb, _Cs cast<QVector>(sdev), _Cs cast<QVector>(variance_proportion));

    G_TASK_END;
};