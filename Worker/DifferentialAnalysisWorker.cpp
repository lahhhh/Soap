#include "DifferentialAnalysisWorker.h"

#include "Custom.h"
#include "WilcoxTest.h"

void DifferentialAnalysisWorker::dense(DifferentialAnalysis::DataType dtype) {

	if (this->comparison_[1] != "ALL") {
		Eigen::ArrayX<bool> group_1_filter = custom::equal(this->metadata_, this->comparison_[1]), group_2_filter;

		if (this->comparison_[2] != "REST") {
			group_2_filter = custom::equal(this->metadata_, this->comparison_[2]);
		}
		else {
			group_2_filter = !group_1_filter;
		}

		int feature_number = this->feature_names_.size();
		Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(feature_number, false);
		Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd log2_fold_change = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd percentage_of_val_in_group_1 = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd percentage_of_val_in_group_2 = Eigen::ArrayXd::Zero(feature_number);
		double n_cell_1 = group_1_filter.count(), n_cell_2 = group_2_filter.count();

	#pragma omp parallel for
		for (int i = 0; i < feature_number; ++i) {
			Eigen::ArrayXd feature_val = this->dense_data_.row(i);
			Eigen::ArrayXd group_1_val = custom::sliced(feature_val, group_1_filter);
			Eigen::ArrayXd group_2_val = custom::sliced(feature_val, group_2_filter);

			double percentage_1 = group_1_val.count() / n_cell_1;
			double percentage_2 = group_2_val.count() / n_cell_2;

			if (percentage_1 < this->minimum_percentage_ && percentage_2 < this->minimum_percentage_) {
				continue;
			}
			percentage_of_val_in_group_1[i] = percentage_1;
			percentage_of_val_in_group_2[i] = percentage_2;
			double denominator = group_2_val.mean(), numerator = group_1_val.mean();
			log2_fold_change[i] = log2((numerator + 1) / (denominator + 1));
			computed[i] = true;
			p_adjusted[i] = wilcox_test(group_1_val, group_2_val);
		}

		this->feature_names_ = custom::sliced(this->feature_names_, computed);
		log2_fold_change = custom::sliced(log2_fold_change, computed);
		p_adjusted = custom::sliced(p_adjusted, computed);
		feature_number = this->feature_names_.size();
		percentage_of_val_in_group_1 = custom::sliced(percentage_of_val_in_group_1, computed);
		percentage_of_val_in_group_2 = custom::sliced(percentage_of_val_in_group_2, computed);
		p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);


		this->res_.data_type_ = dtype;
		this->res_.mat_.set_rownames(this->feature_names_);
		this->res_.mat_.update(METADATA_DE_FEATURE_NAME, this->feature_names_);
		this->res_.mat_.update(METADATA_DE_LOG2_FOLD_CHANGE, custom::cast<QVector>(log2_fold_change));
		this->res_.mat_.update(METADATA_DE_ADJUSTED_P_VALUE, custom::cast<QVector>(p_adjusted));
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_1, custom::cast<QVector>(percentage_of_val_in_group_1));
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_2, custom::cast<QVector>(percentage_of_val_in_group_2));
		this->res_.mat_.update(METADATA_DE_COMPARISON_1, QStringList(feature_number, this->comparison_[1]), CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_COMPARISON_2, QStringList(feature_number, this->comparison_[2]), CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_METADATA_NAME, QStringList(feature_number, this->metadata_name_), CustomMatrix::DataType::QStringFactor);
		
	}
	else {

		QStringList factors = custom::unique(this->metadata_);
		QStringList all_feature_names;
		QVector<double> all_log2_fold_change, all_p_adjusted, all_percentage_1, all_percentage_2;
		QStringList comparison1, comparison2;

		for (auto factor : factors) {
			Eigen::ArrayX<bool> group_1_filter = custom::equal(this->metadata_, factor), group_2_filter;
			if (this->comparison_[2] != "REST") {
				group_2_filter = custom::equal(this->metadata_, this->comparison_[2]);
			}
			else {
				group_2_filter = !group_1_filter;
			}
			int feature_number = this->feature_names_.size();
			Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(feature_number, false);
			Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd log2_fold_change = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd percentage_of_val_in_group_1 = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd percentage_of_val_in_group_2 = Eigen::ArrayXd::Zero(feature_number);
			double n_cell_1 = group_1_filter.count(), n_cell_2 = group_2_filter.count();

		#pragma omp parallel for
			for (int i = 0; i < feature_number; ++i) {
				Eigen::ArrayXd feature_val = this->dense_data_.row(i);
				Eigen::ArrayXd group_1_val = custom::sliced(feature_val, group_1_filter);
				Eigen::ArrayXd group_2_val = custom::sliced(feature_val, group_2_filter);
				double percentage_1 = group_1_val.count() / n_cell_1;
				double percentage_2 = group_2_val.count() / n_cell_2;
				if (percentage_1 < this->minimum_percentage_ && percentage_2 < this->minimum_percentage_) {
					continue;
				}
				percentage_of_val_in_group_1[i] = percentage_1;
				percentage_of_val_in_group_2[i] = percentage_2;
				double denominator = group_2_val.mean(), numerator = group_1_val.mean();

				log2_fold_change[i] = log2((numerator + 1) / (denominator + 1));
				computed[i] = true;
				p_adjusted[i] = wilcox_test(group_1_val, group_2_val);
			};
			QStringList feature_names = custom::sliced(this->feature_names_, computed);
			log2_fold_change = custom::sliced(log2_fold_change, computed);
			p_adjusted = custom::sliced(p_adjusted, computed);
			percentage_of_val_in_group_1 = custom::sliced(percentage_of_val_in_group_1, computed);
			percentage_of_val_in_group_2 = custom::sliced(percentage_of_val_in_group_2, computed);
			p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);
			all_feature_names.append(feature_names);
			all_log2_fold_change.append(custom::cast<QVector>(log2_fold_change));
			all_p_adjusted.append(custom::cast<QVector>(p_adjusted));
			all_percentage_1.append(custom::cast<QVector>(percentage_of_val_in_group_1));
			all_percentage_2.append(custom::cast<QVector>(percentage_of_val_in_group_2));
			comparison1.append(QStringList(feature_names.size(), factor));
			comparison2.append(QStringList(feature_names.size(), this->comparison_[2]));
		}


		this->res_.data_type_ = dtype;
		this->res_.mat_.set_nrow(all_feature_names.size());
		this->res_.mat_.update(METADATA_DE_FEATURE_NAME, all_feature_names);
		this->res_.mat_.update(METADATA_DE_LOG2_FOLD_CHANGE, all_log2_fold_change);
		this->res_.mat_.update(METADATA_DE_ADJUSTED_P_VALUE, all_p_adjusted);
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_1, all_percentage_1);
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_2, all_percentage_2);
		this->res_.mat_.update(METADATA_DE_COMPARISON_1, comparison1, CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_COMPARISON_2, comparison2, CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_METADATA_NAME, QStringList(all_feature_names.size(), this->metadata_name_), CustomMatrix::DataType::QStringFactor);
		
	}
};

void DifferentialAnalysisWorker::sparse(DifferentialAnalysis::DataType dtype) {

	if (this->comparison_[1] != "ALL") {
		Eigen::ArrayX<bool> group_1_filter = custom::equal(this->metadata_, this->comparison_[1]), group_2_filter;

		if (this->comparison_[2] != "REST") {
			group_2_filter = custom::equal(this->metadata_, this->comparison_[2]);
		}
		else {
			group_2_filter = !group_1_filter;
		}
		int feature_number = this->feature_names_.size();
		Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(feature_number, false);
		Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd log2_fold_change = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd percentage_of_val_in_group_1 = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd percentage_of_val_in_group_2 = Eigen::ArrayXd::Zero(feature_number);
		double n_cell_1 = group_1_filter.count(), n_cell_2 = group_2_filter.count();

	#pragma omp parallel for
		for (int i = 0; i < feature_number; ++i) {
			Eigen::ArrayXd feature_val = this->sparse_data_.col(i);
			Eigen::ArrayXd group_1_val = custom::sliced(feature_val, group_1_filter);
			Eigen::ArrayXd group_2_val = custom::sliced(feature_val, group_2_filter);

			double percentage_1 = group_1_val.count() / n_cell_1;
			double percentage_2 = group_2_val.count() / n_cell_2;

			if (percentage_1 < this->minimum_percentage_ && percentage_2 < this->minimum_percentage_) {
				continue;
			}
			percentage_of_val_in_group_1[i] = percentage_1;
			percentage_of_val_in_group_2[i] = percentage_2;
			double denominator = group_2_val.mean(), numerator = group_1_val.mean();
			log2_fold_change[i] = log2((numerator + 1) / (denominator + 1));
			computed[i] = true;
			p_adjusted[i] = wilcox_test(group_1_val, group_2_val);
		}

		this->feature_names_ = custom::sliced(this->feature_names_, computed);
		log2_fold_change = custom::sliced(log2_fold_change, computed);
		p_adjusted = custom::sliced(p_adjusted, computed);
		feature_number = this->feature_names_.size();
		percentage_of_val_in_group_1 = custom::sliced(percentage_of_val_in_group_1, computed);
		percentage_of_val_in_group_2 = custom::sliced(percentage_of_val_in_group_2, computed);
		p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);


		this->res_.data_type_ = dtype;
		this->res_.mat_.set_rownames(this->feature_names_);
		this->res_.mat_.update(METADATA_DE_FEATURE_NAME, this->feature_names_);
		this->res_.mat_.update(METADATA_DE_LOG2_FOLD_CHANGE, custom::cast<QVector>(log2_fold_change));
		this->res_.mat_.update(METADATA_DE_ADJUSTED_P_VALUE, custom::cast<QVector>(p_adjusted));
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_1, custom::cast<QVector>(percentage_of_val_in_group_1));
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_2, custom::cast<QVector>(percentage_of_val_in_group_2));
		this->res_.mat_.update(METADATA_DE_COMPARISON_1, QStringList(feature_number, this->comparison_[1]), CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_COMPARISON_2, QStringList(feature_number, this->comparison_[2]), CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_METADATA_NAME, QStringList(feature_number, this->metadata_name_), CustomMatrix::DataType::QStringFactor);
		
	}
	else {
		QStringList factors = custom::unique(this->metadata_);
		QStringList all_feature_names;
		QVector<double> all_log2_fold_change, all_p_adjusted, all_percentage_1, all_percentage_2;
		QStringList comparison1, comparison2;
		for (auto factor : factors) {
			Eigen::ArrayX<bool> group_1_filter = custom::equal(this->metadata_, factor), group_2_filter;
			if (this->comparison_[2] != "REST") {
				group_2_filter = custom::equal(this->metadata_, this->comparison_[2]);
			}
			else {
				group_2_filter = !group_1_filter;
			}
			int feature_number = this->feature_names_.size();
			Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(feature_number, false);
			Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd log2_fold_change = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd percentage_of_val_in_group_1 = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd percentage_of_val_in_group_2 = Eigen::ArrayXd::Zero(feature_number);
			double n_cell_1 = group_1_filter.count(), n_cell_2 = group_2_filter.count();

		#pragma omp parallel for
			for (int i = 0; i < feature_number; ++i) {
				Eigen::ArrayXd feature_val = this->sparse_data_.col(i);
				Eigen::ArrayXd group_1_val = custom::sliced(feature_val, group_1_filter);
				Eigen::ArrayXd group_2_val = custom::sliced(feature_val, group_2_filter);
				double percentage_1 = group_1_val.count() / n_cell_1;
				double percentage_2 = group_2_val.count() / n_cell_2;

				if (percentage_1 < this->minimum_percentage_ && percentage_2 < this->minimum_percentage_) {
					continue;
				}
				percentage_of_val_in_group_1[i] = percentage_1;
				percentage_of_val_in_group_2[i] = percentage_2;
				double denominator = group_2_val.mean(), numerator = group_1_val.mean();
				log2_fold_change[i] = log2((numerator + 1) / (denominator + 1));
				computed[i] = true;
				p_adjusted[i] = wilcox_test(group_1_val, group_2_val);
			};

			QStringList feature_names = custom::sliced(this->feature_names_, computed);
			log2_fold_change = custom::sliced(log2_fold_change, computed);
			p_adjusted = custom::sliced(p_adjusted, computed);
			percentage_of_val_in_group_1 = custom::sliced(percentage_of_val_in_group_1, computed);
			percentage_of_val_in_group_2 = custom::sliced(percentage_of_val_in_group_2, computed);
			p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);
			all_feature_names.append(feature_names);
			all_log2_fold_change.append(custom::cast<QVector>(log2_fold_change));
			all_p_adjusted.append(custom::cast<QVector>(p_adjusted));
			all_percentage_1.append(custom::cast<QVector>(percentage_of_val_in_group_1));
			all_percentage_2.append(custom::cast<QVector>(percentage_of_val_in_group_2));
			comparison1.append(QStringList(feature_names.size(), factor));
			comparison2.append(QStringList(feature_names.size(), this->comparison_[2]));
		}


		this->res_.data_type_ = dtype;
		this->res_.mat_.set_nrow(all_feature_names.size());
		this->res_.mat_.update(METADATA_DE_FEATURE_NAME, all_feature_names);
		this->res_.mat_.update(METADATA_DE_LOG2_FOLD_CHANGE, all_log2_fold_change);
		this->res_.mat_.update(METADATA_DE_ADJUSTED_P_VALUE, all_p_adjusted);
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_1, all_percentage_1);
		this->res_.mat_.update(METADATA_DE_PERCENTAGE_2, all_percentage_2);
		this->res_.mat_.update(METADATA_DE_COMPARISON_1, comparison1, CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_COMPARISON_2, comparison2, CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_METADATA_NAME, QStringList(all_feature_names.size(), this->metadata_name_), CustomMatrix::DataType::QStringFactor);
		
	}
};

void DifferentialAnalysisWorker::cicero_mode() {

	this->dense_data_ = this->cicero_->regulation_group_normalized_.mat_.array().exp() - 1.0;

	this->dense(DifferentialAnalysis::DataType::Cicero);
};

void DifferentialAnalysisWorker::chromvar_mode() {

	this->dense_data_ = this->chrom_var_->z_;
	if (this->comparison_[1] != "ALL") {
		Eigen::ArrayX<bool> group_1_filter = custom::equal(this->metadata_, this->comparison_[1]), group_2_filter;

		if (this->comparison_[2] != "REST") {
			group_2_filter = custom::equal(this->metadata_, this->comparison_[2]);
		}
		else {
			group_2_filter = !group_1_filter;
		}
		int feature_number = this->feature_names_.size();
		Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(feature_number);
		Eigen::ArrayXd z_diff = Eigen::ArrayXd::Zero(feature_number);
		double n_cell_1 = group_1_filter.count(), n_cell_2 = group_2_filter.count();

	#pragma omp parallel for
		for (int i = 0; i < feature_number; ++i) {
			Eigen::ArrayXd feature_val = this->dense_data_.row(i);
			Eigen::ArrayXd group_1_val = custom::sliced(feature_val, group_1_filter);
			Eigen::ArrayXd group_2_val = custom::sliced(feature_val, group_2_filter);

			double percentage_1 = group_1_val.count() / n_cell_1;
			double percentage_2 = group_2_val.count() / n_cell_2;
			
			double val2 = group_2_val.mean(), val1 = group_1_val.mean();
			p_adjusted[i] = wilcox_test(group_1_val, group_2_val);
			z_diff[i] = val1 - val2;
		}

		p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);

		this->res_.data_type_ = DifferentialAnalysis::DataType::ChromVAR;
		this->res_.mat_.set_rownames(this->feature_names_);
		this->res_.mat_.update(METADATA_DE_FEATURE_NAME, this->feature_names_);
		this->res_.mat_.update(METADATA_DE_Z_SCORE_DIFFERENCE, custom::cast<QVector>(z_diff));
		this->res_.mat_.update(METADATA_DE_ADJUSTED_P_VALUE, custom::cast<QVector>(p_adjusted));
		this->res_.mat_.update(METADATA_DE_COMPARISON_1, QStringList(feature_number, this->comparison_[1]), CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_COMPARISON_2, QStringList(feature_number, this->comparison_[2]), CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_METADATA_NAME, QStringList(feature_number, this->metadata_name_), CustomMatrix::DataType::QStringFactor);
		
	}
	else {
		QStringList factors = custom::unique(this->metadata_);
		QStringList all_feature_names;  
		QVector<double> all_z_diff, all_p_adjusted;
		QStringList comparison1, comparison2;
		for (auto&& factor : factors) {
			Eigen::ArrayX<bool> group_1_filter = custom::equal(this->metadata_, factor), group_2_filter;
			if (this->comparison_[2] != "REST") {
				group_2_filter = custom::equal(this->metadata_, this->comparison_[2]);
			}
			else {
				group_2_filter = !group_1_filter;
			}
			int feature_number = this->feature_names_.size();
			Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(feature_number, false);
			Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(feature_number);
			Eigen::ArrayXd z_diff = Eigen::ArrayXd::Zero(feature_number);
			double n_cell_1 = group_1_filter.count(), n_cell_2 = group_2_filter.count();

		#pragma omp parallel for
			for (int i = 0; i < feature_number; ++i) {
				Eigen::ArrayXd feature_val = this->dense_data_.row(i);
				Eigen::ArrayXd group_1_val = custom::sliced(feature_val, group_1_filter);
				Eigen::ArrayXd group_2_val = custom::sliced(feature_val, group_2_filter);

				double val2 = group_2_val.mean(), val1 = group_1_val.mean();
				p_adjusted[i] = wilcox_test(group_1_val, group_2_val);
				z_diff[i] = val1 - val2;
			};

			p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);
			all_feature_names.append(this->feature_names_);
			all_z_diff.append(custom::cast<QVector>(z_diff));
			all_p_adjusted.append(custom::cast<QVector>(p_adjusted));
			comparison1.append(QStringList(feature_number, factor));
			comparison2.append(QStringList(feature_number, this->comparison_[2]));
		}

		this->res_.data_type_ = DifferentialAnalysis::DataType::ChromVAR;
		this->res_.mat_.set_nrow(all_feature_names.size());
		this->res_.mat_.update(METADATA_DE_FEATURE_NAME, all_feature_names);
		this->res_.mat_.update(METADATA_DE_Z_SCORE_DIFFERENCE, all_z_diff);
		this->res_.mat_.update(METADATA_DE_ADJUSTED_P_VALUE, all_p_adjusted);
		this->res_.mat_.update(METADATA_DE_COMPARISON_1, comparison1, CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_COMPARISON_2, comparison2, CustomMatrix::DataType::QStringFactor);
		this->res_.mat_.update(METADATA_DE_METADATA_NAME, QStringList(all_feature_names.size(), this->metadata_name_), CustomMatrix::DataType::QStringFactor);
		
	}
};

void DifferentialAnalysisWorker::scatac_mode() {

	this->sparse_data_ = custom::normalize(this->single_cell_atac_->counts()->mat_, 10000.0).transpose();

	this->sparse(DifferentialAnalysis::DataType::Peak);
};

void DifferentialAnalysisWorker::bulkrna_mode() {

	this->dense_data_ = custom::normalize(this->bulk_rna_->counts()->mat_.cast<double>(), 10000.0);

	this->dense(DifferentialAnalysis::DataType::Gene);
};

void DifferentialAnalysisWorker::scrna_mode() {

	this->sparse_data_ = custom::normalize(this->single_cell_rna_->counts()->mat_, 10000.0).transpose();

	this->sparse(DifferentialAnalysis::DataType::Gene);
};

void DifferentialAnalysisWorker::scmultiome_mode() {

	this->sparse_data_ = custom::normalize(this->data_field_->counts()->mat_, 10000).transpose();

	auto field_type = this->data_field_->data_type_;

	if (field_type == DataField::DataType::Rna) {

		this->sparse(DifferentialAnalysis::DataType::Gene);
	}
	else if (field_type == DataField::DataType::Atac) {

		this->sparse(DifferentialAnalysis::DataType::Peak);
	}
	else if (field_type == DataField::DataType::Trans) {

		this->sparse(DifferentialAnalysis::DataType::GeneActivity);
	}
};

bool DifferentialAnalysisWorker::work() {

	if (this->mode_ == WorkMode::SingleCellRna) {
		scrna_mode();
	}
	else if (this->mode_ == WorkMode::BulkRna) {
		bulkrna_mode();
	}
	else if (this->mode_ == WorkMode::SingleCellMultiome) {
		scmultiome_mode();
	}
	else if (this->mode_ == WorkMode::ChromVAR) {
		chromvar_mode();
	}
	else if (this->mode_ == WorkMode::SingleCellAtac) {
		scatac_mode();
	}
	else {
		cicero_mode();
	}

	return true;
};

void DifferentialAnalysisWorker::run() {

	G_TASK_LOG("Start differential analysis...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_differential_analysis_ready(this->res_, this->comparison_[1] + " vs. " + this->comparison_[2]);

	G_TASK_LOG("Differential analysis finished.");

	G_TASK_END;
}
