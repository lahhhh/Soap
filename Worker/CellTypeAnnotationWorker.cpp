#include "CellTypeAnnotationWorker.h"

#include <QFile>

#include "ItemDatabase.h"
#include "DenseDouble.h"

void CellTypeAnnotationWorker::run() {

	this->load_database();

	if (!this->filter_data()) {
		G_TASK_END;
	}

	this->assign_type();

	G_TASK_END;
};

void CellTypeAnnotationWorker::load_database() {

	DenseDouble database;

	ItemDatabase::read_item(FILE_HPCA_CELL_TYPE, database);

	this->gene_name_database_ = database.rownames_;
	
	for (auto&& colname : database.colnames_) {
		QStringList types = colname.split('@');
		this->main_type_database_ << types[0];
		this->sub_type_database_ << types[1];
	}

	this->celltype_expression_ = database.mat_;
};

bool CellTypeAnnotationWorker::filter_data() {

	QStringList gene_names = _Cs intersect(this->counts_->rownames_, this->gene_name_database_);
	if (gene_names.size() < 1000) {
		G_TASK_WARN("No enough gene found in database.");
		return false;
	}

	auto query_index = _Cs index_of(gene_names, this->counts_->rownames_);
	Eigen::SparseMatrix<int> query_counts = _Cs row_reordered(this->counts_->mat_, query_index);

	this->query_expression_ = query_counts.toDense().cast<double>();
	this->query_expression_ = (this->query_expression_.array() + 1.0).log();

	auto database_index = _Cs index_of(gene_names, this->gene_name_database_);

	this->celltype_expression_ = this->celltype_expression_(database_index, Eigen::all).eval();

	double mean_database_exp = this->celltype_expression_.colwise().sum().mean();
	
	int n_cell = this->query_expression_.cols();
	int n_gene = gene_names.size();

	for (int i = 0; i < n_cell; ++i) {

		double diff = this->query_expression_.col(i).sum();

		if (diff == 0.0) {
			continue;
		}

		this->query_expression_.col(i).array() *= (mean_database_exp / diff);
	}

	int n_type = this->celltype_expression_.cols();
	for (int i = 0; i < n_type; ++i) {

		double diff = this->celltype_expression_.col(i).sum();

		if (diff == 0.0) {
			continue;
		}

		this->celltype_expression_.col(i).array() *= (mean_database_exp / diff);
	}

	return true;
};

void CellTypeAnnotationWorker::assign_type() {

	int n_cell = this->query_expression_.cols();
	int n_type = this->celltype_expression_.cols();

	QStringList main_type(n_cell);
	QStringList sub_type(n_cell);

	if (this->annotate_by_cluster_) {
		QStringList levels = _Cs unique(this->cluster_);

		int n_level = levels.size();

		for (int i = 0; i < n_level; ++i) {
			auto cluster_index = _Cs match(this->cluster_, levels[i]);
			Eigen::ArrayXd cluster_exp = this->query_expression_(Eigen::all, cluster_index).rowwise().mean();

			double min_dist = (cluster_exp - this->celltype_expression_.col(0).array()).matrix().squaredNorm();
			int min_index = 0;

			for (int j = 1; j < n_type; ++j) {

				double dist = (cluster_exp - this->celltype_expression_.col(j).array()).matrix().squaredNorm();
				if (dist < min_dist) {
					min_dist = dist;
					min_index = j;
				}
			}

			_Cs assign(main_type, this->main_type_database_[min_index], cluster_index);
			_Cs assign(sub_type, this->sub_type_database_[min_index], cluster_index);
		}
	}
	else {

	#pragma omp parallel for
		for (int i = 0; i < n_cell; ++i) {

			Eigen::ArrayXd cell_exp = this->query_expression_.col(i);
			double min_dist = (cell_exp - this->celltype_expression_.col(0).array()).matrix().squaredNorm();
			int min_index = 0;

			for (int j = 1; j < n_type; ++j) {				

				double dist = (cell_exp - this->celltype_expression_.col(j).array()).matrix().squaredNorm();
				if (dist < min_dist) {
					min_dist = dist;
					min_index = j;
				}
			}		

			main_type[i] = this->main_type_database_[min_index];
			sub_type[i] = this->sub_type_database_[min_index];
		}
	}

	emit x_annotation_ready(
		main_type,
		sub_type,
		this->main_type_name_,
		this->sub_type_name_
	);
};