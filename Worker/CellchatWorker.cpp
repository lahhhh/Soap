#include "CellchatWorker.h"
#include "CellchatWorker.h"
#include "FileIO.h"
#include "WilcoxTest.h"
#include "Custom.h"

bool CellchatWorker::load_database() {

	if (this->species_ == soap::Species::Human) {

		this->interactions_.reset(read_sv(FILE_CELLCHAT_HUMAN_INTERACTION));

		this->complex_.reset(read_sv(FILE_CELLCHAT_HUMAN_COMPLEX));

		this->gene_information_.reset(read_sv(FILE_CELLCHAT_HUMAN_GENEINFO));

		this->cofactor_.reset(read_sv(FILE_CELLCHAT_HUMAN_COFACTOR));
	}
	else if (this->species_ == soap::Species::Mouse) {

		this->interactions_.reset(read_sv(FILE_CELLCHAT_MOUSE_INTERACTION));

		this->complex_.reset(read_sv(FILE_CELLCHAT_MOUSE_COMPLEX));

		this->gene_information_.reset(read_sv(FILE_CELLCHAT_MOUSE_GENEINFO));

		this->cofactor_.reset(read_sv(FILE_CELLCHAT_MOUSE_COFACTOR));
	}

	if (this->interactions_ == nullptr ||
		this->complex_ == nullptr ||
		this->gene_information_ == nullptr ||
		this->cofactor_ == nullptr) {
		return false;
	}

	return true;
}

void CellchatWorker::subset_database(const QString& type) {

	this->interactions_->row_slice(custom::equal(this->interactions_->get_qstring("annotation"), type));

};

void CellchatWorker::subset_data() {

	QStringList gene_ligand = this->interactions_->get_qstring("ligand");
	QStringList complex = custom::set_difference(gene_ligand, this->gene_information_->get_qstring("Symbol"));
	gene_ligand = custom::intersect(gene_ligand, this->gene_information_->get_qstring("Symbol"));
	Eigen::ArrayX<bool> filter = custom::in(this->complex_->rownames_, complex);
	QStringList subunits;
	for (const auto& colname : this->complex_->colnames_) {
		subunits << custom::sliced(this->complex_->get_const_qstring_reference(colname), filter);
	}
	subunits.removeAll("");
	gene_ligand = custom::unique(gene_ligand << subunits);

	QStringList gene_receptor = this->interactions_->get_const_qstring_reference("receptor");
	complex = custom::set_difference(gene_receptor, this->gene_information_->get_qstring("Symbol"));
	gene_receptor = custom::intersect(gene_receptor, this->gene_information_->get_qstring("Symbol"));
	filter = custom::in(this->complex_->rownames_, complex);
	subunits.clear();
	for (const auto& colname : this->complex_->colnames_) {
		subunits << custom::sliced(this->complex_->get_const_qstring_reference(colname), filter);
	}
	subunits.removeAll("");
	gene_receptor = custom::unique(gene_receptor << subunits);

	QStringList cofactor;
	cofactor << this->interactions_->get_qstring("agonist") << this->interactions_->get_qstring("antagonist") << this->interactions_->get_qstring("co_A_receptor") << this->interactions_->get_qstring("co_I_receptor");
	cofactor.removeAll("");
	cofactor = custom::unique(cofactor);
	filter = custom::in(this->cofactor_->rownames_, cofactor);
	subunits.clear();
	for (const auto& colname : this->cofactor_->colnames_) {
		subunits << custom::sliced(this->cofactor_->get_qstring(colname), filter);
	}
	subunits.removeAll("");
	QStringList gene_use = custom::unique(gene_ligand << gene_receptor << subunits);

	filter = custom::in(this->normalized_.rownames_, gene_use);
	this->normalized_.row_slice(filter);
};

void CellchatWorker::identify_overexpressed_interactions() {
	QStringList gene_use = this->normalized_.rownames_;
	QStringList feature_significant = this->markers_.get_qstring(METADATA_DE_FEATURE_NAME);

	int n_complex = this->complex_->rows();
	QVector<int> index_significant;

#pragma omp parallel for
	for (int i = 0; i < n_complex; ++i) {
		QStringList subunits;
		for (const auto& colname : this->complex_->colnames_) {
			subunits << this->complex_->get_qstring(colname, i);
		}
		subunits.removeAll("");
		if (custom::intersect_length(subunits, feature_significant) > 0 && custom::set_difference_length(subunits, gene_use) == 0) {
		#pragma omp critical
			{
				index_significant << i;
			}
		}
	}

	QStringList complex = custom::reordered(this->complex_->rownames_, index_significant) << feature_significant;
	index_significant.clear();

	int n_interaction = this->interactions_->rows();

#pragma omp parallel for
	for (int i = 0; i < n_interaction; ++i) {
		if (complex.contains(this->interactions_->get_qstring("ligand", i)) && complex.contains(this->interactions_->get_qstring("receptor", i))) {
		#pragma omp critical
			{
				index_significant << i;
			}
		}
	}
	this->significant_ligand_receptor_ = this->interactions_->row_reordered(index_significant);
};

void CellchatWorker::identify_overexpressed_gene() {
	QStringList factors = custom::unique(this->metadata_);
	QStringList all_gene_names;
	QVector<double> all_log2_fc, all_p_adjusted, allp1, allp2;
	QStringList cluster_name;

	for (auto factor : factors) {
		QStringList gene_names = this->normalized_.rownames_;
		Eigen::ArrayX<bool> group1 = custom::equal(this->metadata_, factor), group2;
		group2 = !group1;
		int gene_number = gene_names.size();

		Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(gene_number, false);
		Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(gene_number);
		Eigen::ArrayXd log2_fold_change = Eigen::ArrayXd::Zero(gene_number);
		Eigen::MatrixXd group_1_matrix = custom::col_sliced(this->normalized_.mat_, group1);
		Eigen::MatrixXd group_2_matrix = custom::col_sliced(this->normalized_.mat_, group2);
		Eigen::ArrayXd percentage_group1 = group_1_matrix.rowwise().count().cast<double>() / group_1_matrix.cols();
		Eigen::ArrayXd percentage_group2 = group_2_matrix.rowwise().count().cast<double>() / group_2_matrix.cols();

	#pragma omp parallel for
		for (int i = 0; i < gene_number; ++i) {
			if (percentage_group1[i] < this->minimum_percentage_ && percentage_group2[i] < this->minimum_percentage_) {
				continue;
			}
			Eigen::ArrayXd g1 = group_1_matrix.row(i);
			Eigen::ArrayXd g2 = group_2_matrix.row(i);
			double denominator = g2.mean(), numerator = g1.mean();
			double log2_fold_change_value;
			log2_fold_change_value = log((numerator + 1) / (denominator + 1));

			if (log2_fold_change_value <= 0)
				continue;

			log2_fold_change[i] = log2_fold_change_value;
			computed[i] = true;
			double p_value_adjusted = wilcox_test(g1, g2);
			p_adjusted[i] = p_value_adjusted;
		}

		gene_names = custom::sliced(gene_names, computed);
		log2_fold_change = custom::sliced(log2_fold_change, computed);
		p_adjusted = custom::sliced(p_adjusted, computed);
		gene_number = gene_names.size();
		percentage_group1 = custom::sliced(percentage_group1, computed);
		percentage_group2 = custom::sliced(percentage_group2, computed);
		p_adjusted = custom::adjust_p_value(p_adjusted, this->p_adjust_method_);
		all_gene_names.append(gene_names);
		all_log2_fc.append(QVector<double>(log2_fold_change.begin(), log2_fold_change.end()));
		all_p_adjusted.append(QVector<double>(p_adjusted.begin(), p_adjusted.end()));
		allp1.append(QVector<double>(percentage_group1.begin(), percentage_group1.end()));
		allp2.append(QVector<double>(percentage_group2.begin(), percentage_group2.end()));
		cluster_name.append(QStringList(gene_names.size(), factor));
	}
	this->markers_.set_rownames(all_gene_names);
	this->markers_.update(METADATA_DE_FEATURE_NAME, all_gene_names);
	this->markers_.update(METADATA_DE_LOG2_FOLD_CHANGE, all_log2_fc);
	this->markers_.update(METADATA_DE_ADJUSTED_P_VALUE, all_p_adjusted);
	this->markers_.update(METADATA_DE_PERCENTAGE_1, allp1);
	this->markers_.update(METADATA_DE_PERCENTAGE_2, allp2);
	this->markers_.update("Identity", cluster_name, CustomMatrix::DataType::QStringFactor);

	// only filter p value here
	this->markers_.row_slice(custom::less_than(this->markers_.get_const_double_reference(METADATA_DE_ADJUSTED_P_VALUE), 0.05));
};

// deprecated
void CellchatWorker::project_data() {
	//int nrow = this->protein_protein_interaction_->rownames_.size(), ncol = this->normalized_.mat_.cols();
	//Eigen::MatrixXd data_in_new_space = Eigen::MatrixXd::Zero(nrow, ncol);

	//for (int i = 0; i < nrow; ++i) {
	//	QString feature = this->protein_protein_interaction_->rownames_[i];
	//	int index = this->normalized_.rownames_.indexOf(feature);
	//	if (index != -1) {
	//		data_in_new_space.row(i) = this->normalized_.mat_.row(index);
	//	}
	//}
	//Eigen::SparseMatrix<double> mat = custom::row_normalize(this->protein_protein_interaction_->mat_.cast<double>(), 1);
	//Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nrow, nrow);
	//double alpha = 0.5;
	//eye = eye - alpha * mat;
	//mat.resize(0, 0);
	////Eigen::MatrixXd this->projected_ = eye.colPivHouseholderQr().solve((1 - alpha) * data_in_new_space);
	//Eigen::MatrixXd projected = eye.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
	//	.solve((1 - alpha) * data_in_new_space).transpose();
	//data_in_new_space.resize(0, 0);
	//QStringList genes_in_both = custom::intersect(this->normalized_.rownames_, this->protein_protein_interaction_->rownames_);
	//this->projected_ = this->normalized_.to_dense();
	//data_in_new_space = this->normalized_.mat_.toDense();
	//for (int i = 0; i < nrow; ++i) {
	//	QString feature = this->protein_protein_interaction_->rownames_[i];
	//	int index = this->projected_.rownames_.indexOf(feature);
	//	if (index != -1) {
	//		this->projected_.mat_.row(index) = projected.row(i);
	//	}
	//}
};

void CellchatWorker::compute_pathway_communication_probability() {
	QStringList pathway_list = this->significant_ligand_receptor_.get_qstring("pathway_name");
	QStringList pathways = custom::unique(pathway_list);

	int size = pathway_list.size(), n_ligand_receptor = this->significant_ligand_receptor_.rows();

	QVector<double> probs;
	for (const auto& path : pathways) {
		Eigen::ArrayXXd pathway_probability = Eigen::ArrayXXd::Zero(this->n_levels_, this->n_levels_);
		for (int i = 0; i < size; ++i) {
			if (pathway_list[i] == path) {
				pathway_probability += this->probability_[i];
			}
		}
		double prob = pathway_probability.sum();
		if (prob > 0) {
			this->significant_pathway_ << path;
			probs << prob;
			this->pathway_probability_[path] = pathway_probability;
		}
	}

	auto order = custom::order(probs, true);
	this->significant_pathway_ = custom::reordered(this->significant_pathway_, order);
};

void CellchatWorker::compute_communication_probability() {
	QStringList gene_ligand = this->significant_ligand_receptor_.get_qstring("ligand");
	QStringList	gene_receptor = this->significant_ligand_receptor_.get_qstring("receptor");
	DenseDouble data_use = this->normalized_.to_dense();
	data_use.mat_ /= data_use.mat_.maxCoeff();

	Eigen::ArrayXXd data_ligand = compute_ligand_receptor_expression(gene_ligand, data_use);
	Eigen::ArrayXXd data_receptor = compute_ligand_receptor_expression(gene_receptor, data_use);
	Eigen::ArrayXXd data_receptor_co_a_receptor = compute_coreceptor_expression(data_use, "A");
	Eigen::ArrayXXd data_receptor_co_i_receptor = compute_coreceptor_expression(data_use, "I");

	data_receptor = data_receptor * data_receptor_co_a_receptor / data_receptor_co_i_receptor;

	int n_ligand = gene_ligand.size(),
		n_receptor = gene_receptor.size(),
		n_cell = data_use.mat_.cols(),
		n_ligand_receptor = this->significant_ligand_receptor_.rows();

	Eigen::ArrayXXd data_ligand_average(n_ligand, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = custom::equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = custom::col_sliced(data_ligand, filter);
		Eigen::ArrayXd means(n_ligand);
		for (int j = 0; j < n_ligand; ++j) {
			means[j] = custom::trimean(exp.row(j));
		}
		data_ligand_average.col(i) = means;
	}

	Eigen::ArrayXXd data_receptor_average(n_receptor, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = custom::equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = custom::col_sliced(data_receptor, filter);
		Eigen::ArrayXd means(n_receptor);
		for (int j = 0; j < n_receptor; ++j) {
			means[j] = custom::trimean(exp.row(j));
		}
		data_receptor_average.col(i) = means;
	}

	Eigen::ArrayXXd data_ligand_binary = (data_ligand > 0).cast<double>();
	Eigen::ArrayXXd data_receptor_binary = (data_receptor > 0).cast<double>();

	Eigen::ArrayXXd data_ligand_average2(n_ligand, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = custom::equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = custom::col_sliced(data_ligand_binary, filter);
		Eigen::ArrayXd means = exp.rowwise().sum() / n_cell;
		data_ligand_average2.col(i) = means;
	}

	Eigen::ArrayXXd data_receptor_average2(n_receptor, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = custom::equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = custom::col_sliced(data_receptor_binary, filter);
		Eigen::ArrayXd means = exp.rowwise().sum() / n_cell;
		data_receptor_average2.col(i) = means;
	}

	QStringList agonists = this->significant_ligand_receptor_.get_qstring("agonist");
	agonists.removeAll("");

	QStringList antagonists = this->significant_ligand_receptor_.get_qstring("antagonist");
	antagonists.removeAll("");

	QVector<int> index_agonist = custom::which(custom::not_equal(this->significant_ligand_receptor_.get_qstring("agonist"), QString("")));
	QVector<int> index_antagonist = custom::which(custom::not_equal(this->significant_ligand_receptor_.get_qstring("antagonist"), QString("")));

	Eigen::ArrayXXi permutation(n_cell, this->n_boot_);
	for (int i = 0; i < this->n_boot_; ++i) {
		QVector<int> tmp = custom::seq_n(0, n_cell);
		std::shuffle(tmp.begin(), tmp.end(), random_engine_);
		permutation.col(i) = custom::cast<Eigen::ArrayX>(tmp);
	}
	Eigen::ArrayX<Eigen::ArrayXXd> prob(n_ligand_receptor), pval(n_ligand_receptor);
	Eigen::ArrayXXd pnull;


	for (int i = 0; i < n_ligand_receptor; ++i) {
		Eigen::VectorXd l = data_ligand_average.row(i), r = data_receptor_average.row(i);
		Eigen::ArrayXXd data_ligand_receptor = (l * r.transpose()).array();
		Eigen::ArrayXXd P1 = data_ligand_receptor / (0.5 + data_ligand_receptor), P2, P3, P4;

		if (P1.sum() == 0) {
			pnull = P1;
			prob[i] = pnull;
			pval[i] = Eigen::ArrayXXd::Ones(this->n_levels_, this->n_levels_);
		}
		else {
			if (index_agonist.contains(i)) {
				Eigen::VectorXd agonist = this->compute_agonist_group_expression(data_use, i, this->metadata_);
				P2 = agonist * agonist.transpose();
			}
			else {
				P2 = Eigen::ArrayXXd::Ones(this->n_levels_, this->n_levels_);
			}
			if (index_antagonist.contains(i)) {
				Eigen::VectorXd antagonist = this->compute_antagonist_group_expression(data_use, i, this->metadata_);
				P3 = antagonist * antagonist.transpose();
			}
			else {
				P3 = Eigen::ArrayXXd::Ones(this->n_levels_, this->n_levels_);
			}
			pnull = P1 * P2 * P3;
			prob[i] = pnull;
			Eigen::ArrayXd datali = data_ligand.row(i);
			Eigen::ArrayXd datari = data_receptor.row(i);
			Eigen::ArrayXd data_coa_i = data_receptor_co_a_receptor.row(i);
			Eigen::ArrayXd data_coi_i = data_receptor_co_i_receptor.row(i);
			Eigen::ArrayXd datal2i = data_ligand_binary.row(i);
			Eigen::ArrayXd datar2i = data_receptor_binary.row(i);
			QList<Eigen::ArrayXXd> Pboot(this->n_boot_);

		#pragma omp parallel for
			for (int k = 0; k < this->n_boot_; ++k) {

				// calculate only 1 pair lr
				Eigen::VectorXd data_ligand_average_b(this->n_levels_),
					data_receptor_average_b(this->n_levels_);

				Eigen::VectorXd data_receptor_average_b_co_A_receptor(this->n_levels_),
					data_receptor_average_b_co_I_receptor(this->n_levels_);

				QStringList new_group = custom::reordered(this->metadata_, Eigen::ArrayXi(permutation.col(k)));

				for (int j = 0; j < this->n_levels_; ++j) {

					Eigen::ArrayX<bool> filter = custom::equal(new_group, this->levels_[j]);

					data_ligand_average_b[j] = custom::trimean(custom::sliced(datali, filter));
					data_receptor_average_b[j] = custom::trimean(custom::sliced(datari, filter));
					data_receptor_average_b_co_A_receptor[j] = custom::trimean(custom::sliced(data_coa_i, filter));
					data_receptor_average_b_co_I_receptor[j] = custom::trimean(custom::sliced(data_coi_i, filter));
				}

				data_receptor_average_b = data_receptor_average_b.array() * data_receptor_average_b_co_A_receptor.array() /
					data_receptor_average_b_co_I_receptor.array();

				Eigen::ArrayXXd data_lr_b = (data_ligand_average_b * data_receptor_average_b.transpose()).array(), P1boot, P2boot, P3boot, Pnullboot;
				P1boot = data_lr_b / (0.5 + data_lr_b);

				if (index_agonist.contains(i)) {
					Eigen::VectorXd agonist = this->compute_agonist_group_expression(data_use, i, new_group);
					P2boot = agonist * agonist.transpose();
				}
				else {
					P2boot = Eigen::ArrayXXd::Ones(this->n_levels_, this->n_levels_);
				}

				if (index_antagonist.contains(i)) {
					Eigen::VectorXd antagonist = this->compute_antagonist_group_expression(data_use, i, new_group);
					P3boot = antagonist * antagonist.transpose();
				}
				else {
					P3boot = Eigen::ArrayXXd::Ones(this->n_levels_, this->n_levels_);
				}
				Pnullboot = P1boot * P2boot * P3boot;
				Pboot[k] = Pnullboot;
			}

			Eigen::ArrayXXd P = Eigen::ArrayXXd::Zero(this->n_levels_, this->n_levels_);
			for (int j = 0; j < this->n_boot_; ++j) {
				P += ((Pboot[j] - pnull) > 0).cast<double>();
			}
			pval[i] = P / this->n_boot_;
		}
	}
	this->counts_ = Eigen::ArrayXXi::Zero(this->n_levels_, this->n_levels_);
	this->weights_ = Eigen::ArrayXXd::Zero(this->n_levels_, this->n_levels_);
	for (int i = 0; i < n_ligand_receptor; ++i) {
		for (int j = 0; j < this->n_levels_; ++j) {
			for (int k = 0; k < this->n_levels_; ++k) {
				if (prob[i](j, k) == 0) {
					pval[i](j, k) = 1;
				}
			}
		}
	}
	for (int i = 0; i < n_ligand_receptor; ++i) {
		for (int j = 0; j < this->n_levels_; ++j) {
			for (int k = 0; k < this->n_levels_; ++k) {
				if (pval[i](j, k) > 0.05) {
					prob[i](j, k) = 0;
				}
			}
		}
	}
	for (int i = 0; i < n_ligand_receptor; ++i) {
		for (int j = 0; j < this->n_levels_; ++j) {
			for (int k = 0; k < this->n_levels_; ++k) {
				if (prob[i](j, k) > 0) {
					this->counts_(j, k) += 1;
					this->weights_(j, k) += prob[i](j, k);
				}
			}
		}
	}
	this->probability_ = prob;
	this->p_value_ = pval;
};

void CellchatWorker::subset_communication() {
	double threshold = 0.05;
	QStringList source, target, interaction_name, interaction_name2, pathway_name, ligand, receptor, annotation, evidence;
	int n_ligand_receptor = this->significant_ligand_receptor_.rows();
	QStringList int_name = this->significant_ligand_receptor_.rownames_;
	QVector<double> pval, prob;
	for (int i = 0; i < n_ligand_receptor; ++i) {
		for (int j = 0; j < this->n_levels_; ++j) {
			for (int k = 0; k < this->n_levels_; ++k) {

				double p_val = this->p_value_[i](j, k);
				double probability = this->probability_[i](j, k);

				if (probability > 0 && p_val < threshold) {
					source << this->levels_[j];
					target << this->levels_[k];
					interaction_name << int_name[i];
					pval << p_val;
					prob << probability;
				}
			}
		}
	}
	auto order = custom::index_of(interaction_name, int_name);
	interaction_name2 = custom::reordered(this->significant_ligand_receptor_.get_qstring("interaction_name_2"), order);
	pathway_name = custom::reordered(this->significant_ligand_receptor_.get_qstring("pathway_name"), order);
	ligand = custom::reordered(this->significant_ligand_receptor_.get_qstring("ligand"), order);
	receptor = custom::reordered(this->significant_ligand_receptor_.get_qstring("receptor"), order);
	annotation = custom::reordered(this->significant_ligand_receptor_.get_qstring("annotation"), order);
	evidence = custom::reordered(this->significant_ligand_receptor_.get_qstring("evidence"), order);

	int n_lr = order.size();

	this->interaction_summary_.set_nrow(n_lr);
	this->interaction_summary_.update("source", source, CustomMatrix::DataType::QStringFactor);
	this->interaction_summary_.update("target", target, CustomMatrix::DataType::QStringFactor);
	this->interaction_summary_.update("p value", pval);
	this->interaction_summary_.update("probability", prob);
	this->interaction_summary_.update("ligand", ligand);
	this->interaction_summary_.update("receptor", receptor);
	this->interaction_summary_.update("interaction name", interaction_name);
	this->interaction_summary_.update("interaction name 2", interaction_name2);
	this->interaction_summary_.update("pathway name", pathway_name);
	this->interaction_summary_.update("annotation", annotation);
	this->interaction_summary_.update("evidence", evidence);

	QList<std::tuple<QString, QString, QString>> levels;

	for (int i = 0; i < n_lr; ++i) {
		levels << std::make_tuple(source[i], target[i], pathway_name[i]);
	}

	levels = custom::unique(levels);

	QStringList c = custom::paste("-", source, target, pathway_name);

	int n_level = levels.size();

	pathway_name.clear();
	QStringList pathway_source, pathway_target;
	QVector<double> pathway_pval, pathway_prob;

	for (int i = 0; i < n_level; ++i) {

		QString l = std::get<0>(levels[i]) + "-" + std::get<1>(levels[i]) + "-" + std::get<2>(levels[i]);

		auto filter = custom::equal(c, l);

		double p = custom::mean(custom::sliced(pval, filter));
		double probability = custom::sum(custom::sliced(prob, filter));

		pathway_source << std::get<0>(levels[i]);
		pathway_target << std::get<1>(levels[i]);
		pathway_name << std::get<2>(levels[i]);
		pathway_pval << p;
		pathway_prob << probability;
	}

	this->pathway_summary_.set_nrow(n_level);
	this->pathway_summary_.update("source", pathway_source, CustomMatrix::DataType::QStringFactor);
	this->pathway_summary_.update("target", pathway_target, CustomMatrix::DataType::QStringFactor);
	this->pathway_summary_.update("pathway name", pathway_name);
	this->pathway_summary_.update("p value", pathway_pval);
	this->pathway_summary_.update("probability", pathway_prob);

};

Eigen::ArrayXd CellchatWorker::compute_antagonist_group_expression(
	const DenseDouble& data_use,
	int index,
	const QStringList& group) {

	QString antagonist = this->significant_ligand_receptor_.get_qstring("antagonist", index);
	QStringList antagonist_index_value = this->cofactor_->get_row(antagonist);
	antagonist_index_value.removeAll("");
	antagonist_index_value = custom::intersect(antagonist_index_value, data_use.rownames_);
	Eigen::ArrayXd data_antagonist(this->n_levels_);
	QVector<int> params = custom::seq_n(0, this->n_levels_);

	if (antagonist_index_value.size() == 1) {
		Eigen::ArrayXd antago = data_use.get_row(antagonist_index_value[0]);
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = custom::equal(group, this->levels_[i]);
			double exp = custom::sliced(antago, filter).mean();
			data_antagonist[i] = exp;
		}
		data_antagonist = 1 + data_antagonist / (0.5 + data_antagonist);
	}
	else if (antagonist_index_value.size() > 1) {
		Eigen::ArrayXXd antago = data_use.get_rows(antagonist_index_value);
		int n_ligand = antago.rows();
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = custom::equal(group, this->levels_[i]);
			Eigen::ArrayXXd exp = custom::col_sliced(antago, filter);
			Eigen::ArrayXd means(n_ligand);
			for (int j = 0; j < n_ligand; ++j) {
				means[j] = custom::trimean(exp.row(j));
			}
			means = 1 + means / (0.5 + means);
			double res = means.prod();
			data_antagonist[i] = res;
		}
	}
	else {
		data_antagonist = Eigen::ArrayXd::Ones(this->n_levels_);
	}
	return data_antagonist;
};

Eigen::ArrayXd CellchatWorker::compute_agonist_group_expression(
	const DenseDouble& data_use,
	int index,
	const QStringList& group) {

	QString agonist = this->significant_ligand_receptor_.get_qstring("agonist", index);
	QStringList agonist_index_value = this->cofactor_->get_row(agonist);
	agonist_index_value.removeAll("");
	agonist_index_value = custom::intersect(agonist_index_value, data_use.rownames_);

	Eigen::ArrayXd data_agonist(this->n_levels_);

	if (agonist_index_value.size() == 1) {
		Eigen::ArrayXd ago = data_use.get_row(agonist_index_value[0]);
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = custom::equal(group, this->levels_[i]);
			double exp = custom::sliced(ago, filter).mean();
			data_agonist[i] = exp;
		}
		data_agonist = 1 + data_agonist / (0.5 + data_agonist);
	}
	else if (agonist_index_value.size() > 1) {
		Eigen::ArrayXXd ago = data_use.get_rows_array(agonist_index_value);
		int n_ligand = ago.rows();
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = custom::equal(group, this->levels_[i]);
			Eigen::ArrayXXd exp = custom::col_sliced(ago, filter);
			Eigen::ArrayXd means(n_ligand);
			for (int j = 0; j < n_ligand; ++j) {
				means[j] = custom::trimean(exp.row(j));
			}
			means = 1 + means / (0.5 + means);
			double res = means.prod();
			data_agonist[i] = res;
		}
	}
	else {
		data_agonist = Eigen::ArrayXd::Ones(this->n_levels_);
	}

	return data_agonist;
};

Eigen::ArrayXXd CellchatWorker::compute_coreceptor_expression(const DenseDouble& data_use, const QString& type) const {
	QStringList coreceptor;
	if (type == "A") {
		coreceptor = this->significant_ligand_receptor_.get_qstring("co_A_receptor");
	}
	else if (type == "I") {
		coreceptor = this->significant_ligand_receptor_.get_qstring("co_I_receptor");
	}
	QVector<int> index_coreceptor = custom::which(custom::not_equal(coreceptor, QString("")));

	int ncol = data_use.mat_.cols(), nrow = coreceptor.size();
	Eigen::ArrayXXd data_coreceptor = Eigen::ArrayXXd::Ones(nrow, ncol);

	for (int i = 0; i < nrow; ++i) {
		if (coreceptor[i] == "") {
			continue;
		}

		QStringList sub = this->cofactor_->get_row(coreceptor[i]);
		sub.removeAll("");
		sub = custom::intersect(sub, data_use.rownames_);
		if (sub.isEmpty()) {
			continue;
		}
		Eigen::ArrayXd tmp;
		if (sub.size() == 1) {
			tmp = data_use.get_row(sub[0]) + 1;
		}
		else {
			tmp = (data_use.get_rows_array(sub) + 1).colwise().prod();
		}
		data_coreceptor.row(i) = tmp;
	}

	return data_coreceptor;
};

Eigen::ArrayXXd CellchatWorker::compute_complex_expression(const DenseDouble& data_use, const QStringList& complex) const {

	int ncol = data_use.mat_.cols(), nrow = complex.size();
	Eigen::ArrayXXd data_complex = Eigen::ArrayXXd::Ones(nrow, ncol);
	for (int i = 0; i < nrow; ++i) {

		QStringList sub = this->complex_->get_row(complex[i]);
		sub.removeAll("");
		sub = custom::intersect(sub, data_use.rownames_);

		if (sub.isEmpty()) {
			continue;
		}

		Eigen::ArrayXd tmp;
		if (sub.size() == 1) {
			tmp = data_use.get_row(sub[0]);
		}
		else {
			tmp = log(data_use.get_rows_array(sub)).colwise().mean().array().exp();
		}

		data_complex.row(i) = tmp;
	}
	return data_complex;
};

Eigen::ArrayXXd CellchatWorker::compute_ligand_receptor_expression(
	const QStringList& gene_ligand_receptor,
	const DenseDouble& data_use) const
{
	int n_ligand_receptor = gene_ligand_receptor.size(), n_cell = data_use.mat_.cols();
	auto index_single_ligand = custom::which(custom::in(gene_ligand_receptor, data_use.rownames_));
	auto data = data_use.row_reordered(custom::reordered(gene_ligand_receptor, index_single_ligand));
	Eigen::ArrayXXd data_lr = Eigen::ArrayXXd::Zero(n_ligand_receptor, n_cell);
	data_lr(index_single_ligand, Eigen::all) = data.mat_;

	auto index_complex_ligand = custom::set_difference(custom::seq_n(0, n_ligand_receptor), index_single_ligand);
	if (index_complex_ligand.size() > 0) {
		QStringList complex = custom::reordered(gene_ligand_receptor, index_complex_ligand);
		Eigen::ArrayXXd data_complex = this->compute_complex_expression(data_use, complex);
		data_lr(index_complex_ligand, Eigen::all) = data_complex;
	}

	return data_lr;
};

void CellchatWorker::filter_levels() {
	auto metadata_table = custom::table(this->metadata_);
	QStringList filtered_levels;
	for (const auto& level : metadata_table.keys()) {
		if (metadata_table[level] < this->minimum_cell_number_) {
			filtered_levels << level;
		}
	}
	if (filtered_levels.size() > 0) {
		auto filter = !custom::in(this->metadata_, filtered_levels);
		this->normalized_ = this->origin_->col_sliced(filter);
		this->metadata_ = custom::sliced(this->metadata_, filter);
	}
	else {
		this->normalized_ = *this->origin_;
	}

	this->levels_ = custom::unique(this->metadata_);
	this->n_levels_ = this->levels_.size();

	this->normalized_ = this->normalized_;
};

bool CellchatWorker::work() {

	this->filter_levels();
	this->random_engine_.seed(this->random_state_);

	if (!this->load_database())
	{
		G_TASK_WARN("Database Loading Failed");
		return false;
	}
	if (this->annotation_type_ != "ALL") {
		this->subset_database(this->annotation_type_);
	}
	this->subset_data();
	this->identify_overexpressed_gene();
	this->identify_overexpressed_interactions();
	//this->project_data();
	this->compute_communication_probability();
	this->compute_pathway_communication_probability();
	this->subset_communication();

	return true;
};

void CellchatWorker::run() {

	G_TASK_LOG("Start cellchat...");

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_cellchat_ready(
		{ this->identity_,
		this->significant_ligand_receptor_.get_qstring("interaction_name"),
		this->significant_pathway_,
		this->levels_,
		this->counts_,
		this->weights_,
		this->probability_,
		this->p_value_,
		this->pathway_probability_,
		this->interaction_summary_,
		this->pathway_summary_
		});

	G_TASK_LOG("Cellchat finished.");

	G_TASK_END;
}
