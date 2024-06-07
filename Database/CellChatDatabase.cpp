#include "CellChatDatabase.h"
#include "FileIO.h"
#include "WilcoxTest.h"
#include "Custom.h"

void CellChatDatabase::load_human_database() {

	if (this->interaction_database_.contains("HUMAN"))return;

	auto ptr = read_sv(FILE_CELLCHAT_HUMAN_INTERACTION, ',');
	this->interaction_database_["HUMAN"] = *ptr;
	delete ptr;

	ptr = read_sv(FILE_CELLCHAT_HUMAN_COMPLEX, ',');
	this->complex_database_["HUMAN"] = *ptr;
	delete ptr;

	ptr = read_sv(FILE_CELLCHAT_HUMAN_GENEINFO, ',');
	this->gene_information_database_["HUMAN"] = *ptr;
	delete ptr;

	ptr = read_sv(FILE_CELLCHAT_HUMAN_COFACTOR, ',');
	this->cofactor_database_["HUMAN"] = *ptr;
	delete ptr;

	//    auto p = read_sparse_int(":/soap/Database/CellChatDatabase/CELLCHAT_PPI_HUMAN");
	//    this->protein_protein_interaction_database_["HUMAN"] = *p;
	//    delete p;
};

void CellChatDatabase::load_mouse_database() {

	if (this->interaction_database_.contains("MOUSE"))return;

	auto ptr = read_sv(FILE_CELLCHAT_MOUSE_INTERACTION, ',');
	this->interaction_database_["MOUSE"] = *ptr;
	delete ptr;

	ptr = read_sv(FILE_CELLCHAT_MOUSE_COMPLEX, ',');
	this->complex_database_["MOUSE"] = *ptr;
	delete ptr;

	ptr = read_sv(FILE_CELLCHAT_MOUSE_GENEINFO, ',');
	this->gene_information_database_["MOUSE"] = *ptr;
	delete ptr;

	ptr = read_sv(FILE_CELLCHAT_MOUSE_COFACTOR, ',');
	this->cofactor_database_["MOUSE"] = *ptr;
	delete ptr;
	//    this->protein_protein_interaction_database_["MOUSE"] = read_sparse_int(":/soap/Database/CellChatDatabase/CELLCHAT_this->protein_protein_interaction_database_MOUSE");
};

void CellChatDatabase::subset_database(const QString& type) {

	this->interactions_->row_slice(_Cs equal(this->interactions_->get_qstring("annotation"), type));

};

// TO DO : reorganize code structure
void CellChatDatabase::subset_data() {

	QStringList gene_ligand = this->interactions_->get_qstring("ligand");
	QStringList complex = _Cs set_difference(gene_ligand, this->gene_information_->get_qstring("Symbol"));

	gene_ligand = _Cs intersect(gene_ligand, this->gene_information_->get_qstring("Symbol"));

	Eigen::ArrayX<bool> filter = _Cs in(this->complex_->rownames_, complex);

	QStringList subunits;
	for (const auto& colname : this->complex_->colnames_) {
		subunits << _Cs sliced(this->complex_->get_const_qstring_reference(colname), filter);
	}
	subunits.removeAll("");
	gene_ligand = _Cs unique(gene_ligand << subunits);

	QStringList gene_receptor = this->interactions_->get_const_qstring_reference("receptor");

	complex = _Cs set_difference(gene_receptor, this->gene_information_->get_qstring("Symbol"));
	gene_receptor = _Cs intersect(gene_receptor, this->gene_information_->get_qstring("Symbol"));
	filter = _Cs in(this->complex_->rownames_, complex);
	subunits.clear();
	for (const auto& colname : this->complex_->colnames_) {
		subunits << _Cs sliced(this->complex_->get_const_qstring_reference(colname), filter);
	}
	subunits.removeAll("");
	gene_receptor = _Cs unique(gene_receptor << subunits);

	QStringList cofactor;
	cofactor << this->interactions_->get_qstring("agonist") << this->interactions_->get_qstring("antagonist") << this->interactions_->get_qstring("co_A_receptor") << this->interactions_->get_qstring("co_I_receptor");
	cofactor.removeAll("");
	cofactor = _Cs unique(cofactor);
	filter = _Cs in(this->cofactor_->rownames_, cofactor);
	subunits.clear();
	for (const auto& colname : this->cofactor_->colnames_) {
		subunits << _Cs sliced(this->cofactor_->get_qstring(colname), filter);
	}
	subunits.removeAll("");
	QStringList gene_use = _Cs unique(gene_ligand << gene_receptor << subunits);

	filter = _Cs in(this->normalized_->rownames_, gene_use);
	this->normalized_->row_slice(filter);
};

void CellChatDatabase::identify_overexpressed_interactions() {
	QStringList gene_use = this->normalized_->rownames_;
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
		if (_Cs intersect_length(subunits, feature_significant) > 0 && _Cs set_difference_length(subunits, gene_use) == 0) {
#pragma omp critical
			index_significant << i;
		}
	}

	QStringList complex = _Cs reordered(this->complex_->rownames_, index_significant) << feature_significant;
	index_significant.clear();

	int n_interaction = this->interactions_->rows();

#pragma omp parallel for
	for (int i = 0; i < n_interaction; ++i) {
		if (complex.contains(this->interactions_->get_qstring("ligand", i)) && complex.contains(this->interactions_->get_qstring("receptor", i))) {
#pragma omp critical
			index_significant << i;
		}
	}
	this->significant_ligand_receptor_ = this->interactions_->row_reordered(index_significant);
};

void CellChatDatabase::identify_overexpressed_gene() {
	QStringList factors = _Cs unique(this->metadata_);
	QStringList all_gene_names;
	QVector<double> all_log2_fc, all_p_adjusted, allp1, allp2;
	QStringList cluster_name;

	for (auto factor : factors) {
		QStringList gene_names = this->normalized_->rownames_;
		Eigen::ArrayX<bool> group1 = _Cs equal(this->metadata_, factor), group2;
		group2 = !group1;
		int gene_number = gene_names.size();

		Eigen::ArrayX<bool> computed = Eigen::ArrayX<bool>::Constant(gene_number, false);
		Eigen::ArrayXd p_adjusted = Eigen::ArrayXd::Zero(gene_number);
		Eigen::ArrayXd log2_fold_change = Eigen::ArrayXd::Zero(gene_number);
		Eigen::MatrixXd group_1_matrix = _Cs col_sliced(this->normalized_->mat_, group1);
		Eigen::MatrixXd group_2_matrix = _Cs col_sliced(this->normalized_->mat_, group2);
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

		gene_names = _Cs sliced(gene_names, computed);
		log2_fold_change = _Cs sliced(log2_fold_change, computed);
		p_adjusted = _Cs sliced(p_adjusted, computed);
		gene_number = gene_names.size();
		percentage_group1 = _Cs sliced(percentage_group1, computed);
		percentage_group2 = _Cs sliced(percentage_group2, computed);
		p_adjusted = _Cs adjust_p_value(p_adjusted, this->p_adjust_method_);
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
};

// TO DO:
void CellChatDatabase::project_data() {
	int nrow = this->protein_protein_interaction_->rownames_.size(), ncol = this->normalized_->mat_.cols();
	Eigen::MatrixXd data_in_new_space = Eigen::MatrixXd::Zero(nrow, ncol);

	for (int i = 0; i < nrow; ++i) {
		QString feature = this->protein_protein_interaction_->rownames_[i];
		int index = this->normalized_->rownames_.indexOf(feature);
		if (index != -1) {
			data_in_new_space.row(i) = this->normalized_->mat_.row(index);
		}
	}
	Eigen::SparseMatrix<double> mat = _Cs row_normalize(this->protein_protein_interaction_->mat_.cast<double>(), 1);
	Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nrow, nrow);
	double alpha = 0.5;
	eye = eye - alpha * mat;
	mat.resize(0, 0);
	//Eigen::MatrixXd this->projected_ = eye.colPivHouseholderQr().solve((1 - alpha) * data_in_new_space);
	Eigen::MatrixXd projected = eye.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
		.solve((1 - alpha) * data_in_new_space).transpose();
	data_in_new_space.resize(0, 0);
	QStringList genes_in_both = _Cs intersect(this->normalized_->rownames_, this->protein_protein_interaction_->rownames_);
	this->projected_ = this->normalized_->to_dense();
	data_in_new_space = this->normalized_->mat_.toDense();
	for (int i = 0; i < nrow; ++i) {
		QString feature = this->protein_protein_interaction_->rownames_[i];
		int index = this->projected_.rownames_.indexOf(feature);
		if (index != -1) {
			this->projected_.mat_.row(index) = projected.row(i);
		}
	}
};

void CellChatDatabase::compute_pathway_communication_probability() {
	QStringList pathway_list = this->significant_ligand_receptor_.get_qstring("pathway_name");
	QStringList pathways = _Cs unique(pathway_list);

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

	auto order = _Cs order(probs, true);
	this->significant_pathway_ = _Cs reordered(this->significant_pathway_, order);
};

void CellChatDatabase::compute_communication_probability() {
	QStringList gene_ligand = this->significant_ligand_receptor_.get_qstring("ligand"), gene_receptor = this->significant_ligand_receptor_.get_qstring("receptor");
	DenseDouble data_use = this->projected_;
	data_use.mat_ = data_use.mat_ / data_use.mat_.maxCoeff();

	Eigen::ArrayXXd data_ligand = compute_ligand_receptor_expression(gene_ligand, data_use);
	Eigen::ArrayXXd data_receptor = compute_ligand_receptor_expression(gene_receptor, data_use);
	Eigen::ArrayXXd data_receptor_co_a_receptor = compute_coreceptor_expression(data_use, "A");
	Eigen::ArrayXXd data_receptor_co_i_receptor = compute_coreceptor_expression(data_use, "I");

	data_receptor = data_receptor * data_receptor_co_a_receptor / data_receptor_co_i_receptor;

	int n_ligand = gene_ligand.size(), n_receptor = gene_receptor.size(), n_cell = data_use.mat_.cols(), n_ligand_receptor = this->significant_ligand_receptor_.rows();

	Eigen::ArrayXXd data_ligand_average(n_ligand, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = _Cs equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = _Cs col_sliced(data_ligand, filter);
		Eigen::ArrayXd means(n_ligand);
		for (int j = 0; j < n_ligand; ++j) {
			means[j] = _Cs trimean(exp.row(j));
		}
		data_ligand_average.col(i) = means;
	}

	Eigen::ArrayXXd dataravg(n_receptor, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = _Cs equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = _Cs col_sliced(data_receptor, filter);
		Eigen::ArrayXd means(n_receptor);
		for (int j = 0; j < n_receptor; ++j) {
			means[j] = _Cs trimean(exp.row(j));
		}
		dataravg.col(i) = means;
	}

	Eigen::ArrayXXd data_ligand_binary = (data_ligand > 0).cast<double>();
	Eigen::ArrayXXd data_receptor_binary = (data_receptor > 0).cast<double>();

	Eigen::ArrayXXd data_ligand_average2(n_ligand, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = _Cs equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = _Cs col_sliced(data_ligand_binary, filter);
		Eigen::ArrayXd means = exp.rowwise().sum() / n_cell;
		data_ligand_average2.col(i) = means;
	}

	Eigen::ArrayXXd dataravg2(n_receptor, this->n_levels_);

#pragma omp parallel for
	for (int i = 0; i < this->n_levels_; ++i) {
		Eigen::ArrayX<bool> filter = _Cs equal(this->metadata_, this->levels_[i]);
		Eigen::ArrayXXd exp = _Cs col_sliced(data_receptor_binary, filter);
		Eigen::ArrayXd means = exp.rowwise().sum() / n_cell;
		dataravg2.col(i) = means;
	}

	QStringList agonists = this->significant_ligand_receptor_.get_qstring("agonist");
	agonists.removeAll("");

	QStringList antagonists = this->significant_ligand_receptor_.get_qstring("antagonist");
	antagonists.removeAll("");

	QVector<int> index_agonist = _Cs which(_Cs not_equal(this->significant_ligand_receptor_.get_qstring("agonist"), QString("")));
	QVector<int> index_antagonist = _Cs which(_Cs not_equal(this->significant_ligand_receptor_.get_qstring("antagonist"), QString("")));

	Eigen::ArrayXXi permutation(n_cell, this->n_boot_);
	for (int i = 0; i < this->n_boot_; ++i) {
		QVector<int> tmp = _Cs seq_n(0, n_cell);
		std::shuffle(tmp.begin(), tmp.end(), random_engine_);
		permutation.col(i) = _Cs cast<Eigen::ArrayX>(tmp);
	}
	Eigen::ArrayX<Eigen::ArrayXXd> prob(n_ligand_receptor), pval(n_ligand_receptor);
	Eigen::ArrayXXd pnull;

	for (int i = 0; i < n_ligand_receptor; ++i) {
		Eigen::VectorXd l = data_ligand_average.row(i), r = dataravg.row(i);
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
			Eigen::ArrayXd datali = data_ligand.row(i), datari = data_receptor.row(i), datal2i = data_ligand_binary.row(i), datar2i = data_receptor_binary.row(i);
			QList<Eigen::ArrayXXd> Pboot(this->n_boot_);

#pragma omp parallel for
			for (int k = 0; k < this->n_boot_; ++k) {
				Eigen::VectorXd data_ligand_average_b(this->n_levels_), data_receptor_average_b(this->n_levels_);
				QStringList new_group = _Cs reordered(this->metadata_, Eigen::ArrayXi(permutation.col(k)));
				for (int j = 0; j < this->n_levels_; ++j) {
					Eigen::ArrayX<bool> filter = _Cs equal(new_group, this->levels_[j]);
					double exp = _Cs trimean(_Cs sliced(datali, filter));
					data_ligand_average_b[j] = exp;
				}
				for (int j = 0; j < this->n_levels_; ++j) {
					Eigen::ArrayX<bool> filter = _Cs equal(new_group, this->levels_[j]);
					double exp = _Cs trimean(_Cs sliced(datari, filter));
					data_receptor_average_b[j] = exp;
				}
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
				P = P + ((Pboot[j] - pnull) >= 0).cast<double>();
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

void CellChatDatabase::subset_communication() {
	double threshold = 0.05;
	QStringList source, target, interaction_name, interaction_name2, pathway_name, ligand, receptor, annotation, evidence;
	int n_ligand_receptor = this->significant_ligand_receptor_.rows();
	QStringList int_name = this->significant_ligand_receptor_.rownames_;
	for (int i = 0; i < n_ligand_receptor; ++i) {
		for (int j = 0; j < this->n_levels_; ++j) {
			for (int k = 0; k < this->n_levels_; ++k) {
				if (this->probability_[i](j, k) > 0 && this->p_value_[i](j, k) < threshold) {
					source << this->levels_[j];
					target << this->levels_[k];
					interaction_name << int_name[i];
				}
			}
		}
	}
	auto order = _Cs index_of(interaction_name, int_name);
	interaction_name2 = _Cs reordered(this->significant_ligand_receptor_.get_qstring("interaction_name2"), order);
	pathway_name = _Cs reordered(this->significant_ligand_receptor_.get_qstring("pathway_name"), order);
	ligand = _Cs reordered(this->significant_ligand_receptor_.get_qstring("ligand"), order);
	receptor = _Cs reordered(this->significant_ligand_receptor_.get_qstring("receptor"), order);
	annotation = _Cs reordered(this->significant_ligand_receptor_.get_qstring("annotation"), order);
	evidence = _Cs reordered(this->significant_ligand_receptor_.get_qstring("evidence"), order);
	// to do
};

Eigen::ArrayXd CellChatDatabase::compute_antagonist_group_expression(const DenseDouble& data_use, int index, const QStringList& group) {
	QString antagonist = this->significant_ligand_receptor_.get_qstring("antagonist", index);
	QStringList antagonist_index_value = this->cofactor_->get_row(antagonist);
	antagonist_index_value.removeAll("");
	antagonist_index_value = _Cs intersect(antagonist_index_value, data_use.rownames_);
	Eigen::ArrayXd data_antagonist(this->n_levels_);
	QVector<int> params = _Cs seq_n(0, this->n_levels_);
	if (antagonist_index_value.size() == 1) {
		Eigen::ArrayXd antago = data_use.get_row(antagonist_index_value[0]);
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = _Cs equal(group, this->levels_[i]);
			double exp = _Cs sliced(antago, filter).mean();
			data_antagonist[i] = exp;
		}
		data_antagonist = 1 + data_antagonist / (0.5 + data_antagonist);
	}
	else if (antagonist_index_value.size() > 1) {
		Eigen::ArrayXXd antago = data_use.get_rows(antagonist_index_value);
		int n_ligand = antago.rows();
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = _Cs equal(group, this->levels_[i]);
			Eigen::ArrayXXd exp = _Cs col_sliced(antago, filter);
			Eigen::ArrayXd means(n_ligand);
			for (int j = 0; j < n_ligand; ++j) {
				means[j] = _Cs trimean(exp.row(j));
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

Eigen::ArrayXd CellChatDatabase::compute_agonist_group_expression(const DenseDouble& data_use, int index, const QStringList& group) {
	QString agonist = this->significant_ligand_receptor_.get_qstring("agonist", index);
	QStringList agonist_index_value = this->cofactor_->get_row(agonist);
	agonist_index_value.removeAll("");
	agonist_index_value = _Cs intersect(agonist_index_value, data_use.rownames_);
	Eigen::ArrayXd data_agonist(this->n_levels_);
	if (agonist_index_value.size() == 1) {
		Eigen::ArrayXd ago = data_use.get_row(agonist_index_value[0]);
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = _Cs equal(group, this->levels_[i]);
			double exp = _Cs sliced(ago, filter).mean();
			data_agonist[i] = exp;
		}
		data_agonist = 1 + data_agonist / (0.5 + data_agonist);
	}
	else if (agonist_index_value.size() > 1) {
		Eigen::ArrayXXd ago = data_use.get_rows_array(agonist_index_value);
		int n_ligand = ago.rows();
		for (int i = 0; i < this->n_levels_; ++i) {
			Eigen::ArrayX<bool> filter = _Cs equal(group, this->levels_[i]);
			Eigen::ArrayXXd exp = _Cs col_sliced(ago, filter);
			Eigen::ArrayXd means(n_ligand);
			for (int j = 0; j < n_ligand; ++j) {
				means[j] = _Cs trimean(exp.row(j));
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

Eigen::ArrayXXd CellChatDatabase::compute_coreceptor_expression(const DenseDouble& data_use, const QString& type) {
	QStringList coreceptor;
	if (type == "A") {
		coreceptor = this->significant_ligand_receptor_.get_qstring("co_A_receptor");
	}
	else if (type == "I") {
		coreceptor = this->significant_ligand_receptor_.get_qstring("co_I_receptor");
	}
	QVector<int> index_coreceptor = _Cs which(_Cs not_equal(coreceptor, QString("")));

	int ncol = data_use.mat_.cols(), nrow = coreceptor.size();
	Eigen::ArrayXXd data_coreceptor = Eigen::ArrayXXd::Ones(nrow, ncol);
	coreceptor.removeAll("");
	if (coreceptor.size() > 0) {
		nrow = coreceptor.size();
		CustomMatrix coreceptor_index = this->cofactor_->row_reordered(coreceptor);
		Eigen::ArrayXXd data_coreceptor_index = Eigen::ArrayXXd::Ones(nrow, ncol);
#pragma omp parallel for
		for (int i = 0; i < nrow; ++i) {
			QStringList sub = coreceptor_index.get_row(i);
			sub.removeAll("");
			sub.removeAll("");
			sub = _Cs intersect(sub, data_use.rownames_);
			Eigen::ArrayXd tmp;
			if (sub.size() == 1) {
				tmp = data_use.get_row(sub[0]) + 1;
			}
			else if (sub.size() > 1) {
				tmp = (data_use.get_rows_array(sub) + 1).colwise().prod();
			}
			data_coreceptor_index.row(i) = tmp;
		}

		data_coreceptor(index_coreceptor, Eigen::all) = data_coreceptor_index;
	}
	return data_coreceptor;
};

Eigen::ArrayXXd CellChatDatabase::compute_complex_expression(const DenseDouble& data_use, const QStringList& complex) {
	CustomMatrix subunits = this->complex_->row_reordered(complex);
	int ncol = data_use.mat_.cols(), nrow = complex.size();
	Eigen::ArrayXXd data_complex(nrow, ncol);
#pragma omp parallel for
	for (int i = 0; i < nrow; ++i) {
		QStringList sub = subunits.get_row(i);
		sub.removeAll("");
		sub = _Cs intersect(sub, data_use.rownames_);
		Eigen::ArrayXd tmp = log(data_use.get_rows_array(sub)).colwise().mean().array().exp();
		data_complex.row(i) = tmp;
	}
	return data_complex;
};

Eigen::ArrayXXd CellChatDatabase::compute_ligand_receptor_expression(const QStringList& gene_ligand_receptor, const DenseDouble& data_use) {
	int n_ligand_receptor = gene_ligand_receptor.size(), n_cell = data_use.mat_.cols();
	auto index_single_ligand = _Cs which(_Cs in(gene_ligand_receptor, data_use.rownames_));
	auto data_average = data_use.row_reordered(_Cs reordered(gene_ligand_receptor, index_single_ligand));
	Eigen::ArrayXXd data_ligand_average = Eigen::ArrayXXd::Zero(n_ligand_receptor, n_cell);
	data_ligand_average(index_single_ligand, Eigen::all) = data_average.mat_;
	auto index_complex_ligand = _Cs set_difference(_Cs seq_n( 0, n_ligand_receptor), index_single_ligand);
	if (index_complex_ligand.size() > 0) {
		QStringList complex = _Cs reordered(gene_ligand_receptor, index_complex_ligand);
		Eigen::ArrayXXd data_complex = this->compute_complex_expression(data_use, complex);
		data_ligand_average(index_complex_ligand, Eigen::all) = data_complex;
	}
	return data_ligand_average;
};

void CellChatDatabase::filter_levels(SparseDouble* normalized, const QStringList& metadata, int minimum_cell_number) {
	auto metadata_table = _Cs table(this->metadata_);
	QStringList filtered_levels;
	for (const auto& level : metadata_table.keys()) {
		if (metadata_table[level] < minimum_cell_number) {
			filtered_levels << level;
		}
	}
	if (filtered_levels.size() > 0) {
		auto filter = !_Cs in(this->metadata_, filtered_levels);
		this->normalized_->col_slice(filter);
		this->metadata_ = _Cs sliced(this->metadata_, filter);
		this->levels_ = _Cs unique(this->metadata_);
		this->n_levels_ = this->levels_.size();
	}
	else {
		this->metadata_ = this->metadata_;
		this->levels_ = _Cs unique(this->metadata_);
		this->n_levels_ = this->levels_.size();
	}
	this->normalized_ = this->normalized_;
};

void CellChatDatabase::clear() {
	this->normalized_ = nullptr;
	this->projected_.clear();
	this->metadata_.clear();
	this->levels_.clear();
	this->markers_.clear();
	this->significant_ligand_receptor_.clear();
	this->probability_.resize(0);
	this->p_value_.resize(0);
	this->pathway_probability_.clear();
	this->counts_.resize(0, 0);
	this->weights_.resize(0, 0);
	this->significant_pathway_.clear();
};

CellChat CellChatDatabase::cellchat(
	const QString& identity, 
	SparseDouble* normalized, 
	soap::Species species, 
	const QStringList& metadata, 
	double minimum_percentage, 
	const QString& p_adjust_method, 
	int random_state, 
	int n_boot, 
	int minimum_cell_number, 
	const QString& annotation_type
) {
	this->filter_levels(normalized, metadata, minimum_cell_number);
	this->minimum_percentage_ = minimum_percentage;
	this->p_adjust_method_ = p_adjust_method;
	this->n_boot_ = n_boot;
	this->random_state_ = random_state;
	random_engine_.seed(random_state);
	if (species == soap::Species::Human) {
		load_human_database();
		this->interactions_ = &this->interaction_database_["HUMAN"];
		this->complex_ = &this->complex_database_["HUMAN"];
		this->gene_information_ = &this->gene_information_database_["HUMAN"];
		this->cofactor_ = &this->cofactor_database_["HUMAN"];
		this->protein_protein_interaction_ = &this->protein_protein_interaction_database_["HUMAN"];
	}
	else {
		load_mouse_database();
		this->interactions_ = &this->interaction_database_["MOUSE"];
		this->complex_ = &this->complex_database_["MOUSE"];
		this->gene_information_ = &this->gene_information_database_["MOUSE"];
		this->cofactor_ = &this->cofactor_database_["MOUSE"];
		this->protein_protein_interaction_ = &this->protein_protein_interaction_database_["MOUSE"];
	}
	if (annotation_type != "ALL") {
		this->subset_database(annotation_type);
	}
	this->subset_data();
	this->identify_overexpressed_gene();

	this->identify_overexpressed_interactions();
	this->project_data();
	this->compute_communication_probability();
	this->compute_pathway_communication_probability();
	CellChat ret(identity, this->significant_ligand_receptor_.get_qstring("interaction_name"), this->significant_pathway_, this->levels_, this->counts_, this->weights_, this->probability_, this->p_value_, this->pathway_probability_);
	this->clear();
	return ret;
}
