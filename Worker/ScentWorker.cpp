#include "ScentWorker.h"

#include "custom.h"

#include <QFile>

#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

bool ScentWorker::calculate_max_sr() {

	Eigen::MatrixXd adj = this->ppi_.mat_.cast<double>();

	Spectra::DenseSymMatProd<double> op(adj);
	Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double> > eigs(op, 1, 3);

	eigs.init();

	double max_sr{ 0.0 };

	int nconv = eigs.compute(Spectra::SortRule::LargestAlge);

	if (eigs.info() == Spectra::CompInfo::Successful){
		max_sr = eigs.eigenvalues()(0);
	}
	else {
		G_TASK_WARN("Eigenvalue computation was not successful.");
		return false;
	}

	if (std::isnan(max_sr) || max_sr <= 1.0) {
		G_TASK_WARN("Max SR Computation Failed.");
		return false;
	}

	this->max_sr_ = log(max_sr);
};

bool ScentWorker::do_integ_ppi() {

	QStringList common_gene_names = _Cs intersect(this->counts_->rownames_, this->ppi_.rownames_);

	if (common_gene_names.size() < 5000) {
		G_TASK_WARN("Too few common genes between PPI and expression matrix");
		return false;
	}

	this->expression_.rownames_ = this->counts_->rownames_;
	this->expression_.colnames_ = this->counts_->colnames_;

	this->expression_.mat_ = this->counts_->mat_.toDense().cast<double>();

	_Cs normalize_in_place(this->expression_.mat_);

	this->expression_.mat_ = log2(this->expression_.mat_.array() + 1.1).eval();

	auto map1_idx = _Cs index_of(common_gene_names, this->ppi_.rownames_);
	this->ppi_.reorder(map1_idx, map1_idx);

	auto map2_idx = _Cs index_of(common_gene_names, this->expression_.rownames_);

	this->expression_.row_reorder(map2_idx);

	int n_vertice = this->ppi_.rows();
	igraph_t graph;
	igraph_matrix_t adjacency;
	igraph_integer_t i, j;

	igraph_matrix_init(&adjacency, n_vertice, n_vertice);
	for (j = 0; j < n_vertice; ++j) {
		for (i = 0; i < n_vertice; ++i) {
			MATRIX(adjacency, i, j) = this->ppi_.mat_(i, j);
		}
	}

	igraph_adjacency(&graph, &adjacency, IGRAPH_ADJ_UNDIRECTED, IGRAPH_LOOPS_ONCE);
	igraph_matrix_destroy(&adjacency);

	igraph_vector_t vertice_id;
	igraph_vector_init(&vertice_id, n_vertice);
	for (int j = 0; j < n_vertice; ++j) {
		VECTOR(vertice_id)[j] = j;
	}

	SETVANV(&graph, "vid", &vertice_id);
	igraph_vector_destroy(&vertice_id);

	igraph_graph_list_t components;
	igraph_graph_list_init(&components, 0);
	igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, 0);
	int n_component = igraph_graph_list_size(&components);

	if (n_component > 1) {

		int max_size{ 0 };
		igraph_t* sub;

		for (i = 0; i < n_component; ++i) {
			igraph_t* component = &VECTOR(components)[i];
			igraph_integer_t size = igraph_vcount(component);
			if (size > max_size) {
				max_size = size;
				sub = component;
			}
		}
		igraph_integer_t n_sub_vertice = igraph_vcount(sub);
		igraph_vector_t sub_vertice_id;
		igraph_vector_init(&sub_vertice_id, n_sub_vertice);
		VANV(sub, "vid", &sub_vertice_id);
		QVector<int> sub_vertice(n_sub_vertice);
		for (int i = 0; i < n_sub_vertice; ++i) {
			sub_vertice[i] = VECTOR(sub_vertice_id)[i];
		}
		igraph_vector_destroy(&sub_vertice_id);

		std::ranges::sort(sub_vertice);

		this->ppi_.reorder(sub_vertice, sub_vertice);
		this->expression_.row_reorder(sub_vertice);
	}

	for (int i = 0; i < n_component; ++i) {
		igraph_destroy(&VECTOR(components)[i]);
	}
	igraph_graph_list_destroy(&components);
	igraph_destroy(&graph);

	if (this->ppi_.rows() < 5000) {
		G_TASK_WARN("Too few valid genes");
		return false;
	}

	return true;
};

void ScentWorker::run() {

	igraph_set_attribute_table(&igraph_cattribute_table);

	if (!this->read_ppi_database()) {
		G_TASK_WARN("Database Loading Failed.");
		G_TASK_END;
	}

	if (!this->do_integ_ppi()) {
		G_TASK_END;
	}

	if (!this->calculate_max_sr()) {
		G_TASK_END;
	}

	this->calculate_sr();

	G_TASK_END;
};

bool ScentWorker::read_ppi_database() {

	QFile file(FILE_SCENT_HUMAN_PPI);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		return false;
	}

	QTextStream in(&file);

	QString line = in.readLine();

	QStringList gene_names = _Cs digest(line);
	gene_names.removeAt(0);

	this->ppi_.rownames_ = this->ppi_.colnames_ = gene_names;

	int n_gene = gene_names.size();

	this->ppi_.mat_ = Eigen::MatrixXi::Zero(n_gene, n_gene);

	line = in.readLine();

	int row{ 0 };

	while (!line.isNull()) {

		auto ss = line.split(',');

		ss.removeAt(0);

		for (int i = 0; i < n_gene; ++i) {
			if (ss[i] == "1") {
				this->ppi_.mat_(row, i) = 1;
			}
		}

		line = in.readLine();

		++row;
	}

	return true;
};

bool ScentWorker::calculate_sr() {

	Eigen::MatrixXd adj = this->ppi_.mat_.cast<double>();

	int n_cell = this->expression_.cols();
	int n_gene = adj.rows();

	Eigen::ArrayXd srv = Eigen::ArrayXd::Zero(n_cell);

#pragma omp parallel for num_threads(10)
	for (int i = 0; i < n_cell; ++i) {
		Eigen::ArrayXd exp_v = this->expression_.mat_.col(i);
		Eigen::ArrayXd sumexp_v = adj * exp_v.matrix();
		Eigen::ArrayXd invp_v = exp_v * sumexp_v;
		double nf = invp_v.sum();
		invp_v /= nf;
		Eigen::MatrixXd pm = (adj.array().colwise() * exp_v).transpose().colwise() / sumexp_v;

		Eigen::ArrayXd sv = Eigen::ArrayXd::Zero(n_gene);
		for (int j = 0; j < n_gene; ++j) {
			for (int i = 0; i < n_gene; ++i) {
				if (pm(i, j) > 0.0) {
					sv[i] -= pm(i, j) * std::log(pm(i, j));
				}
			}
		}

		double sr = (invp_v * sv).sum();
		srv[i] = sr;
	}

	srv /= this->max_sr_;

	emit x_sr_ready(srv);

	return true;
};