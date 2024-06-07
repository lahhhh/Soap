#include "CalculateGeneActivityWorker.h"

#include <zlib.h>

int CalculateGeneActivityWorker::create_index_human() {
	QFile file(FILE_HUMAN_GENE_LOCATION);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	QTextStream in(&file);

	QString line = in.readLine();

	int gene_index = 0;
	while (!line.isNull()) {
		QStringList loc = line.split('\t');

		QString seq_name = _Cs standardize_chromosome_name(loc[0]);
		int start = loc[1].toInt();
		int end = loc[2].toInt();

		if (this->use_tss_) {
			if (loc[4] == "-") {
				start = end - 500;
				end += 500;

				if (start < 1) {
					start = 1;
				}
			}
			else {
				end = start + 500;
				start -= 500;

				if (start < 1) {
					start = 1;
				}
			}
		}

		this->genome_.append_index(seq_name, start, end, gene_index++);
		this->gene_names_ << loc[3];

		line = in.readLine();
	}

	this->genome_.reorganize();

	return gene_index;
};

bool CalculateGeneActivityWorker::create_index() {

	int n_cell = this->fragments_->cell_names_.size();

	if (this->species_ != soap::Species::Human) {
		G_TASK_WARN("Gene Activity Calculation now only support human.");
		return false;
	}

	int n_gene = this->create_index_human();

	this->dense_counts_.resize(n_gene, n_cell);

	return true;
};

void CalculateGeneActivityWorker::find_row(const QString& seq_name, int cell_loc, int start, int end) {

	auto [r1, r2] = this->genome_.find_location(seq_name, start, end);

	if (this->genome_.success(r1)) {
		++this->dense_counts_(r1, cell_loc);
	}

	if (this->genome_.success(r2)) {
		++this->dense_counts_(r2, cell_loc);
	}

};

bool CalculateGeneActivityWorker::calculate_counts() {
	for (const auto& [name, data] : this->fragments_->data_) {

		const int n_cell = data.size();

		for (int i = 0; i < n_cell; ++i) {
			const auto& cell_data = data[i];
			const std::size_t fragments_size = cell_data.first.size();
			for (size_t j = 0; j < fragments_size; ++j) {
				this->find_row(name, i, cell_data.first[j], cell_data.second[j]);
			}
		}
	}
	return true;
}

void CalculateGeneActivityWorker::run() {

	G_TASK_LOG("Creating index...");

	this->create_index();

	G_TASK_LOG("Mapping fragments...");

	if (!this->calculate_counts()) {
		G_TASK_END;
	}


	G_TASK_LOG("Building counts matrix...");

	this->build_matrix();

	emit x_gene_activity_ready(this->counts_);

	G_TASK_LOG("Gene Activity calculation finished.");

	G_TASK_END;
}

void CalculateGeneActivityWorker::build_matrix() {

	this->counts_ = new SparseInt();

	this->counts_->rownames_ = this->gene_names_;

	this->counts_->colnames_ = this->fragments_->cell_names_;

	Eigen::SparseMatrix<int>& counts_matrix = this->counts_->mat_;

	constexpr int min_peak_umi_count = 5;

	auto filtered = _Cs which(this->dense_counts_.array().rowwise().sum() > 5);

	counts_matrix = this->dense_counts_(filtered, Eigen::all).sparseView();
	this->counts_->rownames_ = _Cs reordered(this->counts_->rownames_, filtered);

	this->dense_counts_.resize(0, 0);
}