#pragma once

#include "Identifier.h"

#include "SparseDouble.h"
#include "Cnv.h"

/*
* modified from sciCNV package
* Mahdipour-Shirayeh, Ali et al. ¡°sciCNV: high-throughput paired profiling of transcriptomes and DNA copy number variations at single-cell resolution.¡± 
* Briefings in bioinformatics vol. 23,1 (2022): bbab413. doi:10.1093/bib/bbab413
*/

class ScicnvWorker 
	: public QObject
{
	Q_OBJECT
public:

	ScicnvWorker(
		const SparseDouble& mat,
		const QStringList& metadata,
		const QStringList& reference = {},
		soap::Species species = soap::Species::Human,
		double threshold = 0.4,
		double sharpness = 1.0
	) :
		mat_(mat),
		metadata_(metadata),
		reference_(reference),
		species_(species),
		threshold_(threshold),
		sharpness_(sharpness),
		n_cell_(mat.cols())
	{
		QString current_factor = this->metadata_.front();
		int current_index = 0;
		for (int i = 1; i < n_cell_; ++i)
		{
			if (this->metadata_[i] == current_factor) {
				continue;
			}
			else {
				this->metadata_location_.emplace_back(current_factor, current_index, i - current_index);
				current_factor = this->metadata_[i];
				current_index = i;
			}
		};

		this->metadata_location_.emplace_back(current_factor, current_index, this->n_cell_ - current_index);
	}

	SparseDouble mat_;

	QStringList metadata_;
	QStringList reference_;

	soap::Species species_;

	double threshold_;

	double sharpness_;

	int n_cell_{ 0 };
	int n_gene_{ 0 };
	int resolution_{ 0 };

	std::vector<std::tuple<QString, int, int> > chromosome_location_;
	std::vector<std::tuple<QString, int, int> > metadata_location_;

	Eigen::ArrayXXd moving_average_;
	Eigen::ArrayXXd W_;
	Eigen::ArrayXXd U_;
	Eigen::ArrayXXd V_;

	Eigen::ArrayXd reference_expression_;
	Eigen::ArrayXd reference_moving_average_;

	bool sort_expression_by_chromosome();

	bool filter_data();

	void compute_reference_moving_average(bool sharp);

	void compute_reference();

	void compute_moving_average(bool sharp);

	void compute_W();

	void compute_U();

	void compute_V();

	void compute_gradient();

	void compute_cnv();

	static Eigen::ArrayXd digitalize(const Eigen::ArrayXd& arr);

	static Eigen::ArrayXd slope(const Eigen::ArrayXXd& arr);

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_cnv_ready(CNV*);
	
};
