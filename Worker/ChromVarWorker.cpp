#include "ChromVARWorker.h"

#include "pval.h"

bool ChromVARWorker::work() {

	if (this->species_ != soap::Species::Human) {
		G_TASK_WARN("ChromVAR now only support human genome.");
		return false;
	}

	this->genome_.set_sequence_file(FILE_HUMAN_GRCH38_2BIT);

	if (!this->validate_input()) {
		return false;
	}

	if (!this->add_gc_bias()) {
		return false;
	}

	if (!this->get_background_peaks()) {
		return false;
	}

	if (!this->compute_deviations()) {
		return false;
	}

	return true;
};

void ChromVARWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	emit x_chromvar_ready(this->res_.release());

	G_TASK_END;
};

bool ChromVARWorker::validate_input() {
	
	if (this->motif_position_->peak_locations_.size() != this->atac_counts_->rows()) {
		G_TASK_WARN("Unmatched Motif Location and ATAC data.");
		return false;
	}

	auto fragments_per_peak = custom::row_sum(this->atac_counts_->mat_);
	if (custom::any(fragments_per_peak, 0)) {
		G_TASK_WARN("Detect peak with no fragments.");
		return false;
	}

	this->expectations_ = fragments_per_peak.cast<double>() / this->atac_counts_->mat_.sum();

	return true;
};

bool ChromVARWorker::add_gc_bias() {

	this->peak_sequences_ = this->genome_.get_std_sequence(this->motif_position_->peak_locations_);

	auto n_invalid_peak = std::ranges::count_if(this->peak_sequences_, [](auto&& s) {return s.empty(); });

	if (n_invalid_peak > this->peak_sequences_.size() / 2) {
		G_TASK_WARN("Find two many invalid peak locations.");
		return false;
	}

	int n_peak = this->peak_sequences_.size();
	this->gc_.resize(n_peak);
	this->gc_.setZero();

	for (int i = 0; i < n_peak; ++i) {

		if (this->peak_sequences_[i].empty()) {
			continue;
		}

		int a_count{ 0 };
		int c_count{ 0 };
		int g_count{ 0 };
		int t_count{ 0 };

		for (auto&& c : this->peak_sequences_[i]) {
			switch (c)
			{
			case 'A':
				++a_count;
				break;
			case 'C':
				++c_count;
				break;
			case 'G':
				++g_count;
				break;
			case 'T':
				++t_count;
				break;
			default:
				break;
			}
		}

		this->gc_[i] = (g_count + c_count) / (double)(a_count + g_count + c_count + t_count);
	}

	return true;
};

Eigen::ArrayXi prob_sample_replace(int seed, int size, Eigen::ArrayXd prob) {

	int n_orig = prob.size();
	Eigen::ArrayXi index(size);

	Eigen::ArrayXd HL_dat(n_orig);
	Eigen::ArrayXd alias_tab(n_orig);

	auto h = HL_dat.begin();
	auto h0 = HL_dat.begin();

	auto l = HL_dat.end();
	auto l0 = HL_dat.end();

	for (int i = 0; i < n_orig; ++i) {
		prob[i] *= n_orig;

		if (prob[i] < 1.0) {
			*(h++) = i;
		}
		else {
			*(--l) = i;
		}
	}

	if ((h > h0) && (l < l0)) {
		for (int k = 0; k < n_orig; ++k) {
			int ii = HL_dat[k];
			int jj = *l;
			alias_tab[ii] = jj;
			prob[jj] += (prob[ii] - 1.0);
			if (prob[jj] < 1.0) {
				++l;
			}
			if (l == l0) {
				break;
			}
		}
	}

	for (int i = 0; i < n_orig; ++i) {
		prob[i] += i;
	}

	std::uniform_real_distribution d(0.0, 1.0);
	std::default_random_engine e;
	e.seed(seed);

	for (int i = 0; i < size; ++i) {
		double r_u = d(e) * n_orig;

		int kk = (int)r_u;
		if (kk == n_orig) {
			--kk;
		}
		index[i] = (r_u < prob[kk]) ? kk : alias_tab[kk];
	}

	return index;
}

bool ChromVARWorker::get_background_peaks() {

	auto fragments_per_peak = custom::row_sum(this->atac_counts_->mat_);

	Eigen::ArrayXd intensity = log10(fragments_per_peak.cast<double>());

	int n_peak = intensity.size();

	Eigen::MatrixXd norm_mat(n_peak, 2);

	norm_mat << intensity, this->gc_;

	Eigen::LLT<Eigen::MatrixXd> llt(custom::cov(norm_mat));
	if (llt.info() != Eigen::Success) {
		G_TASK_WARN("LLT decomposition failed!");
		return false;
	}
	Eigen::MatrixXd chol_cov_mat = llt.matrixU();

	Eigen::LLT<Eigen::MatrixXd> llt2(chol_cov_mat.transpose());
	if (llt2.info() != Eigen::Success) {
		G_TASK_WARN("LLT decomposition failed!");
		return false;
	}

	Eigen::MatrixXd trans_norm_mat = llt2.solve(norm_mat.transpose());
	trans_norm_mat.transposeInPlace();
	constexpr int bs = 50;

	Eigen::ArrayXd bin1 = Eigen::ArrayXd::LinSpaced(bs, trans_norm_mat.col(0).minCoeff(), trans_norm_mat.col(0).maxCoeff());
	Eigen::ArrayXd bin2 = Eigen::ArrayXd::LinSpaced(bs, trans_norm_mat.col(1).minCoeff(), trans_norm_mat.col(1).maxCoeff());

	Eigen::MatrixXd bin_data(bs * bs, 2);
	for (int i = 0; i < bs; ++i) {
		bin_data.col(0).segment(i * bs, bs).setConstant(bin1[i]);
		bin_data.col(1).segment(i * bs, bs) = bin2;
	}

	auto bin_dist = custom::euclidean_distance_mt(bin_data, false);
	constexpr double w = 0.1;
	auto bin_p = dnorm(bin_dist, 0, w);
	auto bin_membership = custom::find_nearest(trans_norm_mat, bin_data, false);
	Eigen::ArrayXd bin_density = Eigen::ArrayXd::Zero(bs * bs);
	std::ranges::for_each(bin_membership, [&bin_density](auto t) {++bin_density[t]; });
	constexpr int n_iteration = 50;
	this->background_peaks_.resize(n_peak, n_iteration);

#pragma omp parallel for
	for (int i = 0; i < bs * bs; ++i) {
		auto ix = custom::which(bin_membership == i);
		int n_ix_ele = ix.size();

		if (n_ix_ele == 0) {
			continue;
		}

		Eigen::ArrayXd p_tmp = bin_p.col(i);
		Eigen::ArrayXd p = p_tmp(bin_membership) / bin_density(bin_membership);
		p /= p.sum();
		auto sampled = prob_sample_replace(i, n_iteration * n_ix_ele, p);
		for (int j = 0; j < n_ix_ele; ++j) {
			this->background_peaks_.row(ix[j]) = sampled.segment(j * n_iteration, n_iteration);
		}
	}

	return true;
};

bool ChromVARWorker::compute_deviations() {

	auto motif_match = this->motif_position_->get_motif_matrix();

	Eigen::ArrayXd fragments_per_sample = custom::col_sum(this->atac_counts_->mat_).cast<double>();

	int n_peak = motif_match.rows();
	int n_motif = motif_match.cols();
	int n_cell = this->atac_counts_->cols();

	this->res_.reset(new ChromVAR());
	this->res_->motif_names_ = this->motif_position_->motif_names_;
	this->res_->z_ = Eigen::MatrixXd::Zero(n_motif, n_cell);
	this->res_->dev_ = Eigen::MatrixXd::Zero(n_motif, n_cell);

#pragma omp parallel for
	for (int i = 0; i < n_motif; ++i) {

		auto peak_index = custom::which(motif_match.col(i).array() > 0);

		if (peak_index.isEmpty()) {
			continue;
		}

		int n_motif_peak = peak_index.size();

		Eigen::MatrixXd sampled_deviation;
		Eigen::ArrayXd observed_deviation, expected;
		if (n_motif_peak == 1) {
			
			Eigen::ArrayXd observed = this->atac_counts_->mat_.row(peak_index[0]).cast<double>();
			expected = this->expectations_[peak_index[0]] * fragments_per_sample;
			observed_deviation = (observed - expected) / expected;

			Eigen::MatrixXd sampled = custom::row_reordered(this->atac_counts_->mat_, this->background_peaks_.row(peak_index[0])).toDense().cast<double>();
			Eigen::MatrixXd sampled_expected = this->expectations_(this->background_peaks_.row(peak_index[0])).matrix() * fragments_per_sample.matrix().transpose();
			sampled_deviation = (sampled - sampled_expected).array() / sampled_expected.array();
		}
		else {
			Eigen::ArrayXd observed = custom::row_reorder_and_column_sum(this->atac_counts_->mat_, peak_index).cast<double>();

			expected = this->expectations_(peak_index).sum() * fragments_per_sample;
			observed_deviation = (observed - expected) / expected;

			int n_iterations = this->background_peaks_.cols();

			Eigen::MatrixXi bgp = this->background_peaks_(peak_index, Eigen::all);

			Eigen::SparseMatrix<int> sample_mat;
			{
				std::vector<Eigen::Triplet<int>> triplets;

				for (auto p : peak_index) {
					for (int j = 0; j < n_iterations; ++j) {
						triplets.emplace_back(j, this->background_peaks_(p, j), 1);
					}
				}

				sample_mat.resize(n_iterations, n_peak);
				sample_mat.setFromTriplets(triplets.cbegin(), triplets.cend());
			}
			Eigen::MatrixXd sampled = (sample_mat * this->atac_counts_->mat_).toDense().cast<double>();
			Eigen::MatrixXd sample_expected(n_iterations, n_cell);

			for (int j = 0; j < n_iterations; ++j) {
				sample_expected.row(j) = this->expectations_(bgp.col(j)).sum() * fragments_per_sample;
			}
			sampled_deviation = (sampled - sample_expected).array() / sample_expected.array();
		}

		constexpr double threshold = 1.0;
		auto fail_filter = custom::which(expected < threshold);

		Eigen::ArrayXd mean_sampled_deviation = sampled_deviation.colwise().mean();

		Eigen::ArrayXd sd_sampled_deviation(n_cell);
		for (int i = 0; i < n_cell; ++i) {
			sd_sampled_deviation[i] = custom::sd(sampled_deviation.col(i));
		}

		Eigen::ArrayXd normdev = observed_deviation - mean_sampled_deviation;
		Eigen::ArrayXd z = normdev / sd_sampled_deviation;

		if (!fail_filter.isEmpty()) {
			custom::assign(normdev, 0.0 /*std::nan("0")*/, fail_filter);
			custom::assign(z, 0.0 /*std::nan("0")*/, fail_filter);
		}

		this->res_->z_.row(i) = z;
		this->res_->dev_.row(i) = normdev;
	}

	return true;
};