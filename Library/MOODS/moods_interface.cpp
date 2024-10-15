#include "moods_interface.h"

#include "MOODS/moods.h"
#include "MOODS/moods_misc.h"
#include "MOODS/moods_tools.h"
#include "MOODS/match_types.h"
#include "MOODS/motif.h"
#include "MOODS/scanner.h"

#include "Custom.h"
#include "FileIO.h"
#include "GenomeUtility.h"

std::vector<std::vector<double> > mat_conversion(const Eigen::ArrayXXd& mat) {
	const Eigen::Index n_row = mat.rows(), n_col = mat.cols();
	score_matrix out(n_row);
	for (Eigen::Index i = 0; i < n_row; ++i) {
		out[i] = custom::cast<std::vector>(mat.row(i));
	};
	return out;
}


Eigen::ArrayXXd reverse_complement(const Eigen::ArrayXXd& mat) {
	const Eigen::Index n_row = mat.rows(), n_col = mat.cols();
	Eigen::ArrayXXd ret(n_row, n_col);
	for (Eigen::Index j = 0; j < n_col; ++j) {
		for (Eigen::Index i = 0; i < n_row; ++i) {
			ret(i, j) = mat(n_row - i - 1, n_col - j - 1);
		}
	}
	return ret;
}

double min_delta(const Eigen::ArrayXXd& mat) {
	const Eigen::Index n_row = mat.rows(), n_col = mat.cols();
	double ret = 0;

	double min_delta = std::numeric_limits<double>::infinity();

	for (Eigen::Index j = 0; j < n_col; ++j) {
		double max = -std::numeric_limits<double>::infinity(), second_max = max;
		for (Eigen::Index i = 0; i < n_row; ++i) {
			double value = mat(i, j);
			if (value > max) {
				second_max = max;
				max = value;
			}
			else if (value > second_max) {
				second_max = value;
			}
		}
		min_delta = std::min(min_delta, max - second_max);
	}
	return min_delta;
}

double min_score(const Eigen::ArrayXXd& mat) {
	const Eigen::Index n_row = mat.rows(), n_col = mat.cols();
	double ret = 0;

	for (Eigen::Index j = 0; j < n_col; ++j) {
		double min = std::numeric_limits<double>::infinity();
		for (Eigen::Index i = 0; i < n_row; ++i) {
			min = std::min(mat(i, j), min);
		}
		ret += min;
	}
	return ret;
}

double max_score(const Eigen::ArrayXXd& mat) {
	const Eigen::Index n_row = mat.rows(), n_col = mat.cols();
	double ret = 0;

	for (Eigen::Index j = 0; j < n_col; ++j) {
		double max = -std::numeric_limits<double>::infinity();
		for (Eigen::Index i = 0; i < n_row; ++i) {
			max = std::max(mat(i, j), max);
		}
		ret += max;
	}
	return ret;
}

GenomicRange match_motif(
	const QString& seq_name,
	PatternWeightMatrix pwm,
	const std::string& sequence,
	const double p_cutoff,
	const int window_size
) {
	auto background_frequency = custom::get_nucleic_acid_frequency(sequence);
	std::vector<double> nuc_freqs = custom::cast<std::vector>(background_frequency);

	pwm.from_frequency_to_log2();
	pwm.convert(background_frequency);

	std::vector<double> thresholds(2);
	std::vector<score_matrix> matrices(2);

	matrices[0] = mat_conversion(pwm.weight_.mat_);
	matrices[1] = MOODS::tools::reverse_complement(matrices[0]);
	thresholds[0] = MOODS::tools::threshold_from_p(matrices[0], nuc_freqs, p_cutoff);
	thresholds[1] = thresholds[0];

	MOODS::scan::Scanner scanner = MOODS::scan::Scanner(window_size);
	scanner.set_motifs(matrices, nuc_freqs, thresholds);

	auto results = scanner.scan(sequence);

	GenomicRange ret;

	int motif_width = pwm.weight_.cols();

	for (const auto& match : results[0]) {
		ret.append(seq_name, match.pos, motif_width, '+');
	}
	for (const auto& match : results[1]) {
		ret.append(seq_name, match.pos, motif_width, '-');
	}

	return ret;
}

void match_motif(
	MotifPosition& mp,
	const QStringList& sequences,
	const Eigen::Array4d& background_frequency,
	const double p_cutoff,
	const int window_size
) {
	auto pwms = mp.motifs_;
	for (auto& [name, pwm] : pwms) {
		pwm.from_frequency_to_log2();
		pwm.convert(background_frequency);
	}

	std::vector<double> nuc_freqs = custom::cast<std::vector>(background_frequency);
	const QVector<PatternWeightMatrix> mats = custom::values(pwms);

	auto peak_sequences = custom::sapply(sequences, [](const QString& str) {return str.toStdString(); });
	const size_t n_motif = pwms.size();
	std::vector<double> thresholds(2 * n_motif);
	std::vector<score_matrix> matrices(2 * n_motif);

	for (size_t i = 0; i < n_motif; i++) {
		matrices[i] = mat_conversion(mats[i].weight_.mat_);
		matrices[n_motif + i] = MOODS::tools::reverse_complement(matrices[i]);
		thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p_cutoff);
		thresholds[n_motif + i] = thresholds[i];
	}

	MOODS::scan::Scanner scanner = MOODS::scan::Scanner(window_size);
	scanner.set_motifs(matrices, nuc_freqs, thresholds);

	const std::size_t n_peak = peak_sequences.size();

	std::vector<Eigen::Triplet<std::size_t>> triplets;

#pragma omp parallel for
	for (int i = 0; i < n_peak; ++i) {

		if (peak_sequences[i].empty()) {
			continue;
		}

		auto results = scanner.scan(peak_sequences[i]);

	#pragma omp critical
		{
			for (std::size_t j = 0; j < n_motif; ++j) {
				for (const auto& match : results[j]) {
					mp.append(i, j, match.pos, match.score, '+');
				}
				for (const auto& match : results[n_motif + j]) {
					mp.append(i, j, match.pos, match.score, '-');
				}
			}
		}
	}
};