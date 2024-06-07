#pragma once

#include "PatternWeightMatrix.h"
#include "MotifPosition.h"

void match_motif(
	MotifPosition& mp,
	const QStringList& sequences,
	const Eigen::Array4d& background_frequency,
	const double p_cutoff,
	const int window_size
);

GenomicRange match_motif(
	const QString& seq_name,
	PatternWeightMatrix pwm,
	const std::string& sequence,
	const double p_cutoff,
	const int window_size
);


