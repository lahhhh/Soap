#pragma once

#include "Identifier.h"

#include <QStringList>
#include <QMap>

#include "PatternWeightMatrix.h"
#include "IRange.h"
#include "GenomicRange.h"
#include "Footprint.h"
#include "ChromVAR.h"


class MotifPosition
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(MotifPosition);

	MotifPosition(const std::size_t n_peak, const std::size_t n_motif);

	G_SET_IDENTIFIER("MotifPosition");

	struct match {
		std::size_t peak_index;
		std::size_t position;
		float score;
		char strand;
	};

	DataType data_type_{ DataType::Plain };

	GenomicRange peak_locations_;

	std::map<QString, PatternWeightMatrix> motifs_;

	// length == motif number; this is not ordered
	std::vector<std::vector<match>> motif_information_;

	QStringList peak_names_;
	QStringList motif_names_;

	SOAP_SUBMODULES(Footprint);
	SOAP_SUBMODULES(ChromVAR);

	qsizetype n_motif() const;
	qsizetype n_peak() const;

	Eigen::ArrayXi get_match_count(const int motif_index) const;
	int get_match_count(const int peak_index, const int motif_index) const;
	int get_match_count(const QVector<int>& peak_index, const int motif_index) const;
	int get_match_count(const QString& peak_name, const QString& motif_name) const;

	bool contains(const QString& motif_name) const;

	GenomicRange get_position(const QString& motif_name) const;

	QStringList get_match_peak_position(const QString& motif_name) const;

	MotifPosition& append(
		const std::size_t peak_location,
		const std::size_t motif_location,
		const std::size_t position,
		const float score,
		const char strand
	);

	Eigen::MatrixXi get_motif_matrix() const;

	QStringList peak_to_motif(int peak_index) const;

	QStringList peak_to_tf(int peak_index) const;

	QVector<int> get_motif_match_times(const QVector<int>& peak_index) const;
};

