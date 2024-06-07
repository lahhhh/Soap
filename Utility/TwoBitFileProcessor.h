#pragma once

#include "Identifier.h"

#include <QString>
#include <QMap>
#include <QVector>
#include <QFile>

#include "GenomicRange.h"

/*

	Example:

	TwoBitFileProcessor genome;
	genome.set_sequence_file(FILE_HUMAN_GRCH38_2BIT);
	auto seq = genome.get_sequence(...);
	...

	or 

	TwoBitFileProcessor genome(FILE_HUMAN_GRCH38_2BIT);
	auto seq = genome.get_sequence(...);
	...
*/

// internal
class TwoBitFileProcessor
{
public:

	TwoBitFileProcessor() = default;

	explicit TwoBitFileProcessor(const QString& file_name);

	QString file_name_;

	QFile file_;

	bool byte_swapped_ = false;
	bool loaded_ = false;

	qsizetype sequence_count_ = 0;

	uint64_t file_size_ = 0;

	QMap<QString, uint32_t> sequence_offset_;
	QMap<QString, uint32_t> sequence_offset2_;
	QMap<QString, uint32_t> sequence_dna_size_;
	QMap<QString, uint32_t> sequence_packed_size_;
	QMap<QString, uint32_t> sequence_bytes_;

	QMap<QString, QVector<uint32_t> > sequence_n_block_starts_;
	QMap<QString, QVector<uint32_t> > sequence_n_block_sizes_;
	QMap<QString, QVector<uint32_t> > sequence_mask_block_starts_;
	QMap<QString, QVector<uint32_t> > sequence_mask_block_sizes_;

	void set_sequence_file(const QString& sequenc_file_name);

	bool load_sequence();

	void load_index(QDataStream& ds);

	void release();

	void initialize_sequence();

	//[start, end)
	QString get_sequence(const QString& sequence_name, uint32_t start, uint32_t end);

	QStringList get_sequence(const GenomicRange& genomic_range);

	//[start, end)
	std::string get_std_sequence(const QString& sequence_name, uint32_t start, uint32_t end) ;

	std::vector<std::string> get_std_sequence(const GenomicRange& genomic_range);

	std::string get_std_whole_sequence(const QString& sequence_name);
};

