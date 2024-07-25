#pragma once

#ifndef _Cs
#define _Cs ::custom::
#endif // !_Cs

#include "Identifier.h"

#include <unordered_map>
#include <QStringList>
#include <unordered_set>

#include "GenomicRange.h"
#include "Fragments.h"

namespace custom {

	// A = 0, C = 1, G = 2, T = 3
	int get_hexamer_loc(const char str[6]);

	QString get_hexamer_from_loc(int loc);

	// return start, end 
	std::tuple<QString, int, int, bool> string_to_peak(const QString& peak_name);
	
	GenomicRange stringlist_to_genomic_range(const QStringList& peak_names);

	std::pair<std::vector<std::string>, std::vector<std::size_t> > get_6_nucleotide_frequency(const std::string& sequence);

	std::tuple<bool, std::vector<std::string>, std::vector<std::size_t> > get_6_nucleotide_insertion_frequency(
		std::string sequence_name, 
		const std::string& sequence,
		const std::unordered_set<std::string>& barcodes,
		const QStringList& fragments_files);

	std::pair<std::vector<std::string>, std::vector<std::size_t> > get_6_nucleotide_insertion_frequency(
		QString sequence_name,
		const std::string& sequence,
		const Fragments& fragments
		);

	QStringList get_available_human_genome();

	std::tuple<QString, int, int, char, bool> find_gene_in_genome(const QString& gene_name, const GenomicRange& genome);

	// use when only read once in a job
	QStringList get_human_grch38_sequence(const GenomicRange& genomic_range);

	// use when only read once in a job
	QString get_human_grch38_sequence(const QString& sequence_name, uint32_t start, uint32_t end);

	// use when only read once in a job
	std::vector<std::string> get_human_grch38_std_sequence(const GenomicRange & genomic_range);

	// use when only read once in a job
	std::string get_human_grch38_std_sequence(const QString& sequence_name, uint32_t start, uint32_t end);

	Eigen::Array4d get_nucleic_acid_frequency(const QStringList& sequences);
	Eigen::Array4d get_nucleic_acid_frequency(const std::string& sequence);

	GenomicRange get_hg38_gene_location();
};

