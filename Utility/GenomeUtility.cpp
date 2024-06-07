#include "GenomeUtility.h"

#include <QFileDialog>
#include <zlib.h>

#include "TwoBitFileProcessor.h"
#include "Custom.h"

std::vector<std::string> get_hexamers() {
	constexpr char a[4] = { 'A', 'C', 'G', 'T' };
	const std::string original_oligonucleotide{ "AAAAAA" };
	std::vector<std::string> hexamers;
	hexamers.reserve(4096);
	for (int i1 = 0; i1 < 4; ++i1) {
		for (int i2 = 0; i2 < 4; ++i2) {
			for (int i3 = 0; i3 < 4; ++i3) {
				for (int i4 = 0; i4 < 4; ++i4) {
					for (int i5 = 0; i5 < 4; ++i5) {
						for (int i6 = 0; i6 < 4; ++i6) {
							std::string this_hex = original_oligonucleotide;
							this_hex[0] = a[i1];
							this_hex[1] = a[i2];
							this_hex[2] = a[i3];
							this_hex[3] = a[i4];
							this_hex[4] = a[i5];
							this_hex[5] = a[i6];
							hexamers.push_back(this_hex);
						}
					}
				}
			}
		}
	}
	return hexamers;
}

namespace custom {

	GenomicRange get_hg38_gene_location() {

		QFile file(FILE_HUMAN_GENE_LOCATION);
		file.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream in(&file);

		QString line = in.readLine();

		GenomicRange ret;
		QStringList gene_names;

		int index = 0;
		while (!line.isNull()) {
			QStringList loc = line.split('\t');

			QString seq_name = _Cs standardize_chromosome_name(loc[0]);
			int start = loc[1].toInt();
			int end = loc[2].toInt();

			QString gene_name = loc[3];

			char strand = loc[4][0].toLatin1();

			ret.append(seq_name, start, end - start, strand);
			gene_names << gene_name;

			line = in.readLine();
		}

		ret.metadata_.update(METADATA_GENOMIC_RANGE_GENE_NAME, gene_names);

		return ret;
	};

	Eigen::Array4d get_nucleic_acid_frequency(const std::string& sequence) {

		Eigen::Array4d frequency = Eigen::Array4d::Zero(4);

		const std::size_t size = sequence.size();
		for (std::size_t i = 0; i < size; ++i) {
			char base = sequence[i];
			if (base == 'A') {
				++frequency[0];
			}
			else if (base == 'C') {
				++frequency[1];
			}
			else if (base == 'G') {
				++frequency[2];
			}
			else if (base == 'T') {
				++frequency[3];
			}
		}

		double letter_sum = frequency.sum();
		if (letter_sum > 0.0) {
			return frequency / letter_sum;
		}
		else return frequency;
	}

	Eigen::Array4d get_nucleic_acid_frequency(const QStringList& sequences) {

		Eigen::Array4d frequency = Eigen::Array4d::Zero(4);
		
		for (const auto& sequence : sequences) {
			const qsizetype size = sequence.size();
			for (qsizetype i = 0; i < size; ++i) {
				QChar base = sequence[i];
				if (base == 'A') {
					++frequency[0];
				}
				else if (base == 'C') {
					++frequency[1];
				}
				else if (base == 'G') {
					++frequency[2];
				}
				else if (base == 'T') {
					++frequency[3];
				}
			}
		}
		double letter_sum = frequency.sum();
		if (letter_sum > 0.0) {
			return frequency / letter_sum;
		}
		else return frequency;
	};

	GenomicRange stringlist_to_genomic_range(const QStringList& peak_names) {

		const qsizetype n_peak = peak_names.size();
		
		GenomicRange peak;

		if (n_peak == 0) {
			return peak;
		}

		QVector<int> peak_index;

		for (int i = 0; i < n_peak; ++i) {
			auto [sequence_name, start, end, success] = _Cs string_to_peak(peak_names[i]);

			if (success) {
				peak.append(sequence_name, start, end - start, '*');

				peak_index << i;
			}
		}

		peak.metadata_.update(METADATA_GENOMIC_RANGE_PEAK_INDEX, peak_index);

		peak.finalize();
		return peak;
	};

	std::tuple<QString, int, int, bool> string_to_peak(const QString& peak_name) {
		QStringList sub_string = _Cs multi_split(peak_name, QList<QChar>() << ':' << '-');
		qsizetype split_size = sub_string.size();
		if (split_size == 3) { 
			if (_Cs is_integer(sub_string[1]) && _Cs is_integer(sub_string[2])) {
				int start = sub_string[1].toInt();
				int end = sub_string[2].toInt();
				if (start >= 1 && start < end) {
					return std::make_tuple(_Cs standardize_chromosome_name(sub_string[0]), start, end, true);
				}
			}
		}
		return std::make_tuple("", 0, 0, false);
	};

	QString get_hexamer_from_loc(int loc) {
		if (loc > 4095 || loc < 0) {
			return QString();
		}
		char base[4] = {'A', 'C', 'G', 'T'};
		QString hexamer = "AAAAAA";
		for (int i = 0; i < 6; ++i) {
			hexamer[5 - i] = base[loc % 4];
			loc /= 4;
		}
	};

	std::tuple<QString, int, int, char, bool> find_gene_in_genome(const QString& gene_name, const GenomicRange& genome) {
		const auto& metadata = genome.metadata_;
		if (metadata.contains(METADATA_GENOMIC_RANGE_GENE_NAME) && metadata.data_type_.at(METADATA_GENOMIC_RANGE_GENE_NAME) != CustomMatrix::DataType::QString) {
			return std::make_tuple("", 0, 0, '-', false);
		}
		Eigen::ArrayX<bool> gene_filter = _Cs equal(metadata.get_const_qstring_reference(METADATA_GENOMIC_RANGE_GENE_NAME), gene_name);
		if (gene_filter.count() == 0) {
			return std::make_tuple("", 0, 0, '-', false);
		}
		int index = _Cs get_first_index(gene_filter);
		QString sequence_name = genome.sequence_names_[index];
		char strand = genome.strand_[index];
		int start = std::ranges::min(_Cs sliced(genome.ranges_.start_, gene_filter));
		int end = std::ranges::max(_Cs sliced(genome.get_sequence_end(), gene_filter));
		int width = end - start;

		return std::make_tuple(sequence_name, start, end, strand, true);
	};

	QStringList get_human_grch38_sequence(const GenomicRange& genomic_range) {
		TwoBitFileProcessor genome(FILE_HUMAN_GRCH38_2BIT);
		const qsizetype size = genomic_range.size();
		QStringList ret(size);
		for (qsizetype i = 0; i < size; ++i) {
			auto [sequence_name, start, end, strand] = genomic_range.at(i);
			ret[i] = genome.get_sequence(sequence_name, start, end);
		}
		return ret;
	};

	QString get_human_grch38_sequence(const QString& sequence_name, uint32_t start, uint32_t end) {
		TwoBitFileProcessor genome(FILE_HUMAN_GRCH38_2BIT);
		return genome.get_sequence(sequence_name, start, end);
	};

	std::vector<std::string> get_human_grch38_std_sequence(const GenomicRange& genomic_range) {
		TwoBitFileProcessor genome(FILE_HUMAN_GRCH38_2BIT);
		const qsizetype size = genomic_range.size();
		std::vector<std::string> ret(size);
		for (qsizetype i = 0; i < size; ++i) {
			auto [sequence_name, start, end, strand] = genomic_range.at(i);
			ret[i] = genome.get_std_sequence(sequence_name, start, end);
		}
		return ret;
	};

	std::string get_human_grch38_std_sequence(const QString& sequence_name, uint32_t start, uint32_t end) {
		TwoBitFileProcessor genome(FILE_HUMAN_GRCH38_2BIT);
		return genome.get_std_sequence(sequence_name, start, end);
	};

	std::pair<std::vector<std::string>, std::vector<std::size_t> > get_6_nucleotide_frequency(const std::string& sequence) {
		std::vector<std::string> hexamers = get_hexamers();
		const char* data = sequence.data();
		const std::size_t size = sequence.size();
		std::size_t count = 0;
		std::vector<std::size_t> frequency(4096, 0);
		for (std::size_t i = 0; i < size - 5; ++i) {
			int loc = get_hexamer_loc(data + i);
			if (loc > -1) {
				++frequency[loc];
				++count;
			}
		}
		return std::make_pair(hexamers, frequency);
	};

	std::pair<std::vector<std::string>, std::vector<std::size_t> > get_6_nucleotide_insertion_frequency(
		QString sequence_name,
		const std::string& sequence,
		const Fragments& fragments
	) {
		sequence_name = _Cs standardize_chromosome_name(sequence_name);
		if (!fragments.data_.contains(sequence_name)) {
			return std::pair<std::vector<std::string>, std::vector<std::size_t> >();
		}

		std::vector<int> locations;
		std::vector<std::string> hexamers = get_hexamers();
		std::vector<std::size_t> frequency(4096, 0);

		const std::size_t sequence_size = sequence.size();
		const std::size_t minimum_location = 4, maximum_location = sequence_size - 2;

		const char* data = sequence.data();

		const auto& chr = fragments.data_.at(sequence_name);

		const int n_cell = chr.size();
		for (int i = 0; i < n_cell; ++i) {
			const auto& cell_data = chr[i];
			const std::size_t n_fragments = cell_data.first.size();

			const auto& start_data = cell_data.first;

			for (auto loc : start_data) {
				if (loc >= minimum_location && loc <= maximum_location) {
					loc = get_hexamer_loc(data + loc - 4);
					if (loc > -1) {
						++frequency[loc];
					}
				}
			}

			const auto& end_data = cell_data.second;

			for (auto loc : end_data) {
				if (loc >= minimum_location && loc <= maximum_location) {
					loc = get_hexamer_loc(data + loc - 4);
					if (loc > -1) {
						++frequency[loc];
					}
				}
			}
		}

		return std::make_pair(hexamers, frequency);
	};

	std::tuple<bool, std::vector<std::string>, std::vector<std::size_t> > 
		get_6_nucleotide_insertion_frequency(
		std::string sequence_name,
		const std::string& sequence,
		const std::unordered_set<std::string>& barcodes,
		const QStringList& fragments_files) {

		sequence_name = _Cs standardize_chromosome_name(sequence_name);

		std::vector<int> locations;

		int buffer_length = 256;
		char buffer[256];
		std::string barcode;

		std::vector<std::string> hexamers = get_hexamers();
		std::vector<std::size_t> frequency(4096, 0);

		const std::size_t sequence_size = sequence.size();
		const std::size_t minimum_location = 4, maximum_location = sequence_size - 2;

		const char* data = sequence.data();

		for (const auto& fragments_file : fragments_files) {
			gzFile fragments = gzopen(fragments_file.toUtf8().data(), "rb");

			if (fragments == NULL) {
				return std::make_tuple(false, std::vector<std::string>(), std::vector<std::size_t>());
				continue;
			}

			bool not_end_of_file = false;
			while (not_end_of_file = gzgets(fragments, buffer, buffer_length) != 0)
			{
				if (buffer[0] != '#')break;
			}
			if (!not_end_of_file) {
				gzclose(fragments);
				return std::make_tuple(false, std::vector<std::string>(), std::vector<std::size_t>());
				continue;
			}
						
			do {
				const char* c = buffer, * d, * seq_end;
				while (*c != '\t') {
					++c;
				}
				seq_end = c;

				std::string reads_sequence_name = _Cs standardize_chromosome_name(std::string(buffer, seq_end - buffer));				

				if (sequence_name != reads_sequence_name) {
					continue;
				}

				int start = _Cs atoi_specialized(&c);
				int end = _Cs atoi_specialized(&c);	

				++c;
				d = c;
				while (*c != '\t') {
					++c;
				}
				++c;
				barcode.assign(d, c - d - 1);

				if (barcodes.contains(barcode)) {

					int loc = start;

					if (loc >= minimum_location && loc <= maximum_location) {
						loc = get_hexamer_loc(data + loc - 4);
						if (loc > -1) {
							++frequency[loc];
						}
					}

					loc = end;

					if (loc >= minimum_location && loc <= maximum_location) {
						loc = get_hexamer_loc(data + loc - 4);
						if (loc > -1) {
							++frequency[loc];
						}
					}
				}

			} while (not_end_of_file = gzgets(fragments, buffer, buffer_length) != 0);

			gzclose(fragments);
		}

		return std::make_tuple(true, hexamers, frequency);
	};

	QStringList get_available_human_genome() {
		QDir d(qApp->applicationDirPath() + "/Resources/Genome/Human");
		QStringList filenames = d.entryList();
		QStringList human_genome;
		for (const auto& name : filenames) {
			if (name.endsWith(".2bit")) {
				human_genome << name;
			}
		}
		return human_genome;
	};
};