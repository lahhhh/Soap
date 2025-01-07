#pragma once

#include "Identifier.h"

#include <unordered_set>

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"
#include "TranscriptModel.h"


/*
* modified from velocyto package
* La Manno, Gioele et al. ¡°RNA velocity of single cells.¡± 
* Nature vol. 560,7719 (2018): 494-498. doi:10.1038/s41586-018-0414-6
*/

class VelocytoWorker
	: public QObject
{
	Q_OBJECT
public:
	VelocytoWorker(
		const SingleCellRna* single_cell_rna,
		const QStringList& bam_files
		)
		:
		single_cell_rna_(single_cell_rna),
		bam_files_(bam_files),
		mode_(WorkMode::SingleCellRna)
	{
		this->barcode_tag_[0] = 'C';
		this->barcode_tag_[1] = 'B';

		this->umi_tag_[0] = 'U';
		this->umi_tag_[1] = 'B';
	}

	VelocytoWorker(
		const SingleCellMultiome* single_cell_multiome,
		const QStringList& bam_files
	)
		:
		single_cell_multiome_(single_cell_multiome),
		bam_files_(bam_files),
		mode_(WorkMode::SingleCellMultiome)
	{
		this->barcode_tag_[0] = 'C';
		this->barcode_tag_[1] = 'B';

		this->umi_tag_[0] = 'U';
		this->umi_tag_[1] = 'B';
	}

	const SingleCellRna* single_cell_rna_ = nullptr;

	const SingleCellMultiome* single_cell_multiome_ = nullptr;

	enum class WorkMode{SingleCellRna, SingleCellMultiome};

	WorkMode mode_ = WorkMode::SingleCellMultiome;	

	QStringList bam_files_;
	QStringList gene_names_;
	QStringList cell_names_;

	// no need to care about strand --> from the same RNA molecule
	std::map<std::string, std::vector<TranscriptModel> > transcript_models_;

	std::map<std::string, std::vector<int> > mask_intervals_;

	std::vector<std::unordered_set<int> > umi_set_;

	std::unordered_map<std::string, int> barcode_index_;

	std::unordered_map<std::string, int> gene_name_index_;

	Eigen::MatrixXi spliced_counts_;

	Eigen::MatrixXi unspliced_counts_;

	char barcode_tag_[2];

	char umi_tag_[2];

	std::unique_ptr<VelocytoBase> res_{ nullptr };

	void create_index();

	bool parse_gtf(const char* file_name);

	bool parse_mask(const char* file_name);

	bool filter_transcript_model();

	bool filter_mask_intervals();

	bool count_reads();

	// return 0 : mapped, 1 : found no umi/barcode, 2 : no segment, 3 : masked, 4 : no barcode, 5 : no match, 
	// 6 : no gene, 7 : umi dup, 8 : multi match
	int map_read(const char* read, int alignment_size, const std::string& ref_name);

	bool interval_masked(const std::string& chromosome_name, int32_t start, int32_t end);

	// return value - 0 : no map; 1 : unspliced; 2 : spliced
	int count(const std::string& chromosome_name, const std::vector<int>& segments, std::string* ptr_gene_name);

	// assume that umi length is less than 17bp
	bool get_umi_and_barcode(const char* read, int alignment_size, int* umi, std::string* barcode);

	void generate_count_matrix();

public:

	bool work();

public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_velocyto_ready(VelocytoBase* velocyto_base);	

};

