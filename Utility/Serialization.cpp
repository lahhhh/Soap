#include "Serialization.h"

void sread(std::ifstream& in, GraphSettings& val, const QString& edition) {
	sread(in, val.string_info_, edition);
	sread(in, val.string_list_info_, edition);
	sread(in, val.font_info_, edition);
	sread(in, val.integer_info_, edition);
	sread(in, val.bool_info_, edition);
	sread(in, val.color_info_, edition);
	sread(in, val.palette_, edition);
};

void swrite(std::ofstream& out, const GraphSettings& val, const QString& edition) {
	swrite(out, val.string_info_, edition);
	swrite(out, val.string_list_info_, edition);
	swrite(out, val.font_info_, edition);
	swrite(out, val.integer_info_, edition);
	swrite(out, val.bool_info_, edition);
	swrite(out, val.color_info_, edition);
	swrite(out, val.palette_, edition);
};

void sread(std::ifstream& in, SingleCellMultiome& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.species_, edition);
	sread(in, val.random_state_, edition);

	sread(in, val.string_information_, edition);
	sread(in, val.integer_information_, edition);
	sread(in, val.double_information_, edition);

	sread(in, val.string_vectors_, edition);
	sread(in, val.integer_vectors_, edition);
	sread(in, val.double_vectors_, edition);

	sread(in, SUBMODULES(val, DataField), edition);
	sread(in, SUBMODULES(val, Metadata), edition);
	sread(in, SUBMODULES(val, GSEA), edition);
	sread(in, SUBMODULES(val, CellChat), edition);
	sread(in, SUBMODULES(val, CNV), edition);
	sread(in, SUBMODULES(val, GenomicRange), edition);
	sread(in, SUBMODULES(val, MotifPosition), edition);
	sread(in, SUBMODULES(val, CoverageTrack), edition);
	sread(in, SUBMODULES(val, Fragments), edition);
	sread(in, SUBMODULES(val, VelocytoBase), edition);
	sread(in, SUBMODULES(val, Pando), edition);
	sread(in, SUBMODULES(val, Monocle3), edition);
	sread(in, SUBMODULES(val, Cicero), edition);
};

void swrite(std::ofstream& out, const SingleCellMultiome& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.species_, edition);
	swrite(out, val.random_state_, edition);

	swrite(out, val.string_information_, edition);
	swrite(out, val.integer_information_, edition);
	swrite(out, val.double_information_, edition);

	swrite(out, val.string_vectors_, edition);
	swrite(out, val.integer_vectors_, edition);
	swrite(out, val.double_vectors_, edition);
	swrite(out, SUBMODULES(val, DataField), edition);
	swrite(out, SUBMODULES(val, Metadata), edition);
	swrite(out, SUBMODULES(val, GSEA), edition);
	swrite(out, SUBMODULES(val, CellChat), edition);
	swrite(out, SUBMODULES(val, CNV), edition);
	swrite(out, SUBMODULES(val, GenomicRange), edition);
	swrite(out, SUBMODULES(val, MotifPosition), edition);
	swrite(out, SUBMODULES(val, CoverageTrack), edition);
	swrite(out, SUBMODULES(val, Fragments), edition);
	swrite(out, SUBMODULES(val, VelocytoBase), edition);
	swrite(out, SUBMODULES(val, Pando), edition);
	swrite(out, SUBMODULES(val, Monocle3), edition);
	swrite(out, SUBMODULES(val, Cicero), edition);
};

void sread(std::ifstream& in, DataField& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, SUBMODULES(val, SparseInt), edition);
	sread(in, SUBMODULES(val, SparseDouble), edition);
	sread(in, SUBMODULES(val, Embedding), edition);
	sread(in, SUBMODULES(val, DifferentialAnalysis), edition);
};

void swrite(std::ofstream& out, const DataField& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, SUBMODULES(val, SparseInt), edition);
	swrite(out, SUBMODULES(val, SparseDouble), edition);
	swrite(out, SUBMODULES(val, Embedding), edition);
	swrite(out, SUBMODULES(val, DifferentialAnalysis), edition);
};

void sread(std::ifstream& in, SingleCellAtac& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.species_, edition);
	sread(in, val.random_state_, edition);

	sread(in, val.string_information_, edition);
	sread(in, val.integer_information_, edition);
	sread(in, val.double_information_, edition);

	sread(in, val.string_vectors_, edition);
	sread(in, val.integer_vectors_, edition);
	sread(in, val.double_vectors_, edition);

	sread(in, SUBMODULES(val, Metadata), edition);
	sread(in, SUBMODULES(val, SparseInt), edition);
	sread(in, SUBMODULES(val, SparseDouble), edition);
	sread(in, SUBMODULES(val, Embedding), edition);
	sread(in, SUBMODULES(val, DifferentialAnalysis), edition);
	sread(in, SUBMODULES(val, GenomicRange), edition);
	sread(in, SUBMODULES(val, MotifPosition), edition);
	sread(in, SUBMODULES(val, CoverageTrack), edition);
	sread(in, SUBMODULES(val, Fragments), edition);
	sread(in, SUBMODULES(val, Cicero), edition);
	sread(in, SUBMODULES(val, Monocle3), edition);
};

void swrite(std::ofstream& out, const SingleCellAtac& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.species_, edition);
	swrite(out, val.random_state_, edition);

	swrite(out, val.string_information_, edition);
	swrite(out, val.integer_information_, edition);
	swrite(out, val.double_information_, edition);

	swrite(out, val.string_vectors_, edition);
	swrite(out, val.integer_vectors_, edition);
	swrite(out, val.double_vectors_, edition);

	swrite(out, SUBMODULES(val, Metadata), edition);
	swrite(out, SUBMODULES(val, SparseInt), edition);
	swrite(out, SUBMODULES(val, SparseDouble), edition);
	swrite(out, SUBMODULES(val, Embedding), edition);
	swrite(out, SUBMODULES(val, DifferentialAnalysis), edition);
	swrite(out, SUBMODULES(val, GenomicRange), edition);
	swrite(out, SUBMODULES(val, MotifPosition), edition);
	swrite(out, SUBMODULES(val, CoverageTrack), edition);
	swrite(out, SUBMODULES(val, Fragments), edition);
	swrite(out, SUBMODULES(val, Cicero), edition);
	swrite(out, SUBMODULES(val, Monocle3), edition);
};

void sread(std::ifstream& in, SingleCellRna& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.species_, edition);
	sread(in, val.random_state_, edition);

	sread(in, val.string_information_, edition);
	sread(in, val.integer_information_, edition);
	sread(in, val.double_information_, edition);

	sread(in, val.string_vectors_, edition);
	sread(in, val.integer_vectors_, edition);
	sread(in, val.double_vectors_, edition);

	sread(in, SUBMODULES(val, Metadata), edition);
	sread(in, SUBMODULES(val, SparseInt), edition);
	sread(in, SUBMODULES(val, SparseDouble), edition);
	sread(in, SUBMODULES(val, Embedding), edition);
	sread(in, SUBMODULES(val, DifferentialAnalysis), edition);
	sread(in, SUBMODULES(val, GSEA), edition);
	sread(in, SUBMODULES(val, CellChat), edition);
	sread(in, SUBMODULES(val, CNV), edition);
	sread(in, SUBMODULES(val, VelocytoBase), edition);
	sread(in, SUBMODULES(val, Monocle3), edition);
};

void swrite(std::ofstream& out, const SingleCellRna& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.species_, edition);
	swrite(out, val.random_state_, edition);

	swrite(out, val.string_information_, edition);
	swrite(out, val.integer_information_, edition);
	swrite(out, val.double_information_, edition);

	swrite(out, val.string_vectors_, edition);
	swrite(out, val.integer_vectors_, edition);
	swrite(out, val.double_vectors_, edition);

	swrite(out, SUBMODULES(val, Metadata), edition);
	swrite(out, SUBMODULES(val, SparseInt), edition);
	swrite(out, SUBMODULES(val, SparseDouble), edition);
	swrite(out, SUBMODULES(val, Embedding), edition);
	swrite(out, SUBMODULES(val, DifferentialAnalysis), edition);
	swrite(out, SUBMODULES(val, GSEA), edition);
	swrite(out, SUBMODULES(val, CellChat), edition);
	swrite(out, SUBMODULES(val, CNV), edition);
	swrite(out, SUBMODULES(val, VelocytoBase), edition);
	swrite(out, SUBMODULES(val, Monocle3), edition);
};

void sread(std::ifstream& in, BulkRna& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.species_, edition);
	sread(in, val.random_state_, edition);

	sread(in, val.string_information_, edition);
	sread(in, val.integer_information_, edition);
	sread(in, val.double_information_, edition);

	sread(in, val.string_vectors_, edition);
	sread(in, val.integer_vectors_, edition);
	sread(in, val.double_vectors_, edition);

	sread(in, SUBMODULES(val, Metadata), edition);
	sread(in, SUBMODULES(val, DenseInt), edition);
	sread(in, SUBMODULES(val, DenseDouble), edition);
	sread(in, SUBMODULES(val, Embedding), edition);
	sread(in, SUBMODULES(val, DifferentialAnalysis), edition);
	sread(in, SUBMODULES(val, GSEA), edition);
};

void swrite(std::ofstream& out, const BulkRna& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.species_, edition);
	swrite(out, val.random_state_, edition);

	swrite(out, val.string_information_, edition);
	swrite(out, val.integer_information_, edition);
	swrite(out, val.double_information_, edition);

	swrite(out, val.string_vectors_, edition);
	swrite(out, val.integer_vectors_, edition);
	swrite(out, val.double_vectors_, edition);

	swrite(out, SUBMODULES(val, Metadata), edition);
	swrite(out, SUBMODULES(val, DenseInt), edition);
	swrite(out, SUBMODULES(val, DenseDouble), edition);
	swrite(out, SUBMODULES(val, Embedding), edition);
	swrite(out, SUBMODULES(val, DifferentialAnalysis), edition);
	swrite(out, SUBMODULES(val, GSEA), edition);
};

void sread(std::ifstream& in, Cicero& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.connections_, edition);
	sread(in, val.regulation_groups_, edition);
	sread(in, val.peak_names_, edition);
	sread(in, val.regulation_group_counts_, edition);
	sread(in, val.regulation_group_normalized_, edition);

	sread(in, SUBMODULES(val, Embedding), edition);
	sread(in, SUBMODULES(val, DifferentialAnalysis), edition);
};

void swrite(std::ofstream& out, const Cicero& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.connections_, edition);
	swrite(out, val.regulation_groups_, edition);
	swrite(out, val.peak_names_, edition);
	swrite(out, val.regulation_group_counts_, edition);
	swrite(out, val.regulation_group_normalized_, edition);

	swrite(out, SUBMODULES(val, Embedding), edition);
	swrite(out, SUBMODULES(val, DifferentialAnalysis), edition);
};

void sread(std::ifstream& in, Monocle3& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.pr_graph_, edition);
	sread(in, val.cell_graph_, edition);

	sread(in, val.original_embedding_, edition);
	sread(in, val.cell_included_, edition);

	sread(in, val.cell_embedding_, edition);
	sread(in, val.pr_embedding_, edition);
	sread(in, val.cell_graph_weights_, edition);
	sread(in, val.pseudo_time_, edition);
};

void swrite(std::ofstream& out, const Monocle3& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.pr_graph_, edition);
	swrite(out, val.cell_graph_, edition);

	swrite(out, val.original_embedding_, edition);
	swrite(out, val.cell_included_, edition);

	swrite(out, val.cell_embedding_, edition);
	swrite(out, val.pr_embedding_, edition);
	swrite(out, val.cell_graph_weights_, edition);
	swrite(out, val.pseudo_time_, edition);
};

void sread(std::ifstream& in, Pando& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.mat_, edition);
	sread(in, SUBMODULES(val, GeneName), edition);
};

void swrite(std::ofstream& out, const Pando& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.mat_, edition);
	swrite(out, SUBMODULES(val, GeneName), edition);
};

void sread(std::ifstream& in, VelocytoBase& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, SUBMODULES(val, SparseInt), edition);
	sread(in, SUBMODULES(val, VelocityEstimate), edition);
	sread(in, SUBMODULES(val, ScveloEstimate), edition);
};

void swrite(std::ofstream& out, const VelocytoBase& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, SUBMODULES(val, SparseInt), edition);
	swrite(out, SUBMODULES(val, VelocityEstimate), edition);
	swrite(out, SUBMODULES(val, ScveloEstimate), edition);
};

void sread(std::ifstream& in, StringVector& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.data_, edition);
};

void swrite(std::ofstream& out, const StringVector& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.data_, edition);
};

void sread(std::ifstream& in, ScveloEstimate& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.graph_, edition);
	sread(in, val.graph_neg_, edition);
	sread(in, val.self_probability_, edition);
};

void swrite(std::ofstream& out, const ScveloEstimate& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.graph_, edition);
	swrite(out, val.graph_neg_, edition);
	swrite(out, val.self_probability_, edition);
};

void sread(std::ifstream& in, VelocityEstimate& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.projected_, edition);
	sread(in, val.current_, edition);
	sread(in, val.deltaE_, edition);
	sread(in, val.gene_names_, edition);
};

void swrite(std::ofstream& out, const VelocityEstimate& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.projected_, edition);
	swrite(out, val.current_, edition);
	swrite(out, val.deltaE_, edition);
	swrite(out, val.gene_names_, edition);
};

void sread(std::ifstream& in, CoverageTrack& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.insertion_matrix_, edition);
	sread(in, val.level_name_, edition);
	sread(in, val.levels_, edition);
	sread(in, val.annotation_, edition);
};

void swrite(std::ofstream& out, const CoverageTrack& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.insertion_matrix_, edition);
	swrite(out, val.level_name_, edition);
	swrite(out, val.levels_, edition);
	swrite(out, val.annotation_, edition);
};

void sread(std::ifstream& in, GenomicRange& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.sequence_names_, edition);
	sread(in, val.ranges_, edition);
	sread(in, val.strand_, edition);
	sread(in, val.sequence_information_, edition);
	sread(in, val.metadata_, edition);
};

void swrite(std::ofstream& out, const GenomicRange& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.sequence_names_, edition);
	swrite(out, val.ranges_, edition);
	swrite(out, val.strand_, edition);
	swrite(out, val.sequence_information_, edition);
	swrite(out, val.metadata_, edition);
};

void sread(std::ifstream& in, Fragments& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.data_, edition);
	sread(in, val.cell_names_, edition);
};

void swrite(std::ofstream& out, const Fragments& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.data_, edition);
	swrite(out, val.cell_names_, edition);
};

void sread(std::ifstream& in, SeqInfo& val, const QString& edition) {
	sread(in, val.sequence_names_, edition);
	sread(in, val.sequence_lengths_, edition);
	sread(in, val.is_circular_, edition);
	sread(in, val.genome_, edition);
};

void swrite(std::ofstream& out, const SeqInfo& val, const QString& edition) {
	swrite(out, val.sequence_names_, edition);
	swrite(out, val.sequence_lengths_, edition);
	swrite(out, val.is_circular_, edition);
	swrite(out, val.genome_, edition);
};

void sread(std::ifstream& in, IRange& val, const QString& edition) {
	sread(in, val.start_, edition);
	sread(in, val.width_, edition);
};

void swrite(std::ofstream& out, const IRange& val, const QString& edition) {
	swrite(out, val.start_, edition);
	swrite(out, val.width_, edition);
};

void sread(std::ifstream& in, CellChat& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.identity_, edition);
	sread(in, val.ligand_receptor_index_, edition);
	sread(in, val.pathway_index_, edition);
	sread(in, val.levels_, edition);
	sread(in, val.counts_, edition);
	sread(in, val.weights_, edition);
	sread(in, val.ligand_receptor_probability_, edition);
	sread(in, val.ligand_receptor_p_value_, edition);
	sread(in, val.pathway_probability_, edition);
	sread(in, val.interaction_summary_, edition);
	sread(in, val.pathway_summary_, edition);
};

void swrite(std::ofstream& out, const CellChat& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.identity_, edition);
	swrite(out, val.ligand_receptor_index_, edition);
	swrite(out, val.pathway_index_, edition);
	swrite(out, val.levels_, edition);
	swrite(out, val.counts_, edition);
	swrite(out, val.weights_, edition);
	swrite(out, val.ligand_receptor_probability_, edition);
	swrite(out, val.ligand_receptor_p_value_, edition);
	swrite(out, val.pathway_probability_, edition);
	swrite(out, val.interaction_summary_, edition);
	swrite(out, val.pathway_summary_, edition);
};

void sread(std::ifstream& in, CNV& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.mat_, edition);
	sread(in, val.chromosome_info_, edition);
	sread(in, val.cluster_info_, edition);
};

void swrite(std::ofstream& out, const CNV& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.mat_, edition);
	swrite(out, val.chromosome_info_, edition);
	swrite(out, val.cluster_info_, edition);
};

void sread(std::ifstream& in, MotifPosition& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.motifs_, edition);
	sread(in, val.peak_locations_, edition);
	sread(in, val.motif_information_, edition);
	sread(in, val.peak_names_, edition);
	sread(in, val.motif_names_, edition);
	sread(in, SUBMODULES(val, Footprint), edition);
	sread(in, SUBMODULES(val, ChromVAR), edition);
};

void swrite(std::ofstream& out, const MotifPosition& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.motifs_, edition);
	swrite(out, val.peak_locations_, edition);
	swrite(out, val.motif_information_, edition);
	swrite(out, val.peak_names_, edition);
	swrite(out, val.motif_names_, edition);
	swrite(out, SUBMODULES(val, Footprint), edition);
	swrite(out, SUBMODULES(val, ChromVAR), edition);
};

void sread(std::ifstream& in, MotifPosition::match& val, const QString& edition) {
	sread(in, val.peak_index, edition);
	sread(in, val.position, edition);
	sread(in, val.score, edition);
	sread(in, val.strand, edition);
};

void swrite(std::ofstream& out, const MotifPosition::match& val, const QString& edition) {
	swrite(out, val.peak_index, edition);
	swrite(out, val.position, edition);
	swrite(out, val.score, edition);
	swrite(out, val.strand, edition);
};

void sread(std::ifstream& in, ChromVAR& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.motif_names_, edition);
	sread(in, val.z_, edition);
	sread(in, val.dev_, edition);
	sread(in, SUBMODULES(val, DifferentialAnalysis), edition);
	sread(in, SUBMODULES(val, Embedding), edition);
};

void swrite(std::ofstream& out, const ChromVAR& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.motif_names_, edition);
	swrite(out, val.z_, edition);
	swrite(out, val.dev_, edition);
	swrite(out, SUBMODULES(val, DifferentialAnalysis), edition);
	swrite(out, SUBMODULES(val, Embedding), edition);
};

void sread(std::ifstream& in, Footprint& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.motif_, edition);
	sread(in, val.insertion_matrix_, edition);
	sread(in, val.expected_insertions_, edition);
	sread(in, val.motif_location_, edition);
	sread(in, SUBMODULES(val, GeneName), edition);
};

void swrite(std::ofstream& out, const Footprint& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.motif_, edition);
	swrite(out, val.insertion_matrix_, edition);
	swrite(out, val.expected_insertions_, edition);
	swrite(out, val.motif_location_, edition);
	swrite(out, SUBMODULES(val, GeneName), edition);
};

void sread(std::ifstream& in, GeneName& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.data_, edition);
	sread(in, SUBMODULES(val, Enrichment), edition);
};

void swrite(std::ofstream& out, const GeneName& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.data_, edition);
	swrite(out, SUBMODULES(val, Enrichment), edition);
};

void sread(std::ifstream& in, DifferentialAnalysis& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.mat_, edition);
	sread(in, SUBMODULES(val, Enrichment), edition);
};

void swrite(std::ofstream& out, const DifferentialAnalysis& val, const QString& edition) {

	swrite(out, val.data_type_, edition);
	swrite(out, val.mat_, edition);
	swrite(out, SUBMODULES(val, Enrichment), edition);
};

void sread(std::ifstream& in, GSEA& val, const QString& edition) {
	sread(in, val.data_type_, edition);

	sread(in, val.mat_, edition);
	sread(in, val.database_, edition);
	sread(in, val.comparison_, edition);
	sread(in, val.pathways_, edition);
	sread(in, val.correlations_, edition);
	sread(in, val.bar_points_, edition);
	sread(in, val.bar_colors_, edition);
	sread(in, val.gene_location_, edition);
	sread(in, val.point_x_, edition);
	sread(in, val.point_y_, edition);
};

void swrite(std::ofstream& out, const GSEA& val, const QString& edition) {
	swrite(out, val.data_type_, edition);

	swrite(out, val.mat_, edition);
	swrite(out, val.database_, edition);
	swrite(out, val.comparison_, edition);
	swrite(out, val.pathways_, edition);
	swrite(out, val.correlations_, edition);
	swrite(out, val.bar_points_, edition);
	swrite(out, val.bar_colors_, edition);
	swrite(out, val.gene_location_, edition);
	swrite(out, val.point_x_, edition);
	swrite(out, val.point_y_, edition);
};

void sread(std::ifstream& in, SparseInt& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.rownames_, edition);
	sread(in, val.colnames_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const SparseInt& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.rownames_, edition);
	swrite(out, val.colnames_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, SparseDouble& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.rownames_, edition);
	sread(in, val.colnames_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const SparseDouble& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.rownames_, edition);
	swrite(out, val.colnames_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, DenseInt& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.rownames_, edition);
	sread(in, val.colnames_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const DenseInt& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.rownames_, edition);
	swrite(out, val.colnames_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, DenseDouble& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.rownames_, edition);
	sread(in, val.colnames_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const DenseDouble& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.rownames_, edition);
	swrite(out, val.colnames_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, Embedding& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.data_, edition);
};

void swrite(std::ofstream& out, const Embedding& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.data_, edition);
};

void sread(std::ifstream& in, NumericMatrix& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.data_, edition);
};

void swrite(std::ofstream& out, const NumericMatrix& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.data_, edition);
};

void sread(std::ifstream& in, Metadata& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const Metadata& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, Enrichment& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const Enrichment& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, DataFrame& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.mat_, edition);
};

void swrite(std::ofstream& out, const DataFrame& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.mat_, edition);
};

void sread(std::ifstream& in, PatternWeightMatrix& val, const QString& edition) {
	sread(in, val.data_type_, edition);
	sread(in, val.weight_, edition);
	sread(in, val.motif_name_, edition);
	sread(in, val.transcriptional_factor_name_, edition);
};

void swrite(std::ofstream& out, const PatternWeightMatrix& val, const QString& edition) {
	swrite(out, val.data_type_, edition);
	swrite(out, val.weight_, edition);
	swrite(out, val.motif_name_, edition);
	swrite(out, val.transcriptional_factor_name_, edition);
};

void sread(std::ifstream& in, QColor& val, const QString& edition) {

	int r, g, b, a;
	sread(in, r, edition);
	sread(in, g, edition);
	sread(in, b, edition);
	sread(in, a, edition);

	val = QColor(r, g, b, a);
};

void swrite(std::ofstream& out, const QColor& val, const QString& edition) {

	int r, g, b, a;
	val.getRgb(&r, &g, &b);

	a = val.alpha();
	swrite(out, r, edition);
	swrite(out, g, edition);
	swrite(out, b, edition);
	swrite(out, a, edition);
};

void sread(std::ifstream& in, QString& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size, edition);
	if (size > 0) {

		// for COW, can not read directly

		std::vector<ushort> buffer(size);
		in.read(reinterpret_cast<char*>(buffer.data()), size * sizeof(ushort));

		val = QString::fromUtf16(buffer.data(), size);
	}
};

void swrite(std::ofstream& out, const QString& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size, edition);
	if (size > 0) {
		out.write(reinterpret_cast<const char*>(val.utf16()), size * sizeof(ushort));
	}
};

void sread(std::ifstream& in, std::string& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size, edition);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(char));
};

void swrite(std::ofstream& out, const std::string& val, const QString& edition) {
	std::size_t size = val.size();
	swrite(out, size, edition);

	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(char));
};

void sread(std::ifstream& in, QFont& val, const QString& edition) {
	QString family, style_name;
	int point_size, pixel_size;

	sread(in, family, edition);
	sread(in, style_name, edition);
	sread(in, point_size, edition);
	sread(in, pixel_size, edition);

	val.setFamily(family);
	val.setStyleName(style_name);
	val.setPointSize(point_size);
	val.setPixelSize(pixel_size);
};

void swrite(std::ofstream& out, const QFont& val, const QString& edition) {
	swrite(out, val.family(), edition);
	swrite(out, val.styleName(), edition);
	swrite(out, val.pointSize(), edition);
	swrite(out, val.pixelSize(), edition);
};

void sread(std::ifstream& in, igraph_t& val, const QString& edition) {

	igraph_empty(&val, 0, false);

	sread(in, val.n, edition);
	sread(in, val.directed, edition);
	sread(in, val.from, edition);
	sread(in, val.to, edition);
	sread(in, val.oi, edition);
	sread(in, val.ii, edition);
	sread(in, val.os, edition);
	sread(in, val.is, edition);
};

void swrite(std::ofstream& out, const igraph_t& val, const QString& edition) {
	swrite(out, val.n, edition);
	swrite(out, val.directed, edition);
	swrite(out, val.from, edition);
	swrite(out, val.to, edition);
	swrite(out, val.oi, edition);
	swrite(out, val.ii, edition);
	swrite(out, val.os, edition);
	swrite(out, val.is, edition);
};

void sread(std::ifstream& in, CustomMatrix& val, const QString& edition) {
	sread(in, val.rownames_, edition);
	sread(in, val.colnames_, edition);
	sread(in, val.data_type_, edition);
	sread(in, val.string_factors_, edition);
	sread(in, val.integer_factors_, edition);

	int size = val.colnames_.size();

	for (int i = 0; i < size; ++i) {

		QString name;
		sread(in, name, edition);

		CustomMatrix::DataType data_type = val.data_type_[name];
		if (data_type == CustomMatrix::DataType::DoubleNumeric) {
			QVector<double>* ptr = new QVector<double>();
			sread(in, *ptr, edition);
			val.data_[name] = ptr;
		}
		else if (data_type == CustomMatrix::DataType::IntegerNumeric || data_type == CustomMatrix::DataType::IntegerFactor) {
			QVector<int>* ptr = new QVector<int>();
			sread(in, *ptr, edition);
			val.data_[name] = ptr;
		}
		else if (data_type == CustomMatrix::DataType::QStringFactor || data_type == CustomMatrix::DataType::QString) {
			QStringList* ptr = new QStringList();
			sread(in, *ptr, edition);
			val.data_[name] = ptr;
		}
	}
};

void swrite(std::ofstream& out, const CustomMatrix& val, const QString& edition) {

	swrite(out, val.rownames_, edition);
	swrite(out, val.colnames_, edition);
	swrite(out, val.data_type_, edition);
	swrite(out, val.string_factors_, edition);
	swrite(out, val.integer_factors_, edition);

	for (const QString& colname : val.colnames_) {
		swrite(out, colname, edition);
		CustomMatrix::DataType data_type = val.data_type_.at(colname);

		if (data_type == CustomMatrix::DataType::DoubleNumeric) {
			const QVector<double>& ref = *static_cast<QVector<double>*>(val.data_.at(colname));
			swrite(out, ref, edition);
		}
		else if (data_type == CustomMatrix::DataType::IntegerNumeric || data_type == CustomMatrix::DataType::IntegerFactor) {
			const QVector<int>& ref = *static_cast<QVector<int>*>(val.data_.at(colname));
			swrite(out, ref, edition);
		}
		else if (data_type == CustomMatrix::DataType::QStringFactor || data_type == CustomMatrix::DataType::QString) {
			const QStringList& ref = *static_cast<QStringList*>(val.data_.at(colname));
			swrite(out, ref, edition);
		}
	}
};

void sread(std::ifstream& in, igraph_vector_int_t& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size, edition);
	igraph_vector_int_init(&val, size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, VECTOR(val)[i], edition);
	}
};

void swrite(std::ofstream& out, const igraph_vector_int_t& val, const QString& edition) {
	std::size_t size = igraph_vector_int_size(&val);
	swrite(out, size, edition);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, VECTOR(val)[i], edition);
	}
};

void sread(std::ifstream& in, igraph_vector_t& val, const QString& edition) {
	std::size_t size{ 0 };
	sread(in, size, edition);
	igraph_vector_init(&val, size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, VECTOR(val)[i], edition);
	}
};

void swrite(std::ofstream& out, const igraph_vector_t& val, const QString& edition) {
	std::size_t size = igraph_vector_size(&val);
	swrite(out, size, edition);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, VECTOR(val)[i], edition);
	}
};