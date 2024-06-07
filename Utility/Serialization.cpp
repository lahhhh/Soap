#include "Serialization.h"

void sread(std::ifstream& in, GraphSettings& val) {
	sread(in, val.string_info_);
	sread(in, val.string_list_info_);
	sread(in, val.font_info_);
	sread(in, val.integer_info_);
	sread(in, val.bool_info_);
	sread(in, val.color_info_);
	sread(in, val.palette_);
};

void swrite(std::ofstream& out, const GraphSettings& val) {
	swrite(out, val.string_info_);
	swrite(out, val.string_list_info_);
	swrite(out, val.font_info_);
	swrite(out, val.integer_info_);
	swrite(out, val.bool_info_);
	swrite(out, val.color_info_);
	swrite(out, val.palette_);
};

void sread(std::ifstream& in, SingleCellMultiome& val) {
	sread(in, val.data_type_);
	sread(in, val.species_);
	sread(in, val.random_state_);

	sread(in, val.string_information_);
	sread(in, val.integer_information_);
	sread(in, val.double_information_);

	sread(in, val.string_vectors_);
	sread(in, val.integer_vectors_);
	sread(in, val.double_vectors_);

	sread(in, SUBMODULES(val, DataField));
	sread(in, SUBMODULES(val, Metadata));
	sread(in, SUBMODULES(val, GSEA));
	sread(in, SUBMODULES(val, CellChat));
	sread(in, SUBMODULES(val, CNV));
	sread(in, SUBMODULES(val, GenomicRange));
	sread(in, SUBMODULES(val, MotifPosition));
	sread(in, SUBMODULES(val, CoverageTrack));
	sread(in, SUBMODULES(val, Fragments));
	sread(in, SUBMODULES(val, VelocytoBase));
	sread(in, SUBMODULES(val, Pando));
	sread(in, SUBMODULES(val, Monocle3));
	sread(in, SUBMODULES(val, Cicero));
};

void swrite(std::ofstream& out, const SingleCellMultiome& val) {
	swrite(out, val.data_type_);
	swrite(out, val.species_);
	swrite(out, val.random_state_);

	swrite(out, val.string_information_);
	swrite(out, val.integer_information_);
	swrite(out, val.double_information_);

	swrite(out, val.string_vectors_);
	swrite(out, val.integer_vectors_);
	swrite(out, val.double_vectors_);

	swrite(out, SUBMODULES(val, DataField));
	swrite(out, SUBMODULES(val, Metadata));
	swrite(out, SUBMODULES(val, GSEA));
	swrite(out, SUBMODULES(val, CellChat));
	swrite(out, SUBMODULES(val, CNV));
	swrite(out, SUBMODULES(val, GenomicRange));
	swrite(out, SUBMODULES(val, MotifPosition));
	swrite(out, SUBMODULES(val, CoverageTrack));
	swrite(out, SUBMODULES(val, Fragments));
	swrite(out, SUBMODULES(val, VelocytoBase));
	swrite(out, SUBMODULES(val, Pando));
	swrite(out, SUBMODULES(val, Monocle3));
	swrite(out, SUBMODULES(val, Cicero));
};

void sread(std::ifstream& in, DataField& val) {
	sread(in, val.data_type_);
	sread(in, SUBMODULES(val, SparseInt));
	sread(in, SUBMODULES(val, SparseDouble));
	sread(in, SUBMODULES(val, Embedding));
	sread(in, SUBMODULES(val, DifferentialAnalysis));
};

void swrite(std::ofstream& out, const DataField& val) {
	swrite(out, val.data_type_);
	swrite(out, SUBMODULES(val, SparseInt));
	swrite(out, SUBMODULES(val, SparseDouble));
	swrite(out, SUBMODULES(val, Embedding));
	swrite(out, SUBMODULES(val, DifferentialAnalysis));
};

void sread(std::ifstream& in, SingleCellAtac& val) {
	sread(in, val.data_type_);
	sread(in, val.species_);
	sread(in, val.random_state_);

	sread(in, val.string_information_);
	sread(in, val.integer_information_);
	sread(in, val.double_information_);

	sread(in, val.string_vectors_);
	sread(in, val.integer_vectors_);
	sread(in, val.double_vectors_);

	sread(in, SUBMODULES(val, Metadata));
	sread(in, SUBMODULES(val, SparseInt));
	sread(in, SUBMODULES(val, SparseDouble));
	sread(in, SUBMODULES(val, Embedding));
	sread(in, SUBMODULES(val, DifferentialAnalysis));
	sread(in, SUBMODULES(val, GenomicRange));
	sread(in, SUBMODULES(val, MotifPosition));
	sread(in, SUBMODULES(val, CoverageTrack));
	sread(in, SUBMODULES(val, Fragments));
	sread(in, SUBMODULES(val, Cicero));
	sread(in, SUBMODULES(val, Monocle3));
};

void swrite(std::ofstream& out, const SingleCellAtac& val) {
	swrite(out, val.data_type_);
	swrite(out, val.species_);
	swrite(out, val.random_state_);

	swrite(out, val.string_information_);
	swrite(out, val.integer_information_);
	swrite(out, val.double_information_);

	swrite(out, val.string_vectors_);
	swrite(out, val.integer_vectors_);
	swrite(out, val.double_vectors_);

	swrite(out, SUBMODULES(val, Metadata));
	swrite(out, SUBMODULES(val, SparseInt));
	swrite(out, SUBMODULES(val, SparseDouble));
	swrite(out, SUBMODULES(val, Embedding));
	swrite(out, SUBMODULES(val, DifferentialAnalysis));
	swrite(out, SUBMODULES(val, GenomicRange));
	swrite(out, SUBMODULES(val, MotifPosition));
	swrite(out, SUBMODULES(val, CoverageTrack));
	swrite(out, SUBMODULES(val, Fragments));
	swrite(out, SUBMODULES(val, Cicero));
	swrite(out, SUBMODULES(val, Monocle3));
};

void sread(std::ifstream& in, SingleCellRna& val) {
	sread(in, val.data_type_);
	sread(in, val.species_);
	sread(in, val.random_state_);

	sread(in, val.string_information_);
	sread(in, val.integer_information_);
	sread(in, val.double_information_);

	sread(in, val.string_vectors_);
	sread(in, val.integer_vectors_);
	sread(in, val.double_vectors_);

	sread(in, SUBMODULES(val, Metadata));
	sread(in, SUBMODULES(val, SparseInt));
	sread(in, SUBMODULES(val, SparseDouble));
	sread(in, SUBMODULES(val, Embedding));
	sread(in, SUBMODULES(val, DifferentialAnalysis));
	sread(in, SUBMODULES(val, GSEA));
	sread(in, SUBMODULES(val, CellChat));
	sread(in, SUBMODULES(val, CNV));
	sread(in, SUBMODULES(val, VelocytoBase));
	sread(in, SUBMODULES(val, Monocle3));
};

void swrite(std::ofstream& out, const SingleCellRna& val) {
	swrite(out, val.data_type_);
	swrite(out, val.species_);
	swrite(out, val.random_state_);

	swrite(out, val.string_information_);
	swrite(out, val.integer_information_);
	swrite(out, val.double_information_);

	swrite(out, val.string_vectors_);
	swrite(out, val.integer_vectors_);
	swrite(out, val.double_vectors_);

	swrite(out, SUBMODULES(val, Metadata));
	swrite(out, SUBMODULES(val, SparseInt));
	swrite(out, SUBMODULES(val, SparseDouble));
	swrite(out, SUBMODULES(val, Embedding));
	swrite(out, SUBMODULES(val, DifferentialAnalysis));
	swrite(out, SUBMODULES(val, GSEA));
	swrite(out, SUBMODULES(val, CellChat));
	swrite(out, SUBMODULES(val, CNV));
	swrite(out, SUBMODULES(val, VelocytoBase));
	swrite(out, SUBMODULES(val, Monocle3));
};

void sread(std::ifstream& in, BulkRna& val) {
	sread(in, val.data_type_);
	sread(in, val.species_);
	sread(in, val.random_state_);

	sread(in, val.string_information_);
	sread(in, val.integer_information_);
	sread(in, val.double_information_);

	sread(in, val.string_vectors_);
	sread(in, val.integer_vectors_);
	sread(in, val.double_vectors_);

	sread(in, SUBMODULES(val, Metadata));
	sread(in, SUBMODULES(val, DenseInt));
	sread(in, SUBMODULES(val, DenseDouble));
	sread(in, SUBMODULES(val, Embedding));
	sread(in, SUBMODULES(val, DifferentialAnalysis));
	sread(in, SUBMODULES(val, GSEA));
};

void swrite(std::ofstream& out, const BulkRna& val) {
	swrite(out, val.data_type_);
	swrite(out, val.species_);
	swrite(out, val.random_state_);

	swrite(out, val.string_information_);
	swrite(out, val.integer_information_);
	swrite(out, val.double_information_);

	swrite(out, val.string_vectors_);
	swrite(out, val.integer_vectors_);
	swrite(out, val.double_vectors_);

	swrite(out, SUBMODULES(val, Metadata));
	swrite(out, SUBMODULES(val, DenseInt));
	swrite(out, SUBMODULES(val, DenseDouble));
	swrite(out, SUBMODULES(val, Embedding));
	swrite(out, SUBMODULES(val, DifferentialAnalysis));
	swrite(out, SUBMODULES(val, GSEA));
};

void sread(std::ifstream& in, Cicero& val) {
	sread(in, val.data_type_);

	sread(in, val.connections_);
	sread(in, val.regulation_groups_);
	sread(in, val.peak_names_);
	sread(in, val.regulation_group_counts_);
	sread(in, val.regulation_group_normalized_);

	sread(in, SUBMODULES(val, Embedding));
	sread(in, SUBMODULES(val, DifferentialAnalysis));
};

void swrite(std::ofstream& out, const Cicero& val) {
	swrite(out, val.data_type_);

	swrite(out, val.connections_);
	swrite(out, val.regulation_groups_);
	swrite(out, val.peak_names_);
	swrite(out, val.regulation_group_counts_);
	swrite(out, val.regulation_group_normalized_);

	swrite(out, SUBMODULES(val, Embedding));
	swrite(out, SUBMODULES(val, DifferentialAnalysis));
};

void sread(std::ifstream& in, Monocle3& val) {
	sread(in, val.data_type_);

	sread(in, val.pr_graph_);
	sread(in, val.cell_graph_);

	sread(in, val.original_embedding_);
	sread(in, val.cell_included_);

	sread(in, val.cell_embedding_);
	sread(in, val.pr_embedding_);
	sread(in, val.cell_graph_weights_);
	sread(in, val.pseudo_time_);
};

void swrite(std::ofstream& out, const Monocle3& val) {
	swrite(out, val.data_type_);

	swrite(out, val.pr_graph_);
	swrite(out, val.cell_graph_);

	swrite(out, val.original_embedding_);
	swrite(out, val.cell_included_);

	swrite(out, val.cell_embedding_);
	swrite(out, val.pr_embedding_);
	swrite(out, val.cell_graph_weights_);
	swrite(out, val.pseudo_time_);
};

void sread(std::ifstream& in, Pando& val) {
	sread(in, val.data_type_);

	sread(in, val.mat_);
	sread(in, SUBMODULES(val, GeneName));
};

void swrite(std::ofstream& out, const Pando& val) {
	swrite(out, val.data_type_);

	swrite(out, val.mat_);
	swrite(out, SUBMODULES(val, GeneName));
};

void sread(std::ifstream& in, VelocytoBase& val) {
	sread(in, val.data_type_);

	sread(in, SUBMODULES(val, SparseInt));
	sread(in, SUBMODULES(val, VelocityEstimate));
	sread(in, SUBMODULES(val, ScveloEstimate));
};

void swrite(std::ofstream& out, const VelocytoBase& val) {
	swrite(out, val.data_type_);

	swrite(out, SUBMODULES(val, SparseInt));
	swrite(out, SUBMODULES(val, VelocityEstimate));
	swrite(out, SUBMODULES(val, ScveloEstimate));
};

void sread(std::ifstream& in, StringVector& val) {
	sread(in, val.data_type_);
	sread(in, val.data_);
};

void swrite(std::ofstream& out, const StringVector& val) {
	swrite(out, val.data_type_);
	swrite(out, val.data_);
};

void sread(std::ifstream& in, ScveloEstimate& val) {
	sread(in, val.data_type_);

	sread(in, val.graph_);
	sread(in, val.graph_neg_);
	sread(in, val.self_probability_);
};

void swrite(std::ofstream& out, const ScveloEstimate& val) {
	swrite(out, val.data_type_);

	swrite(out, val.graph_);
	swrite(out, val.graph_neg_);
	swrite(out, val.self_probability_);
};

void sread(std::ifstream& in, VelocityEstimate& val) {
	sread(in, val.data_type_);

	sread(in, val.projected_);
	sread(in, val.current_);
	sread(in, val.deltaE_);
	sread(in, val.gene_names_);
};

void swrite(std::ofstream& out, const VelocityEstimate& val) {
	swrite(out, val.data_type_);

	swrite(out, val.projected_);
	swrite(out, val.current_);
	swrite(out, val.deltaE_);
	swrite(out, val.gene_names_);
};

void sread(std::ifstream& in, CoverageTrack& val) {
	sread(in, val.data_type_);

	sread(in, val.insertion_matrix_);
	sread(in, val.level_name_);
	sread(in, val.levels_);
	sread(in, val.annotation_);
};

void swrite(std::ofstream& out, const CoverageTrack& val) {
	swrite(out, val.data_type_);

	swrite(out, val.insertion_matrix_);
	swrite(out, val.level_name_);
	swrite(out, val.levels_);
	swrite(out, val.annotation_);
};

void sread(std::ifstream& in, GenomicRange& val) {
	sread(in, val.data_type_);

	sread(in, val.sequence_names_);
	sread(in, val.ranges_);
	sread(in, val.strand_);
	sread(in, val.sequence_information_);
	sread(in, val.metadata_);
};

void swrite(std::ofstream& out, const GenomicRange& val) {
	swrite(out, val.data_type_);

	swrite(out, val.sequence_names_);
	swrite(out, val.ranges_);
	swrite(out, val.strand_);
	swrite(out, val.sequence_information_);
	swrite(out, val.metadata_);
};

void sread(std::ifstream& in, Fragments& val) {
	sread(in, val.data_type_);

	sread(in, val.data_);
	sread(in, val.cell_names_);
};

void swrite(std::ofstream& out, const Fragments& val) {
	swrite(out, val.data_type_);

	swrite(out, val.data_);
	swrite(out, val.cell_names_);
};

void sread(std::ifstream& in, SeqInfo& val) {
	sread(in, val.sequence_names_);
	sread(in, val.sequence_lengths_);
	sread(in, val.is_circular_);
	sread(in, val.genome_);
};

void swrite(std::ofstream& out, const SeqInfo& val) {
	swrite(out, val.sequence_names_);
	swrite(out, val.sequence_lengths_);
	swrite(out, val.is_circular_);
	swrite(out, val.genome_);
};

void sread(std::ifstream& in, IRange& val) {
	sread(in, val.start_);
	sread(in, val.width_);
};

void swrite(std::ofstream& out, const IRange& val) {
	swrite(out, val.start_);
	swrite(out, val.width_);
};

void sread(std::ifstream& in, CellChat& val) {
	sread(in, val.data_type_);

	sread(in, val.identity_);
	sread(in, val.ligand_receptor_index_);
	sread(in, val.pathway_index_);
	sread(in, val.levels_);
	sread(in, val.counts_);
	sread(in, val.weights_);
	sread(in, val.ligand_receptor_probability_);
	sread(in, val.ligand_receptor_p_value_);
	sread(in, val.pathway_probability_);
};

void swrite(std::ofstream& out, const CellChat& val) {
	swrite(out, val.data_type_);

	swrite(out, val.identity_);
	swrite(out, val.ligand_receptor_index_);
	swrite(out, val.pathway_index_);
	swrite(out, val.levels_);
	swrite(out, val.counts_);
	swrite(out, val.weights_);
	swrite(out, val.ligand_receptor_probability_);
	swrite(out, val.ligand_receptor_p_value_);
	swrite(out, val.pathway_probability_);
};

void sread(std::ifstream& in, CNV& val) {
	sread(in, val.data_type_);

	sread(in, val.mat_);
	sread(in, val.chromosome_info_);
	sread(in, val.cluster_info_);
};

void swrite(std::ofstream& out, const CNV& val) {
	swrite(out, val.data_type_);

	swrite(out, val.mat_);
	swrite(out, val.chromosome_info_);
	swrite(out, val.cluster_info_);
};

void sread(std::ifstream& in, MotifPosition& val) {
	sread(in, val.data_type_);

	sread(in, val.motifs_);
	sread(in, val.peak_locations_);
	sread(in, val.motif_information_);
	sread(in, val.peak_names_);
	sread(in, val.motif_names_);
	sread(in, SUBMODULES(val, Footprint));
	sread(in, SUBMODULES(val, ChromVAR));
};

void swrite(std::ofstream& out, const MotifPosition& val) {
	swrite(out, val.data_type_);

	swrite(out, val.motifs_);
	swrite(out, val.peak_locations_);
	swrite(out, val.motif_information_);
	swrite(out, val.peak_names_);
	swrite(out, val.motif_names_);
	swrite(out, SUBMODULES(val, Footprint));
	swrite(out, SUBMODULES(val, ChromVAR));
};

void sread(std::ifstream& in, MotifPosition::match& val) {
	sread(in, val.peak_index);
	sread(in, val.position);
	sread(in, val.score);
	sread(in, val.strand);
};

void swrite(std::ofstream& out, const MotifPosition::match& val) {
	swrite(out, val.peak_index);
	swrite(out, val.position);
	swrite(out, val.score);
	swrite(out, val.strand);
};

void sread(std::ifstream& in, ChromVAR& val) {
	sread(in, val.data_type_);

	sread(in, val.motif_names_);
	sread(in, val.z_);
	sread(in, val.dev_);
	sread(in, SUBMODULES(val, DifferentialAnalysis));
	sread(in, SUBMODULES(val, Embedding));
};

void swrite(std::ofstream& out, const ChromVAR& val) {
	swrite(out, val.data_type_);

	swrite(out, val.motif_names_);
	swrite(out, val.z_);
	swrite(out, val.dev_);
	swrite(out, SUBMODULES(val, DifferentialAnalysis));
	swrite(out, SUBMODULES(val, Embedding));
};

void sread(std::ifstream& in, Footprint& val) {
	sread(in, val.data_type_);

	sread(in, val.motif_);
	sread(in, val.insertion_matrix_);
	sread(in, val.expected_insertions_);
	sread(in, val.motif_location_);
	sread(in, SUBMODULES(val, GeneName));
};

void swrite(std::ofstream& out, const Footprint& val) {
	swrite(out, val.data_type_);

	swrite(out, val.motif_);
	swrite(out, val.insertion_matrix_);
	swrite(out, val.expected_insertions_);
	swrite(out, val.motif_location_);
	swrite(out, SUBMODULES(val, GeneName));
};

void sread(std::ifstream& in, GeneName& val) {
	sread(in, val.data_type_);

	sread(in, val.data_);
	sread(in, SUBMODULES(val, Enrichment));
};

void swrite(std::ofstream& out, const GeneName& val) {
	swrite(out, val.data_type_);

	swrite(out, val.data_);
	swrite(out, SUBMODULES(val, Enrichment));
};

void sread(std::ifstream& in, DifferentialAnalysis& val) {
	sread(in, val.data_type_);

	sread(in, val.mat_);
	sread(in, SUBMODULES(val, Enrichment));
};

void swrite(std::ofstream& out, const DifferentialAnalysis& val) {
	swrite(out, val.data_type_);

	swrite(out, val.mat_);
	swrite(out, SUBMODULES(val, Enrichment));
};

void sread(std::ifstream& in, GSEA& val) {
	sread(in, val.data_type_);

	sread(in, val.mat_);
	sread(in, val.database_);
	sread(in, val.comparison_);
	sread(in, val.pathways_);
	sread(in, val.correlations_);
	sread(in, val.bar_points_);
	sread(in, val.bar_colors_);
	sread(in, val.gene_location_);
	sread(in, val.point_x_);
	sread(in, val.point_y_);
};

void swrite(std::ofstream& out, const GSEA& val) {
	swrite(out, val.data_type_);

	swrite(out, val.mat_);
	swrite(out, val.database_);
	swrite(out, val.comparison_);
	swrite(out, val.pathways_);
	swrite(out, val.correlations_);
	swrite(out, val.bar_points_);
	swrite(out, val.bar_colors_);
	swrite(out, val.gene_location_);
	swrite(out, val.point_x_);
	swrite(out, val.point_y_);
};

void sread(std::ifstream& in, SparseInt& val) {
	sread(in, val.data_type_);
	sread(in, val.rownames_);
	sread(in, val.colnames_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const SparseInt& val) {
	swrite(out, val.data_type_);
	swrite(out, val.rownames_);
	swrite(out, val.colnames_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, SparseDouble& val) {
	sread(in, val.data_type_);
	sread(in, val.rownames_);
	sread(in, val.colnames_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const SparseDouble& val) {
	swrite(out, val.data_type_);
	swrite(out, val.rownames_);
	swrite(out, val.colnames_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, DenseInt& val) {
	sread(in, val.data_type_);
	sread(in, val.rownames_);
	sread(in, val.colnames_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const DenseInt& val) {
	swrite(out, val.data_type_);
	swrite(out, val.rownames_);
	swrite(out, val.colnames_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, DenseDouble& val) {
	sread(in, val.data_type_);
	sread(in, val.rownames_);
	sread(in, val.colnames_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const DenseDouble& val) {
	swrite(out, val.data_type_);
	swrite(out, val.rownames_);
	swrite(out, val.colnames_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, Embedding& val) {
	sread(in, val.data_type_);
	sread(in, val.data_);
};

void swrite(std::ofstream& out, const Embedding& val) {
	swrite(out, val.data_type_);
	swrite(out, val.data_);
};

void sread(std::ifstream& in, NumericMatrix& val) {
	sread(in, val.data_type_);
	sread(in, val.data_);
};

void swrite(std::ofstream& out, const NumericMatrix& val) {
	swrite(out, val.data_type_);
	swrite(out, val.data_);
};

void sread(std::ifstream& in, Metadata& val) {
	sread(in, val.data_type_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const Metadata& val) {
	swrite(out, val.data_type_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, Enrichment& val) {
	sread(in, val.data_type_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const Enrichment& val) {
	swrite(out, val.data_type_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, DataFrame& val) {
	sread(in, val.data_type_);
	sread(in, val.mat_);
};

void swrite(std::ofstream& out, const DataFrame& val) {
	swrite(out, val.data_type_);
	swrite(out, val.mat_);
};

void sread(std::ifstream& in, PatternWeightMatrix& val) {
	sread(in, val.data_type_);
	sread(in, val.weight_);
	sread(in, val.motif_name_);
	sread(in, val.transcriptional_factor_name_);
};

void swrite(std::ofstream& out, const PatternWeightMatrix& val) {
	swrite(out, val.data_type_);
	swrite(out, val.weight_);
	swrite(out, val.motif_name_);
	swrite(out, val.transcriptional_factor_name_);
};

void sread(std::ifstream& in, QColor& val) {

	int r, g, b, a;
	sread(in, r);
	sread(in, g);
	sread(in, b);
	sread(in, a);

	val = QColor(r, g, b, a);
};

void swrite(std::ofstream& out, const QColor& val) {

	int r, g, b, a;
	val.getRgb(&r, &g, &b);

	a = val.alpha();
	swrite(out, r);
	swrite(out, g);
	swrite(out, b);
	swrite(out, a);
};

void sread(std::ifstream& in, QString& val) {
	std::size_t size{ 0 };
	sread(in, size);
	if (size > 0) {

		// for COW, can not read directly

		std::vector<ushort> buffer(size);
		in.read(reinterpret_cast<char*>(buffer.data()), size * sizeof(ushort));

		val = QString::fromUtf16(buffer.data(), size);
	}
};

void swrite(std::ofstream& out, const QString& val) {
	std::size_t size = val.size();
	swrite(out, size);
	if (size > 0) {
		out.write(reinterpret_cast<const char*>(val.utf16()), size * sizeof(ushort));
	}
};

void sread(std::ifstream& in, std::string& val) {
	std::size_t size{ 0 };
	sread(in, size);
	val.resize(size);
	in.read(reinterpret_cast<char*>(val.data()), size * sizeof(char));
};

void swrite(std::ofstream& out, const std::string& val) {
	std::size_t size = val.size();
	swrite(out, size);

	out.write(reinterpret_cast<const char*>(val.data()), size * sizeof(char));
};

void sread(std::ifstream& in, QFont& val) {
	QString family, style_name;
	int point_size, pixel_size;

	sread(in, family);
	sread(in, style_name);
	sread(in, point_size);
	sread(in, pixel_size);

	val.setFamily(family);
	val.setStyleName(style_name);
	val.setPointSize(point_size);
	val.setPixelSize(pixel_size);
};

void swrite(std::ofstream& out, const QFont& val) {
	swrite(out, val.family());
	swrite(out, val.styleName());
	swrite(out, val.pointSize());
	swrite(out, val.pixelSize());
};

void sread(std::ifstream& in, igraph_t& val) {

	igraph_empty(&val, 0, false);

	sread(in, val.n);
	sread(in, val.directed);
	sread(in, val.from);
	sread(in, val.to);
	sread(in, val.oi);
	sread(in, val.ii);
	sread(in, val.os);
	sread(in, val.is);
};

void swrite(std::ofstream& out, const igraph_t& val) {
	swrite(out, val.n);
	swrite(out, val.directed);
	swrite(out, val.from);
	swrite(out, val.to);
	swrite(out, val.oi);
	swrite(out, val.ii);
	swrite(out, val.os);
	swrite(out, val.is);
};

void sread(std::ifstream& in, CustomMatrix& val) {
	sread(in, val.rownames_);
	sread(in, val.colnames_);
	sread(in, val.data_type_);
	sread(in, val.string_factors_);
	sread(in, val.integer_factors_);

	int size = val.colnames_.size();

	for (int i = 0; i < size; ++i) {

		QString name;
		sread(in, name);

		CustomMatrix::DataType data_type = val.data_type_[name];
		if (data_type == CustomMatrix::DataType::DoubleNumeric) {
			QVector<double>* ptr = new QVector<double>();
			sread(in, *ptr);
			val.data_[name] = ptr;
		}
		else if (data_type == CustomMatrix::DataType::IntegerNumeric || data_type == CustomMatrix::DataType::IntegerFactor) {
			QVector<int>* ptr = new QVector<int>();
			sread(in, *ptr);
			val.data_[name] = ptr;
		}
		else if (data_type == CustomMatrix::DataType::QStringFactor || data_type == CustomMatrix::DataType::QString) {
			QStringList* ptr = new QStringList();
			sread(in, *ptr);
			val.data_[name] = ptr;
		}
	}
};

void swrite(std::ofstream& out, const CustomMatrix& val) {

	swrite(out, val.rownames_);
	swrite(out, val.colnames_);
	swrite(out, val.data_type_);
	swrite(out, val.string_factors_);
	swrite(out, val.integer_factors_);

	for (const QString& colname : val.colnames_) {
		swrite(out, colname);
		CustomMatrix::DataType data_type = val.data_type_.at(colname);

		if (data_type == CustomMatrix::DataType::DoubleNumeric) {
			const QVector<double>& ref = *static_cast<QVector<double>*>(val.data_.at(colname));
			swrite(out, ref);
		}
		else if (data_type == CustomMatrix::DataType::IntegerNumeric || data_type == CustomMatrix::DataType::IntegerFactor) {
			const QVector<int>& ref = *static_cast<QVector<int>*>(val.data_.at(colname));
			swrite(out, ref);
		}
		else if (data_type == CustomMatrix::DataType::QStringFactor || data_type == CustomMatrix::DataType::QString) {
			const QStringList& ref = *static_cast<QStringList*>(val.data_.at(colname));
			swrite(out, ref);
		}
	}
};

void sread(std::ifstream& in, igraph_vector_int_t& val) {
	std::size_t size{ 0 };
	sread(in, size);
	igraph_vector_int_init(&val, size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, VECTOR(val)[i]);
	}
};

void swrite(std::ofstream& out, const igraph_vector_int_t& val) {
	std::size_t size = igraph_vector_int_size(&val);
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, VECTOR(val)[i]);
	}
};

void sread(std::ifstream& in, igraph_vector_t& val) {
	std::size_t size{ 0 };
	sread(in, size);
	igraph_vector_init(&val, size);
	for (std::size_t i = 0; i < size; ++i) {
		sread(in, VECTOR(val)[i]);
	}
};

void swrite(std::ofstream& out, const igraph_vector_t& val) {
	std::size_t size = igraph_vector_size(&val);
	swrite(out, size);
	for (std::size_t i = 0; i < size; ++i) {
		swrite(out, VECTOR(val)[i]);
	}
};