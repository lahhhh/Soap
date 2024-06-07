#pragma once

#include "VariableItem.h"

#include "qCustomPlot.h"
#include "EmbeddingItem.h"
#include "MetadataItem.h"
#include "SparseDoubleItem.h"
#include "SparseIntItem.h"
#include "NoteItem.h"
#include "DifferentialAnalysisItem.h"
#include "GenomicRangeItem.h"
#include "MotifPositionItem.h"
#include "CoverageTrackItem.h"
#include "FragmentsItem.h"
#include "CiceroItem.h"
#include "Monocle3Item.h"

#include "CoveragePlotWorker.h"
#include "AtacLandscapePlotWorker.h"

#include "SingleCellAtac.h"

class SingleCellAtacItem
	: public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(SingleCellAtac, "scATAC");

	void __set_menu() override;

	void __check_data() override;

	G_QUICK_ACCESS_ITEM(Metadata, metadata);
	G_QUICK_ACCESS_ITEM_TYPE(SparseInt, counts, Counts);
	G_QUICK_ACCESS_ITEM_TYPE(SparseInt, gene_activity_counts, GeneActivity);
	G_QUICK_ACCESS_ITEM_TYPE(SparseDouble, normalized, Normalized);
	G_QUICK_ACCESS_ITEM_TYPE(SparseDouble, gene_activity_normalized, GeneActivity);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, pca, Pca);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, umap, Umap);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, tsne, Tsne);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, harmony, Harmony);
	G_QUICK_ACCESS_ITEM(Fragments, fragments);
	G_QUICK_ACCESS_ITEM(MotifPosition, motif_position);
	G_QUICK_ACCESS_ITEM(Cicero, cicero);

	void col_slice(const Eigen::ArrayX<bool>& col_slice);

	Eigen::ArrayX<bool> get_parameter_filter(const QStringList& parameters);

	void update_quality_control_information();

public slots:

	void __s_rename() override;

	void __s_delete_this() override;

private slots:

	void s_statistics();

	void s_calculate_quality_metrics();
	void s_receive_qc(QMap<QString, QList<double>>);

	void s_filter_by_parameters();
	void s_filter_by_features();
	void s_filter_by_quality_metrics();

	void s_show_quality_matrics();
	void s_show_fragments_length_distribution();

	void s_tfidf_default();
	void s_receive_tfidf(SparseDouble*);

	void s_normalize_gene_activity();
	void s_receive_normalized_gene_activity(SparseDouble*);

	void s_svd_default();
	void s_svd_custom();
	void s_receive_svd(Eigen::MatrixXd, QVector<double>);

	void s_umap_default();
	void s_umap_custom();
	void s_receive_umap(Eigen::MatrixXd);

	void s_tsne_default();
	void s_tsne_custom();
	void s_receive_tsne(Eigen::MatrixXd);

	void s_leiden_default();
	void s_leiden_custom();
	void s_receive_leiden(QVector<int> cluster);

	void s_louvain_default();
	void s_louvain_custom();
	void s_receive_louvain(std::vector<int> cluster);

	void s_modified_louvain_default();
	void s_modified_louvain_custom();
	void s_receive_modified_louvain(std::vector<int> cluster);

	void s_smart_local_moving_default();
	void s_smart_local_moving_custom();
	void s_receive_slm(std::vector<int> cluster);

	void s_find_dap();
	void s_receive_differential_analysis(DifferentialAnalysis da, QString name);

	void s_tss_plot();
	void s_receive_tss_plot_data(Eigen::ArrayXXd tss_matrix, QVector<double> tss_vector);

	void s_coverage_plot();
	void s_receive_coverage_plot_data(COVERAGE_PLOT_ELEMENTS);

	void s_create_coverage_track();
	void s_receive_coverage_track(CoverageTrack* d);

	void s_call_peaks_by_macs();
	void s_receive_macs_peaks(GenomicRange);

	void s_sample();

	void s_set_random_state();

	void s_set_species();

	void s_add_new_metadata();
	void s_combine_existed_metadata();
	void s_edit_metadata();

	void s_show_atac_landscape();
	void s_receive_atac_landscape_plot(ATAC_LANDSCAPE_PLOT_ELEMENTS);

	void s_calculate_gene_activity();
	void s_receive_gene_activity(SparseInt*);

	void s_fast_annotation();
	void s_receive_fast_annotation(
		QStringList main_type,
		QStringList sub_type,
		QString main_type_name,
		QString sub_type_name
	);

	void s_find_motifs();
	void s_receive_motif_location(MotifPosition);

	void s_integrate();
	void s_receive_integrated_data(SingleCellAtac* data, QList<const SingleCellAtac*> items);

	void s_load_fragments();
	void s_receive_loaded_fragments(Fragments*);

	void s_monocle3();
	void s_receive_monocle3(Monocle3*);
};

