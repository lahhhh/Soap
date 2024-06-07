#pragma once

#include "VariableItem.h"

#include "CustomPlot.h"
#include "EmbeddingItem.h"
#include "MetadataItem.h"
#include "SparseDoubleItem.h"
#include "SparseIntItem.h"
#include "GSEAItem.h"
#include "CellChatItem.h"
#include "CNVItem.h"
#include "VelocytoBaseItem.h"
#include "NoteItem.h"
#include "DifferentialAnalysisItem.h"
#include "Monocle3Item.h"

#include "SingleCellRna.h"

class SingleCellRnaItem 
	: public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(SingleCellRna, "SingleCellRna");

	void __set_menu() override;

	void __check_data() override;

	G_QUICK_ACCESS_ITEM(Metadata, metadata);
	G_QUICK_ACCESS_ITEM_TYPE(SparseInt, counts, Counts);
	G_QUICK_ACCESS_ITEM_TYPE(SparseDouble, normalized, Normalized);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, pca, Pca);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, umap, Umap);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, tsne, Tsne);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, harmony, Harmony);
	G_QUICK_ACCESS_ITEM(VelocytoBase, velocyto_base);

	void slice(const Eigen::ArrayX<bool>& row_slice, const Eigen::ArrayX<bool>& col_slice);

	void col_slice(const Eigen::ArrayX<bool>& col_slice);

	std::pair<Eigen::ArrayX<bool>, Eigen::ArrayX<bool> > get_parameter_filter(const QStringList& parameters);

public slots:

	void __s_rename() override;

	void __s_delete_this() override;

private slots:

	void s_fast_annotation();
	void s_receive_fast_annotation(
		QStringList main_type,
		QStringList sub_type,
		QString main_type_name,
		QString sub_type_name
	);

	void s_recalculate_quality_parameters();

	void s_bubble_plot();

	void s_distribution_plot();

	void s_sample();

	void s_statistics();

	void s_set_random_state();

	void s_set_species();

	void s_add_metadata();
	void s_combine_existed_metadata();
	void s_edit_metadata();

	void s_view_quality();

	void s_log_normalize_default();
	void s_log_normalize_custom();
	void s_receive_normalize(SparseDouble*);

	void s_filter_by_parameters();
	void s_filter_by_features();

	void s_umap_default();
	void s_umap_custom();
	void s_receive_umap(Eigen::MatrixXd mat);

	void s_tsne_default();
	void s_tsne_custom();
	void s_receive_tsne(Eigen::MatrixXd mat);

	void s_pca_default();
	void s_pca_custom();
	void s_receive_pca(Eigen::MatrixXd, QVector<double>);

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

	void s_scrublet();
	void s_receive_scrublet(const Eigen::ArrayXd& original_score, const Eigen::ArrayXd& simulate_score);
	void s_set_scrublet_threshold();

	void s_monocle3();
	void s_receive_monocle3(Monocle3*);

	void s_velocyto();
	void s_receive_velocyto(VelocytoBase* velocyto_base);

	void s_find_deg();
	void s_receive_differential_analysis(DifferentialAnalysis da, QString name);

	void s_gsea();
	void s_receive_gsea(GSEA gsea);

	void s_cellchat_default();
	void s_receive_cellchat(CellChat cellchat);

	void s_scicnv();
	void s_receive_scicnv(CNV* cnv);

	void s_infercnv();
	void s_receive_infercnv(CNV* cnv);

	void s_integrate();
	void s_receive_integrated_data(SingleCellRna* data, QList<const SingleCellRna*> items);

	void s_harmony();
	void s_receive_harmony(Eigen::MatrixXd);

};
