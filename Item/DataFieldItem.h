#pragma once

#include "VariableItem.h"

#include "SingleCellMultiome.h"

#include "SparseIntItem.h"
#include "SparseDoubleItem.h"
#include "EmbeddingItem.h"
#include "DifferentialAnalysisItem.h"

class DataFieldItem : 
    public VariableItem
{
public:

    G_ITEM_CONSTRUCTION(DataField, "Field", SingleCellMultiome);

    void __set_menu() override;

    void __check_data() override;

    SparseIntItem* counts();
    SparseDoubleItem* normalized();
    EmbeddingItem* pca();
    EmbeddingItem* tsne();
    EmbeddingItem* umap();
    EmbeddingItem* harmony();

    void rna_view_quality();
    void atac_view_quality();

    bool rna_re_estimate_quality_parameters();
    bool atac_re_estimate_quality_parameters();

private slots:

    void s_filter_by_feature();

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

    void s_pca_default();
    void s_pca_custom();
    void s_receive_pca(Eigen::MatrixXd, QVector<double>);

    void s_svd_default();
    void s_svd_custom();
    void s_receive_svd(Eigen::MatrixXd, QVector<double>);

    void s_tsne_default();
    void s_tsne_custom();
    void s_receive_tsne(Eigen::MatrixXd);

    void s_umap_default();
    void s_umap_custom();
    void s_receive_umap(Eigen::MatrixXd);

    void s_harmony();
    void s_receive_harmony(Eigen::MatrixXd);

    void s_log_normalize_default();
    void s_log_normalize_custom();
    void s_receive_normalize(SparseDouble*);

    void s_tfidf_default();
    void s_receive_tfidf(SparseDouble*);

    void s_view_quality();

    void s_distribution_plot();

    void s_bubble_plot();

    void s_heatmap_plot();    
    void s_cell_heatmap_plot();

    void s_correlation_heatmap_plot();

    void s_find_deg();
    void s_find_dap();
    void s_find_dag();

    void s_receive_differential_analysis(DifferentialAnalysis da, QString name);

    void s_re_estimate_quality_parameters();
};



