#pragma once

#include "VariableItem.h"

#include "CustomPlot.h"

#include "EmbeddingItem.h"
#include "MetadataItem.h"
#include "DenseDoubleItem.h"
#include "DenseIntItem.h"
#include "GSEAItem.h"
#include "NoteItem.h"
#include "DifferentialAnalysisItem.h"

class BulkRnaItem
	: public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(BulkRna, "RNA-Seq");

	G_QUICK_ACCESS_ITEM(Metadata, metadata);
	G_QUICK_ACCESS_ITEM_TYPE(DenseInt, counts, Counts);
	G_QUICK_ACCESS_ITEM_TYPE(DenseDouble, normalized, Normalized);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, pca, Pca);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, umap, Umap);
	G_QUICK_ACCESS_ITEM_TYPE(Embedding, tsne, Tsne);

	void __set_menu() override;

	void __check_data() override;

public slots:

	void __s_rename() override;

	void __s_delete_this() override;

private slots:

	void s_statistics();

	void s_normalize1();
	void s_fpkm();
	void s_tpm();
	void s_receive_normalize(DenseDouble);

	void s_pca_default();
	void s_pca_custom();
	void s_receive_pca(Eigen::MatrixXd, QVector<double>, QVector<double>);

	void s_umap_default();
	void s_umap_custom();
	void s_receive_umap(Eigen::MatrixXd mat);

	void s_tsne_default();
	void s_tsne_custom();
	void s_receive_tsne(Eigen::MatrixXd mat);

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

	void s_find_deg();
	void s_receive_differential_analysis(DifferentialAnalysis da, QString name);

	void s_gsea();
	void s_receive_gsea(GSEA gsea);

	void s_integrate();
	void s_receive_integrated_data(BulkRna* data, QList<const BulkRna*> items);

	void s_set_species();

	void s_set_random_state();
};

