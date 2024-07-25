#pragma once

#include "VariableItem.h"

#include "EmbeddingItem.h"
#include "DifferentialAnalysisItem.h"
#include "CoveragePlotWorker.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

class CiceroItem :
	public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(Cicero, "Cicero", SingleCellAtac, SingleCellMultiome);

	void __set_menu() override;

	void __check_data() override;

	void __s_update_interface() override;

private slots:

	void s_differential_analysis();
	void s_receive_differential_analysis(DifferentialAnalysis da, QString name);

	void s_show_group();

	void s_show_group_transcription_factor_possibility();

	void s_umap();
	void s_receive_umap(Eigen::MatrixXd);

	void s_show_ccan_coverage();
	void s_receive_coverage_plot_data(COVERAGE_PLOT_ELEMENTS);

	void s_scan_downstream_ccan();

	void s_ccan_embedding_feature_plot();

	void s_reset_threshold();
};
