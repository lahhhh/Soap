#pragma once

#include "VariableItem.h"

#include "DifferentialAnalysisItem.h"
#include "EmbeddingItem.h"

#include "MotifPosition.h"

class ChromVARItem :
	public VariableItem
{
public:
	G_ITEM_CONSTRUCTION(ChromVAR, "ChromVAR", MotifPosition);

	void __set_menu() override;

	void __check_data() override;

private slots:

	void s_differential_analysis();
	void s_receive_differential_analysis(DifferentialAnalysis da, QString name);

	void s_umap();
	void s_receive_umap(Eigen::MatrixXd);

	void s_embedding_plot();
};
