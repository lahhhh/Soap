#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

#include "FeatureHandler.h"

class EmbeddingItem : 
	public VariableItem
{
public:
	
	G_ITEM_CONSTRUCTION(Embedding, "Matrix : Embedding", DataField, SingleCellRna, SingleCellAtac, ChromVAR, Cicero, BulkRna);

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

	void show(
		const Eigen::ArrayXd& data,
		const QString& title,
		bool scale,
		const QString& legend_title
	);

	void show(
		const QStringList& data,
		const QString& name
	);

private slots:

	void s_raw_plot();

	void s_feature_plot();

	void s_highlight();

	void s_multiple_numeric_feature_plot();

	void s_chromvar_plot();

	void s_cicero_plot();

	void s_view();

	void s_show_depth_correlation();

	void s_deviation_plot();
};
