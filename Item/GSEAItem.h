#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "BulkRna.h"
#include "SingleCellMultiome.h"

class GSEAItem : 
	public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(GSEA, "GSEA", SingleCellRna, SingleCellMultiome, BulkRna);

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

	void single_mountain_plot(const QString& path_name, const QString& custom_name);

	void mountain_plot_patch(
		QCustomPlot* draw_area,
		QCPLayoutGrid* layout,
		const QString& path_name,
		const QString& path_label
	);

private slots:

	void s_show_significant();

	void s_show_selected();

	void s_mountain_plot();

	void s_extract_pathway_names();
};
