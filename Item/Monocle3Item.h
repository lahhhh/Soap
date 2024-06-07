#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

class Monocle3Item :
	public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(Monocle3, "Lineage Tracing", SingleCellRna, SingleCellMultiome);

	void __set_menu() override;

	void __s_update_interface() override;

	Eigen::ArrayXd extract_general_graph_ordering(const Eigen::ArrayXi& root_pr_nodes);

	Eigen::ArrayXi get_root_nodes();

	void pseudo_time_feature_plot_single_cell_rna();
	void pseudo_time_feature_plot_single_cell_multiome();

private slots:

	void s_order_cells();

	void s_feature_plot();

	void s_pseudo_time_graph();

	void s_pseudo_time_feature_plot();
};
