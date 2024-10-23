#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

class CellChatItem :
	public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(CellChat, "Cell-Cell Interaction : CellChat", SingleCellRna, SingleCellMultiome);

	void __set_menu() override;

private slots:

	void s_show_interaction_numbers();
	void s_interaction_number_heatmap();
	void s_interaction_number_circos();

	void s_show_interaction_weights();
	void s_interaction_weight_heatmap();
	void s_interaction_weight_circos();

	void s_show_interaction_probability();
	void s_interaction_probability_heatmap();
	void s_interaction_probability_circos();

	void s_show_interaction_p_value();
	void s_interaction_p_value_heatmap();
	void s_interaction_p_value_circos();

	void s_show_pathway_probability();
	void s_pathway_probability_heatmap();
	void s_pathway_probability_circos();

	void s_show_interaction_summary();
	void s_extract_interaction_summary();

	void s_show_pathway_summary();
	void s_extract_pathway_summary();
};
