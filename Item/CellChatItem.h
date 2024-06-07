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

	void s_show_interaction_weights();

	void s_show_interaction_probability();

	void s_show_interaction_p_value();

	void s_show_pathway_probability();
};
