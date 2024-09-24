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
};
