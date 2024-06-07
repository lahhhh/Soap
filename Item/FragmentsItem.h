#pragma once

#include "VariableItem.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class FragmentsItem : 
    public VariableItem
{
public:

    G_ITEM_CONSTRUCTION(Fragments, "Fragments", SingleCellAtac, SingleCellMultiome);

    void __set_menu() override;

    void __s_update_interface() override;

};

