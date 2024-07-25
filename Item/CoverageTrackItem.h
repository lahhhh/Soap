#pragma once

#include "VariableItem.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class CoverageTrackItem 
    : public VariableItem
{
public:

    G_ITEM_CONSTRUCTION(CoverageTrack, "Coverage Track", SingleCellAtac, SingleCellMultiome);

    void __set_menu() override;

    void __show_this() override;

    void __s_update_interface() override;

    void s_show_coverage();

};

