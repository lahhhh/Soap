#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class SparseDoubleItem :
	public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(SparseDouble, "Matrix : SparseDouble", SingleCellRna, SingleCellAtac, DataField);

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

private slots:

	void s_search_row();

	void s_search_column();

	void s_search_block();

	void s_correlation_plot();
};

