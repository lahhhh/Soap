#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class SparseIntItem :
	public VariableItem
{

public:
	G_ITEM_CONSTRUCTION(SparseInt, "Matrix : SparseInt", SingleCellRna, SingleCellAtac, DataField, VelocytoBase);

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

public slots:

	void __s_delete_this() override;

private slots:

	void s_search_row();
	void s_search_column();
	void s_search_block();

	void s_extract_as_single_cell_rna();
	void s_receive_single_cell_rna(SingleCellRna*);
};

