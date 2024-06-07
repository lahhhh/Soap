#pragma once

#include "VariableItem.h"

#include "BulkRna.h"

class DenseDoubleItem :
	public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(DenseDouble, "Matrix : DenseDouble", BulkRna);

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

	void update_shape();

private slots:

	void s_search_row();

	void s_search_column();

	void s_search_block();

	void s_correlation_plot();
};



