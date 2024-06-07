#pragma once

#include "VariableItem.h"

#include "NumericMatrix.h"

class NumericMatrixItem 
	: public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(NumericMatrix, "Numeric Matrix");

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

};

