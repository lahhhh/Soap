#pragma once

#include "VariableItem.h"

#include "StringVector.h"

class StringVectorItem : 
	public VariableItem
{

public:
	G_ITEM_CONSTRUCTION(StringVector, "String Vector");

	void __set_menu() override;

	void __show_this() override;

public slots:

	void __s_update_interface() override;

private slots:

	void s_edit_manually();
	void s_edit_sort();
	void s_edit_remove_empty_elements();

	void s_convert_to_genenames();

	void s_intersect_with();
};


