#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "BulkRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class MetadataItem : 
	public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(Metadata, "Matrix : Meta", SingleCellRna, SingleCellAtac, SingleCellMultiome, BulkRna);

	void __set_menu() override;

	void __show_this() override;

public slots:

	void __s_update_interface() override;

	void __s_delete_this() override;

private slots:

	void s_extract();

	void s_rownames_as_column();
	void s_delete_metadata();
	void s_change_data_type();
	void s_change_data_name();

	void s_slice_data();
	void s_map_values();
	void s_replace_first();
	void s_replace_last();
	void s_replace_all();
	void s_remove_first();
	void s_remove_last();
	void s_remove_all();
	void s_remove_prefix();
	void s_remove_suffix();
	void s_remove_until_first();
	void s_remove_until_last();
	void s_remove_from_first();
	void s_remove_from_last();
	void s_add_prefix();
	void s_add_suffix();

	void s_view();

	void s_compare_plot();
};
