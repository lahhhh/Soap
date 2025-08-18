#pragma once

#include "VariableItem.h"

#include "CustomMatrix.h"
#include "qCustomPlot.h"

#include "DataFrame.h"

class DataFrameItem : 
	public VariableItem
{
public:
	
	G_ITEM_CONSTRUCTION(DataFrame, "DataFrame");

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

	// style : 0 - by content, 1 - by order
	// type : 0 - stop, 1 - fill, 2 - fill if no original
	static QString metadata_insert(
		const CustomMatrix& from,
		CustomMatrix& to,
		int style,
		int type,
		const QStringList& from_index,
		const QStringList& to_index
	);

	void col_slice(const Eigen::ArrayX<bool>& filter, bool in_place);
	void row_slice(const Eigen::ArrayX<bool>& filter, bool in_place);

private slots:

	void s_view();

	void s_change_data_type();
	void s_change_data_name();

	void s_choose_columns_remained();
	void s_delete_column();

	void s_rownames_as_column();
	void s_column_as_rowname();
	void s_first_row_as_colname();
	void s_first_column_as_rowname();
	void s_delete_first_column();
	void s_delete_first_row();

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

	void s_promote_to_metadata_by_rowname();
	void s_promote_to_metadata_by_metadata();

	void s_convert_to_genome_annotation();
	void s_convert_to_numeric_matrix();
	void s_convert_to_dense_double();
	void s_convert_to_bulkrna();

	void s_filter_by_feature();
};
