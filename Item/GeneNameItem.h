#pragma once

#include "VariableItem.h"

#include "SingleCellMultiome.h"

class GeneNameItem : 
	public VariableItem
{
public:
	
	G_ITEM_CONSTRUCTION(GeneName, "Gene Name", Footprint, Pando);

	void __set_menu() override;

	void __check_data() override;

	void __show_this() override;

public slots:

	void __s_update_interface() override;

private slots:

	void s_enrich_go();

	void s_enrich_kegg();

	void s_receive_enrichment(const CustomMatrix& matrix, QString name);

	void s_edit_manually();
	void s_edit_sort();
	void s_edit_unique();
	void s_remove_empty_elements();
	void s_remove_elements_by_regular_expression();
	void s_remove_elements_by_chromosome_location();

	void s_convert_to_string_vector();
};


