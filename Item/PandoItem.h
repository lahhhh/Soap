#pragma once

#include "VariableItem.h"

#include "SingleCellMultiome.h"

class PandoItem :
	public VariableItem
{
	
public:
	G_ITEM_CONSTRUCTION(Pando, "Pando", SingleCellMultiome);

	void __set_menu() override;

	void __check_data() override;

	void __show_this() override;

	void __s_update_interface() override;

private slots:

	void s_show_significant();

	void s_abstract();

	void s_show_factor();
	void s_show_factor_coef();

	void s_show_gene_name();
	void s_show_gene_coef();

	void s_show_factor_gene_strength();

	void s_extract_gene_name();
};


