#pragma once

#include "VariableItem.h"

#include "DifferentialAnalysis.h"
#include "GeneName.h"

class EnrichmentItem : 
	public VariableItem
{
public:	
	G_ITEM_CONSTRUCTION(Enrichment, "DataFrame : Enrichment", DifferentialAnalysis, GeneName);

	void __set_menu() override;

	void __show_this() override;

	void __s_update_interface() override;

	void barplot(int path_number = 10);

	void barplot(const QStringList& pathways);

private slots:

	void s_show_significant();

	void s_barplot_show_selected();

	void s_barplot_default();

	void s_barplot();
};
