#pragma once

#include "VariableItem.h"

#include "EnrichmentItem.h"

#include "SingleCellRna.h"
#include "BulkRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class DifferentialAnalysisItem 
	: public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(DifferentialAnalysis, "DataFrame : DifferentialAnalysis", SingleCellRna, SingleCellAtac, DataField, ChromVAR, Cicero, BulkRna);

	void __set_menu() override;

	void __check_data() override;

	void __show_this() override;

	void __s_update_interface() override;

	void show_significant_chrom_var();

private slots:

	void s_volcano_plot();

	void s_heatmap_plot();

	void s_show_significant();

	void s_enrich_go();

	void s_enrich_kegg();

	void s_enrich_motif();

	void s_receive_enrichment(const CustomMatrix& matrix, QString name);
};


