#pragma once

#include "VariableItem.h"

#include "SparseIntItem.h"
#include "VelocityEstimateItem.h"
#include "ScveloEstimateItem.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

class VelocytoBaseItem : 
	public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(VelocytoBase, "Lineage Tracing : Velocyto", SingleCellRna, SingleCellMultiome);

	void __set_menu() override;

	void __check_data() override;

private slots:

	void s_estimate_velocity();
	void s_receive_estimate(VelocityEstimate*);

	void s_estimate_scvelo();
	void s_receive_scvelo(ScveloEstimate*);
};

