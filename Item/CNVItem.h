#pragma once

#include "VariableItem.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

class CNVItem : 
	public VariableItem
{
public:

	G_ITEM_CONSTRUCTION(CNV, "CNV Detection", SingleCellRna, SingleCellMultiome);

	void draw_cnv(
		const Eigen::MatrixXd& cnv_mat,
		const std::vector<std::tuple<QString, int, int>>& chromosome_location,
		const std::vector<std::tuple<QString, int, int>>& metadata_location
	);

	void __set_menu() override;

	void __s_update_interface() override;

private slots:

	void s_view();
};
