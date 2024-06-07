#pragma once

#include "VariableItem.h"

#include "VelocytoBase.h"

#include "CustomPlot.h"

class VelocityEstimateItem : 
	public VariableItem
{

public:

	G_ITEM_CONSTRUCTION(VelocityEstimate, "Velocity", VelocytoBase);

	void __set_menu() override;

	void show_embedding_velocity(int mode);

	std::pair<QCustomPlot*, QCPAxisRect*> draw_feature_plot(
		const Eigen::MatrixXd& embedding,
		const QStringList& embedding_names,
		const QStringList& graph_settings
	);

private slots:

	void s_show_embedding_velocity_by_stream();
	void s_show_embedding_velocity_by_grid();

	void s_receive_velocity_grid(VELO_GRID_PLOT_ELEMENTS);

	void s_receive_velocity_stream(STREAM_PLOT_ELEMENTS);
};

