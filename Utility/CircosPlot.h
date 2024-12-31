#pragma once

#include "Identifier.h"

#include "GraphSettings.h"

void circos_plot(
	QCustomPlot* draw_area,
	QCPAxisRect* axis_rect,
	const QStringList& levels,
	const QList<QColor>& colors,
	const Eigen::MatrixXd& data,
	bool color_by_source = true
);