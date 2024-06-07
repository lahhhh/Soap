#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>

#include "Monocle3.h"
#include "CustomPlot.h"

class Monocle3ChooseRootDialog
	: public QDialog
{

public:
	Monocle3ChooseRootDialog(Monocle3* data);

	Monocle3* data_{ nullptr };

	QVector<int> choosed_nodes_;

	QPushButton* finish_button_{ nullptr };
	QPushButton* cancel_button_{ nullptr };

	QCustomPlot* draw_area_{ nullptr };
	QCPAxisRect* axis_rect_{ nullptr };
	QCPGraph* move_graph_{ nullptr };
	std::map<int, QCPGraph*> choosed_graph_;

	bool is_accepted_{ false };

	static QVector<int> get_response(Monocle3* data);

private slots:

	void accept();
	void reject();

	void s_plot_mouse_press(QMouseEvent*);

	void s_plot_mouse_move(QMouseEvent*);
	
};

