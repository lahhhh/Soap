#pragma once

#include "Identifier.h"

#include <QMainWindow>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>

#include "qCustomPlot.h"
#include "PlotsSuite.h"

class PlotWindow : public QMainWindow
{
	Q_OBJECT
public:
	PlotWindow(PlotsSuite* plot_suite, QString title, int width, int height, QWidget* parent = nullptr);

	PlotWindow(const QImage& pic, QString title, QWidget* parent = nullptr);

	QWidget* main_interface_;

	QVBoxLayout* main_layout_;

	QCustomPlot* plot_;


	static void show_plot(PlotsSuite* plot_suite, QString title, int width, int height, QWidget* parent = nullptr);

	static void show_plot(const QImage& pic, QString title, QWidget* parent = nullptr);

};
