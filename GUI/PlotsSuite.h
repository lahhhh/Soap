#pragma once

#include "Identifier.h"

#include "qCustomPlot.h"
#include "GraphSettings.h"

class PlotsSuite : public QObject
{
    Q_OBJECT
public:
    PlotsSuite();
    ~PlotsSuite();

    QCustomPlot * current_plot_;
    QList<QCustomPlot * > plots_;

    GraphSettings graph_settings_;


    int current_plot_id_;
    int maximum_plot_id_;

    void update(QCustomPlot *);

    void clear();

    void prepare();

public slots:
    void s_setting_activate(bool activated);

signals:
    void x_plot_prepared();

};
