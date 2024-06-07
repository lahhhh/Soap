#include "PlotsSuite.h"

PlotsSuite::PlotsSuite() : 
    current_plot_id_(0), 
    maximum_plot_id_(0)
{
}

void PlotsSuite::prepare(){
    QCustomPlot * plot = new QCustomPlot();

    plot->clearPlottables();
    plot->clearGraphs();
    plot->plotLayout()->clear();

    this->current_plot_ = plot;
    this->current_plot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    this->plots_.append(plot);

    emit x_plot_prepared();
}

void PlotsSuite::update(QCustomPlot * plot){
    this->current_plot_ = plot;
    this->current_plot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    this->plots_.append(plot);

    this->maximum_plot_id_ ++;
    this->current_plot_id_ = this->maximum_plot_id_;

    emit x_plot_prepared();
}

PlotsSuite::~PlotsSuite(){
    for(auto ptr : this->plots_){
        delete ptr;
    }
};

void PlotsSuite::clear(){
    for(auto ptr : this->plots_){
        delete ptr;
    }
    this->plots_.clear();

    this->current_plot_id_ = 0;
    this->maximum_plot_id_ = 0;
    prepare();
};

void PlotsSuite::s_setting_activate(bool activated){
    this->graph_settings_.set_activated(activated);
};
