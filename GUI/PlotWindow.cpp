#include "PlotWindow.h"

PlotWindow::PlotWindow(PlotsSuite * plot_suite, QString title, int width, int height, QWidget *parent) : 
    QMainWindow(parent)
{
    this->main_interface_ = new QWidget(this);
    this->setCentralWidget(this->main_interface_);

    this->main_layout_ = new QVBoxLayout;

    this->main_interface_->setLayout(this->main_layout_);

    this->plot_ = new QCustomPlot(this);

    this->plot_->plotLayout()->clear();

    this->plot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    this->plot_->setBackgroundScaledMode(Qt::IgnoreAspectRatio);

    this->plot_->setBackground(plot_suite->plots_[plot_suite->current_plot_id_]->toPixmap());

    this->main_layout_->addWidget(this->plot_);
    this->main_layout_->setContentsMargins(0, 0, 0, 0);
    this->resize(width, height);
    this->setAttribute(Qt::WA_DeleteOnClose);
    this->setWindowTitle(title);
    this->setWindowIcon(QIcon(FILE_SOAP_ICON_JPG));
    this->move(100, 100);

    this->plot_->replot();

    this->show();

}

PlotWindow::PlotWindow(const QImage& pic, QString title, QWidget* parent) : QMainWindow(parent) {
    this->main_interface_ = new QWidget(this);
    this->setCentralWidget(this->main_interface_);

    this->main_layout_ = new QVBoxLayout;

    this->main_interface_->setLayout(this->main_layout_);

    this->plot_ = new QCustomPlot(this);

    this->plot_->plotLayout()->clear();

    this->plot_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    this->plot_->setBackgroundScaledMode(Qt::IgnoreAspectRatio);

    this->plot_->setBackground(QPixmap::fromImage(pic));

    this->main_layout_->addWidget(this->plot_);
    this->main_layout_->setContentsMargins(0, 0, 0, 0);
    this->resize(pic.width(), pic.height());
    this->setAttribute(Qt::WA_DeleteOnClose);
    this->setWindowTitle(title);
    this->setWindowIcon(QIcon(FILE_SOAP_ICON_JPG));
    this->move(100, 100);

    this->plot_->replot();

    this->show();
};


void PlotWindow::show_plot(PlotsSuite * plot_suite, QString title, int width, int height, QWidget *parent){
    PlotWindow * plot = new PlotWindow(plot_suite, title, width, height, parent);
};

void PlotWindow::show_plot(const QImage& pic, QString title, QWidget* parent) {
    PlotWindow* plot = new PlotWindow(pic, title, parent);
};
