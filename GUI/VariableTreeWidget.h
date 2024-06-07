#pragma once

#include "Identifier.h"

#include <QTreeWidget>
#include <QCursor>

#include "qCustomPlot.h"
#include "InformationTextBrowser.h"
#include "PlotsSuite.h"
#include "SignalEmitter.h"

class VariableTreeWidget : 
    public QTreeWidget
{
    Q_OBJECT
public:
    VariableTreeWidget(
        QWidget * parent, 
        PlotsSuite *draw_suite,
        InformationTextBrowser * information_area, 
        SignalEmitter * signal_emitter
    );

    PlotsSuite * draw_suite_;
    InformationTextBrowser * information_area_;
    SignalEmitter * signal_emitter_;

private slots:

    void s_item_double_clicked(QTreeWidgetItem *, int column);

    void s_item_right_clicked(const QPoint&);
};

