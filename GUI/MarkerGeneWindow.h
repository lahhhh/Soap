#pragma once

#include "Identifier.h"

#include <QMainWindow>
#include <QTableWidget>
#include <QHBoxLayout>
#include <QPushButton>

#include "CellMarkerTreeWidget.h"

class MarkerGeneWindow : 
    public QMainWindow
{
    Q_OBJECT
public:
    explicit MarkerGeneWindow(
        const QString & cell_type, 
        const QStringList & marker_names, 
        CellMarkerTreeWidget::MarkerWidgetType type,
        QWidget *parent = nullptr);

    void set_matrix_table();

    QWidget * main_interface_;

    QHBoxLayout main_layout_;

    QTableWidget * matrix_table_;

    QStringList marker_names_;

public slots:

    void s_marker_clicked1();

    void s_marker_clicked2();

    void s_marker_clicked3();

signals:

    void x_show_gene(QString gene_name);

    void x_to_active(QString gene_name);

    void x_to_alternative(QString gene_name);
};

