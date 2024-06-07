#pragma once

#include "Identifier.h"

#include <QTreeWidget>

class CellMarkerTreeWidget : 
    public QTreeWidget
{

    Q_OBJECT

public:

    enum class MarkerWidgetType { MarkerWidget, ActiveWidget, AlternativeWidget };

    CellMarkerTreeWidget(CellMarkerTreeWidget::MarkerWidgetType type, QWidget * parent);

    void set_types(const QMap<QString, QStringList>& type_to_markers);

    void safe_clear();

    QMap<QString, QStringList> type_to_markers_;

    QStringList types_;

    QList<QStringList> markers_;

    CellMarkerTreeWidget::MarkerWidgetType type_;
    

public slots:

    void s_type_clicked_1();

    void s_type_clicked_2();

    void s_type_clicked_3();

    void s_type_clicked_4();

    void QTreeWidgetItemDoubleClicked(QTreeWidgetItem *, int column);

    void s_receive_type(const QString & type, const QStringList & markers);

    void s_receive_marker(const QString & gene);

signals:

    void x_show_gene(QString gene);

    void x_show_type(QString type, QStringList markers);

    void x_show_complex(QString title, QStringList markers);

    void x_type_to_active(QString type, QStringList markers);

    void x_type_to_alternative(QString type, QStringList markers);

    void x_gene_to_active(QString gene);

    void x_gene_to_alternative(QString gene);
};
