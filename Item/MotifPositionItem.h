#pragma once

#include "VariableItem.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

#include "FootprintItem.h"
#include "ChromVARItem.h"

class MotifPositionItem : 
    public VariableItem
{

public:

    G_ITEM_CONSTRUCTION(MotifPosition, "Motif Position", SingleCellAtac, SingleCellMultiome);

    void __set_menu() override;

    void __check_data() override;

    void __show_this() override;

    void __s_update_interface() override;

    void footprint_plot_patch(
        QCustomPlot* draw_area,
        QCPLayoutGrid* layout,
        const QString& transcriptional_factor_name,
        const QString& factor_name,
        const QStringList& levels,
        const QStringList& factors,
        bool show_in_group,
        const Eigen::ArrayX<bool>& show_filter,
        const QList<QColor>& colors
    );

private slots:

    void s_show_transcriptional_factor_footprinting();
    void s_receive_transcriptional_factor_footprinting(Footprint);

    void s_difference_footprinting();
    void s_batch_task();

    void s_chromvar();
    void s_receive_chromvar(ChromVAR*);

    void s_show_motif();

    void s_show_motif_location();
    void s_show_motif_binding_location();

    void s_show_motif_matrix();

    void s_multiple_footprint_plot();

    void s_show_motif_overlapping_rate();

    void s_show_typical_binding_sequence();
};

