#pragma once

#include "VariableItem.h"

#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class GenomicRangeItem :
    public VariableItem
{
public:

    G_ITEM_CONSTRUCTION(GenomicRange, "Genomic Range", SingleCellAtac, SingleCellMultiome);

    void __set_menu() override;

    void __show_this() override;

    void __s_update_interface() override;

    static QStringList ucsc_to_ensembl(const QStringList& ucsc_name);

    static QStringList ensembl_to_ucsc(const QStringList& ensembl_name);

private slots:

    void s_recalculate_counts();

    void s_change_sequence_name_style();

    void s_recalculated_counts(SparseInt*);
};

