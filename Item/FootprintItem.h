#pragma once

#include "VariableItem.h"

#include "MotifPosition.h"

class FootprintItem : 
    public VariableItem
{
public:

    G_ITEM_CONSTRUCTION(Footprint, "Footprint", MotifPosition);

    void __set_menu() override;

    void __check_data() override;

    void __show_this() override;

private slots:

    void s_show();

    void s_receive_correlated_gene(QString, QStringList);
};

