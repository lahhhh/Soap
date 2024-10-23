#pragma once

#include "Identifier.h"

#include <unordered_map>

#include "SignalEmitter.h"
#include "SingleCellMultiome.h"

class LoadFragmentsWorker: 
    public QObject
{
    Q_OBJECT
public:

    LoadFragmentsWorker(
        const QStringList& barcodes,
        const QStringList& fragments_files
    );

    QStringList barcodes_;

    QStringList fragments_files_;

    std::unordered_map<QString, int> barcodes_index_;

    int n_cell_;

    std::unique_ptr<Fragments> fragments_;

    void create_index();

    bool load_fragments();

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_fragments_ready(Fragments*);
    
};

