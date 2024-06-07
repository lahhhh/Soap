#include "cellchatworker.h"
#include "CellChatDatabase.h"

CellChatWorker::CellChatWorker(const QString& identity, const SparseDouble& normalized, soap::Species species, 
    const QStringList& metadata, double minimum_percentage, const QString& p_adjust_method, int random_state,
    int n_boot, int minimum_cell_number, const QString& annotation_type) :
    identity_(identity),
    normalized_(normalized), 
    species_(species), 
    metadata_(metadata), 
    minimum_percentage_(minimum_percentage), 
    p_adjust_method_(p_adjust_method), 
    random_state_(random_state), 
    n_boot_(n_boot), 
    minimum_cell_number_(minimum_cell_number), 
    annotation_type_(annotation_type)
{

}

void CellChatWorker::run(){
    CellChatDatabase db;

    emit x_cellchat_ready(
        db.cellchat(
            this->identity_, 
            &this->normalized_,
            this->species_, 
            this->metadata_, 
            this->minimum_percentage_, 
            this->p_adjust_method_, 
            this->random_state_, 
            this->n_boot_, 
            this->minimum_cell_number_,
            this->annotation_type_)
    );
    G_TASK_END;
}
