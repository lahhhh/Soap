#include "CellChat.h"

CellChat::CellChat(
    const QString & identity, 
    const QStringList & ligand_receptor_index, 
    const QStringList & pathway_index, 
    const QStringList & levels, 
    const Eigen::ArrayXXi & counts, 
    const Eigen::ArrayXXd & weights,
    const Eigen::ArrayX<Eigen::ArrayXXd> & ligand_receptor_probability, 
    const Eigen::ArrayX<Eigen::ArrayXXd> & ligand_receptor_p_value, 
    const QMap<QString, Eigen::ArrayXXd> & pathway_probability
):
    identity_(identity), 
    ligand_receptor_index_(ligand_receptor_index), 
    pathway_index_(pathway_index), 
    levels_(levels), 
    counts_(counts), 
    weights_(weights), 
    ligand_receptor_probability_(ligand_receptor_probability), 
    ligand_receptor_p_value_(ligand_receptor_p_value),  
    pathway_probability_(pathway_probability)
{

}