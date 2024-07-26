#pragma once

#include "Identifier.h"

#include "SparseInt.h"
#include "SparseDouble.h"
#include "Cnv.h"

/*
* modified from inferCNV of the Trinity CTAT Project.  https://github.com/broadinstitute/inferCNV
*/

class InferCnvWorker 
    : public QObject
{
    Q_OBJECT

public:

    InferCnvWorker(
        const SparseInt& counts,
        const QStringList& metadata,
        const QStringList& reference,
        soap::Species species = soap::Species::Human
    ) :
        counts_(counts),
        metadata_(metadata),
        ref_group_names_(reference),
        species_(species)
    {};

    bool prepare_0();

    bool remove_insufficiently_expressed_genes_1();

    bool normalization_by_sequencing_depth_2();

    bool log_transformation_3();

    bool subtract_average_reference_5_8();

    bool smooth_6();

    bool center_7();

    bool invert_log_transform_10();

    bool compute_tumor_subcluster_11();

    bool hmm_13();

    bool generate_reports_14();

    void finalize();

    //bool bayesian_14();

    //bool bayesian_filter_15();

    //bool convert_hmm_state_16();

    //bool denoise_17();

private:

    soap::Species species_;

    SparseInt counts_;

    QStringList metadata_;
    QStringList ref_group_names_;
    QStringList obs_group_names_;

    int min_counts_per_cell_{ 100 };
    int min_cells_per_gene_{ 3 };
    int window_length_{ 101 };
    int k_nn_{ 20 };
    int random_state_{ 0 };

    double max_centered_threshold_{ 3.0 };
    double z_score_filter_{ 0.8 };
    double cut_off_{ 0.1 };

    // data

    QStringList chrs_;

    QVector<int> subclusters_;
    QVector<int> starts_;
    QVector<int> stops_;

    std::map<QString, QVector<int>> ref_grouped_cell_indices_;
    std::map<QString, QVector<int>> obs_grouped_cell_indices_;

    Eigen::MatrixXd expr_;
    Eigen::MatrixXi hmm_data_;

    // fake data

    QList<std::tuple<QString, double, int>> fake_chr_info_;

    QStringList fake_chrs_;

    std::map<QString, QVector<int>> fake_ref_grouped_cell_indices_;
    std::map<QString, QVector<int>> fake_obs_grouped_cell_indices_;

    Eigen::MatrixXd fake_expr_;

    bool scale_data_{ false };
    bool use_bounds_{ true };
    bool hmm_{ false };

    QString analysis_mode_{ "subclusters" };
    QString tumor_subcluster_partition_method_{ "leiden" };

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_cnv_ready(CNV*);
    
};

