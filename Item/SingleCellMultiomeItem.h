#pragma once

#include "VariableItem.h"

#include <string>
#include <set>

#include "qCustomPlot.h"

#include "SingleCellMultiome.h"

#include "EmbeddingItem.h"
#include "MetadataItem.h"
#include "SparseDoubleItem.h"
#include "SparseIntItem.h"
#include "GSEAItem.h"
#include "CellChatItem.h"
#include "CNVItem.h"
#include "DataFieldItem.h"
#include "GenomicRangeItem.h"
#include "MotifPositionItem.h"
#include "CoverageTrackItem.h"
#include "FragmentsItem.h"
#include "VelocytoBaseItem.h"
#include "NoteItem.h"
#include "PandoItem.h"
#include "Monocle3Item.h"
#include "CiceroItem.h"

#include "CoveragePlotWorker.h"
#include "AtacLandscapePlotWorker.h"

class SingleCellMultiomeItem :
	public VariableItem
{

public:
	G_ITEM_CONSTRUCTION(SingleCellMultiome, "SingleCellMultiome");

	G_QUICK_ACCESS_ITEM(Metadata, metadata);
	G_QUICK_ACCESS_ITEM(Fragments, fragments);
	G_QUICK_ACCESS_ITEM(MotifPosition, motif_position);
	G_QUICK_ACCESS_ITEM(Cicero, cicero);
	G_QUICK_ACCESS_ITEM(VelocytoBase, velocyto_base);

	DataFieldItem* create_field_item(DataField::DataType type, const QString& name);
	DataFieldItem* get_field_item(DataField::DataType type);

	DataFieldItem* get_rna_field();
	DataFieldItem* get_atac_field();
	DataFieldItem* get_trans_field();

	void __set_menu() override;

	void __check_data() override;

	void update_quality_control_information(DataField::DataType type);

	void rna_update_quality_control_information();

	void atac_update_quality_control_information();

	void slice(const Eigen::ArrayX<bool>& row_slice, const Eigen::ArrayX<bool>& col_slice, bool in_place);

	void col_slice(const Eigen::ArrayX<bool>& filter, bool in_place);

	std::pair<Eigen::ArrayX<bool>, Eigen::ArrayX<bool> > get_parameter_filter(const QStringList& parameters);

public slots:

	void __s_rename() override;

	void __s_delete_this() override;

private slots:

	void s_fast_annotation();
	void s_receive_fast_annotation(
		QStringList main_type,
		QStringList sub_type,
		QString main_type_name,
		QString sub_type_name
	);

	void s_statistics();

	void s_load_fragments();
	void s_receive_loaded_fragments(Fragments*);

	void s_find_motifs();
	void s_receive_motif_location(MotifPosition);

	void s_cicero();
	void s_receive_cicero(Cicero*);

	void s_calculate_gene_activity();
	void s_receive_gene_activity(SparseInt*);

	void s_call_peaks_by_macs();
	void s_receive_macs_peaks(GenomicRange);

	void s_column_slice(Eigen::ArrayX<bool> filter, bool in_place);

	void s_filter_by_parameters();
	void s_filter_by_features();
	void s_filter_by_quality_metrics();

	void s_calculate_quality_metrics();
	void s_receive_qc(QMap<QString, QList<double>>);

	void s_show_quality_matrics();

	void s_show_fragments_length_distribution();

	void s_tss_plot();
	void s_receive_tss_plot_data(Eigen::ArrayXXd tss_matrix, QVector<double> tss_vector);

	void s_coverage_plot();
	void s_receive_coverage_plot_data(COVERAGE_PLOT_ELEMENTS);

	void s_create_coverage_track();
	void s_receive_coverage_track(CoverageTrack*);

	void s_sample();

	void s_set_random_state();

	void s_set_species();

	void s_add_new_metadata();
	void s_combine_existed_metadata();
	void s_edit_metadata();

	void s_rna_scrublet();
	void s_receive_scrublet(Eigen::ArrayXd original_score, Eigen::ArrayXd simulate_score);
	void s_set_scrublet_threshold();

	void s_gsea_in_database();
	void s_gsea_from_input();
	void s_receive_gsea(GSEA gsea);

	void s_scicnv();
	void s_receive_scicnv(CNV* cnv);

	void s_infercnv();
	void s_receive_infercnv(CNV* cnv);

	void s_integrate();
	void s_receive_integrated_data(SingleCellMultiome* data, QList<const SingleCellMultiome*> items);

	void s_velocyto();
	void s_receive_velocyto(VelocytoBase* velocyto_base);

	void s_monocle3();
	void s_receive_monocle3(Monocle3*);

	void s_pando();
	void s_receive_pando(Pando pando);

	void s_show_atac_landscape();
	void s_receive_atac_landscape_plot(ATAC_LANDSCAPE_PLOT_ELEMENTS);

	void s_scent();
	void s_receive_scent(Eigen::ArrayXd);

};
