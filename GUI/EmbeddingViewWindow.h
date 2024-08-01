#pragma once

#include "Identifier.h"

#include "PlotsSuite.h"
#include "SingleCellRna.h"
#include "SingleCellMultiome.h"
#include "CellMarkerTreeWidget.h"
#include "InformationTextBrowser.h"
#include "Switch.h"
#include "SignalEmitter.h"
#include "FeatureHandler.h"

class EmbeddingViewWindow :
	public QMainWindow
{

	Q_OBJECT

public:
	EmbeddingViewWindow(
		SingleCellRna* single_cell_rna,
		Embedding* embedding,
		SignalEmitter* signal_emitter
	);

	static void view(
		SingleCellRna* single_cell_rna,
		Embedding* embedding,
		SignalEmitter* signal_emitter
	);

	EmbeddingViewWindow(
		SingleCellMultiome* single_cell_multiome,
		Embedding* embedding,
		SignalEmitter* signal_emitter
	);

	static void view(
		SingleCellMultiome* single_cell_multiome,
		Embedding* embedding,
		SignalEmitter* signal_emitter
	);

	EmbeddingViewWindow(
		SingleCellAtac* single_cell_atac,
		Embedding* embedding,
		SignalEmitter* signal_emitter
	);

	static void view(
		SingleCellAtac* single_cell_atac,
		Embedding* embedding,
		SignalEmitter* signal_emitter
	);

	~EmbeddingViewWindow();

private:

	void set_layout();

	void set_property();

	void preload();

	void set_left_layout();

	void set_middle_layout();

	void set_right_layout();

	Eigen::ArrayX<bool> get_selected();

	FeatureHandler handler_;

	SignalEmitter* signal_emitter_{ nullptr };
	Embedding* embedding_{ nullptr };

	QMap<QString, QMap<QString, QMap<QString, QStringList>>> cell_marker_database_;
	QMap<QString, QMap<QString, QStringList>> pathway_content_;

	bool plot_interacting_{ false };
	QVector<double> point_x_;
	QVector<double> point_y_;

	QMap<QCustomPlot*, QCPAxisRect*> plot_to_axis_rect_;
	QMap<QCustomPlot*, QCPGraph*> plot_to_graph_;
	QMap<QCustomPlot*, QCPCurve*> plot_active_shape_;

	QMap<QCustomPlot*, QList<QCPCurve*>> plot_shapes_;

	QList<QMetaObject::Connection> connections_;

	//--------------------------------//

	PlotsSuite* draw_suite_{ nullptr };
	QCustomPlot* draw_area_{ nullptr };

	QWidget* main_interface_{ nullptr };

	QHBoxLayout* main_layout_{ nullptr };
	QWidget* left_layout_{ nullptr };
	QWidget* middle_layout_{ nullptr };
	QVBoxLayout* right_layout_{ nullptr };

	QLabel* cell_marker_area_label_{ nullptr };
	QLabel* cell_marker_database_label_{ nullptr };
	QComboBox* cell_marker_database_box_{ nullptr };
	QLabel* cell_marker_tissue_label_{ nullptr };
	QComboBox* cell_marker_tissue_box_{ nullptr };
	CellMarkerTreeWidget* cell_marker_tree_widget_{ nullptr };

	QLabel* gene_label_{ nullptr };
	QLabel* gene_name_label_{ nullptr };
	QLineEdit* gene_line_edit_{ nullptr };
	QPushButton* gene_view_button_{ nullptr };
	QPushButton* gene_alternative_button_{ nullptr };
	QPushButton* gene_active_button_{ nullptr };

	QLabel* metadata_label_{ nullptr };
	QLabel* metadata_name_label_{ nullptr };
	QComboBox* metadata_box_{ nullptr };
	QPushButton* metadata_view_button_{ nullptr };

	QLabel* pathway_label_{ nullptr };
	QLabel* database_label_{ nullptr };
	QComboBox* database_box_{ nullptr };
	QPushButton* database_show_button_{ nullptr };
	QLabel* pathway_name_label_{ nullptr };
	QLineEdit* pathway_line_edit_{ nullptr };
	QPushButton* pathway_view_button_{ nullptr };
	QPushButton* pathway_alternative_button_{ nullptr };
	QPushButton* pathway_active_button_{ nullptr };

	QLabel* data_source_label_{ nullptr };
	QLabel* normalize_label_{ nullptr };
	Switch* normalize_switch_{ nullptr };
	QLabel* gene_activity_label_{ nullptr };
	Switch* gene_activity_switch_{ nullptr };

	QLabel* alternative_label_{ nullptr };
	CellMarkerTreeWidget* alternative_tree_widget_{ nullptr };

	QLabel* active_label_{ nullptr };
	QPushButton* active_clear_button_{ nullptr };
	QPushButton* active_refresh_button_{ nullptr };
	CellMarkerTreeWidget* active_tree_widget_{ nullptr };

	QPushButton* select_cell_button_{ nullptr };
	QPushButton* graph_setting_button_{ nullptr };
	Switch* graph_setting_switch_{ nullptr };
	QPushButton* previous_picture_button_{ nullptr };
	QPushButton* next_picture_button_{ nullptr };
	QPushButton* pop_picture_button_{ nullptr };
	QPushButton* clear_picture_button_{ nullptr };
	QPushButton* save_picture_button_{ nullptr };

	QMenu* save_picture_menu_{ nullptr };
	QAction* save_png_action_{ nullptr };
	QAction* save_jpg_action_{ nullptr };
	QAction* save_bmp_action_{ nullptr };
	QAction* save_pdf_action_{ nullptr };
	QAction* save_pdf_and_png_action_{ nullptr };

	InformationTextBrowser* information_area_{ nullptr };

	//----------------------------------------------//

	void show_gene(const QString& gene_name);

	void show_type(const QString& type_name, const QStringList& markers);

	void clear_interactive();

	void final_connect();

signals:

	void x_gene_to_alternative(QString gene_name);

	void x_gene_to_active(QString gene_name);

	void x_type_to_alternative(QString typeName, QStringList markers);

	void x_type_to_active(QString typeName, QStringList markers);

private slots:

	void s_check_data(void* data, soap::VariableType type, void* item);

	void s_set_graph_settings();

	void s_refresh_active_items();

	void s_clear_active_items();

	void s_select_cell();

	void s_plot_mouse_press(QMouseEvent*);

	void s_plot_mouse_double_click(QMouseEvent*);

	void s_plot_mouse_move(QMouseEvent*);

	void s_pathway_to_alternative();

	void s_pathway_to_active();

	void s_explore_pathway();

	void s_view_metadata();

	void s_view_pathway();

	void s_view_gene();
	void s_alternative_gene();
	void s_active_gene();

	void s_update_tissue();

	void s_update_cell_marker();

	void s_show_gene(const QString& gene_name);
	void s_show_type(const QString& type_name, const QStringList& markers);

	void s_save_png();
	void s_save_bmp();
	void s_save_jpg();
	void s_save_pdf();
	void s_save_pdf_and_png();

	void s_refresh_plot();
	void s_previous_plot();
	void s_next_plot();
	void s_clear_plot();
	void s_pop_plot();
	void s_new_plot();
};
