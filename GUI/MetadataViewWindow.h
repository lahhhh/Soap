#pragma once

#include "Identifier.h"

#include "SignalEmitter.h"
#include "PlotsSuite.h"
#include "InformationTextBrowser.h"

#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"
#include "DataFrame.h"
#include "Switch.h"
#include "FeatureHandler.h"

class MetadataViewWindow 
	: public QMainWindow
{
	Q_OBJECT
public:
	MetadataViewWindow(
		SingleCellRna* single_cell_rna,
		SignalEmitter* signal_emitter
	);
	MetadataViewWindow(
		SingleCellMultiome* single_cell_multiome,
		SignalEmitter* signal_emitter
	);
	MetadataViewWindow(
		SingleCellAtac* single_cell_atac,
		SignalEmitter* signal_emitter
	);
	MetadataViewWindow(
		DataFrame* data_frame,
		SignalEmitter* signal_emitter
	);

	~MetadataViewWindow();

	void set_layout();

	void set_property();

	// ---------------------------------------//

	SignalEmitter* signal_emitter_{ nullptr };

	PlotsSuite* draw_suite_{ nullptr };

	QCustomPlot* draw_area_{ nullptr };

	QLabel* feature_label_{ nullptr };
	QLineEdit* feature_line_edit_{ nullptr };

	QLabel* normalize_label_{ nullptr };
	Switch* normalize_switch_{ nullptr };

	QLabel* gene_activity_label_{ nullptr };
	Switch* gene_activity_switch_{ nullptr };

	QLabel* main_group_label_{ nullptr };
	QComboBox* main_group_box_{ nullptr };

	QLabel* sub_group_label_{ nullptr };
	QComboBox* sub_group_box_{ nullptr };

	QPushButton* graph_setting_button_{ nullptr };
	Switch* graph_setting_switch_{ nullptr };

	QPushButton* refresh_picture_button_{ nullptr };
	QPushButton* previous_picture_button_{ nullptr };
	QPushButton* next_picture_button_{ nullptr };
	QPushButton* pop_picture_button_{ nullptr };
	QPushButton* clear_picture_button_{ nullptr };
	QPushButton* save_picture_button_{ nullptr };

	QMenu* menu_save_picture_{ nullptr };
	QAction* action_save_png_{ nullptr };
	QAction* action_save_jpg_{ nullptr };
	QAction* action_save_bmp_{ nullptr };
	QAction* action_save_pdf_{ nullptr };
	QAction* action_save_pdf_and_png_{ nullptr };

	QLabel* information_label_{ nullptr };
	InformationTextBrowser* information_area_{ nullptr };

	QHBoxLayout* main_layout_{ nullptr };
	QVBoxLayout* left_layout_{ nullptr };

	QWidget* main_interface_{ nullptr };

	QStringList valid_features_{ nullptr };

	FeatureHandler handler_;

	//------------------------------------//

	void s_check_data(void* data, soap::VariableType type, void* item);

	void plot();

	void subplot();

	void mainplot();

	void no_group_plot();

	static void view(SingleCellRna* data, SignalEmitter* signal_emitter);
	static void view(SingleCellAtac* data, SignalEmitter* signal_emitter);
	static void view(SingleCellMultiome* data, SignalEmitter* signal_emitter);
	static void view(DataFrame* data, SignalEmitter* signal_emitter);

private slots:

	void s_set_graph_settings();

	void s_check_feature_name();

	void s_save_pdf_and_png();

	void s_save_png();

	void s_save_bmp();

	void s_save_jpg();

	void s_save_pdf();

	void s_refresh_plot();

	void s_previous_plot();

	void s_next_plot();

	void s_clear_plot();

	void s_pop_plot();

	void s_new_plot();
};
