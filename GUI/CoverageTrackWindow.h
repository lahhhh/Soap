#pragma once

#include "Identifier.h"

#include "PlotsSuite.h"
#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

#include "CellMarkerTreeWidget.h"
#include "InformationTextBrowser.h"
#include "Switch.h"
#include "SignalEmitter.h"

#include "CoveragePlotWorker.h"

class CoverageTrackWindow : 
	public QMainWindow
{

	Q_OBJECT

public:
	CoverageTrackWindow(
		SingleCellMultiome* single_cell_multiome,
		CoverageTrack* coverage_track,
		SignalEmitter* signal_emitter
	);

	CoverageTrackWindow(
		SingleCellAtac* single_cell_atac,
		CoverageTrack* coverage_track,
		SignalEmitter* signal_emitter
	);

	CoverageTrackWindow(
		CoverageTrack* coverage_track,
		SignalEmitter* signal_emitter
	);

	~CoverageTrackWindow();

	static void view(
		SingleCellMultiome* single_cell_multiome,
		CoverageTrack* coverage_track,
		SignalEmitter* signal_emitter
	);

	static void view(
		SingleCellAtac* single_cell_atac,
		CoverageTrack* coverage_track,
		SignalEmitter* signal_emitter
	);

	static void view(
		CoverageTrack* coverage_track,
		SignalEmitter* signal_emitter
	);

private:	

	CoverageTrack* coverage_track_{ nullptr };

	SignalEmitter* signal_emitter_ = nullptr;

	PlotsSuite* draw_suite_ = nullptr;

	QCustomPlot* draw_area_ = nullptr;

	QWidget* main_interface_ = nullptr;

	QVBoxLayout* main_layout_ = nullptr;

	QHBoxLayout* top_layout_ = nullptr;
	QHBoxLayout* middle_layout_ = nullptr;
	QHBoxLayout* bottom_layout_ = nullptr;

	QVBoxLayout* figure_layout_ = nullptr;
	QVBoxLayout* information_layout_ = nullptr;

	InformationTextBrowser* information_area_ = nullptr;

	QPushButton* graph_setting_button_ = nullptr;
	Switch* graph_setting_switch_ = nullptr;
	QPushButton* previous_picture_button_ = nullptr;
	QPushButton* next_picture_button_ = nullptr;
	QPushButton* pop_picture_button_ = nullptr;
	QPushButton* clear_picture_button_ = nullptr;
	QPushButton* save_picture_button_ = nullptr;

	QMenu* save_picture_menu_ = nullptr;
	QAction* save_png_action_ = nullptr;
	QAction* save_jpg_action_ = nullptr;
	QAction* save_bmp_action_ = nullptr;
	QAction* save_pdf_action_ = nullptr;
	QAction* save_pdf_and_png_action_ = nullptr;

	QPushButton* left_button_ = nullptr;
	QPushButton* right_button_ = nullptr;

	QLabel* search_gene_location_label_ = nullptr;
	QLineEdit* search_gene_location_line_ = nullptr;
	QPushButton* search_gene_location_button_ = nullptr;
	QLabel* gene_location_label_ = nullptr;

	QLabel* region_label_ = nullptr;
	QLineEdit* region_line_ = nullptr;
	QPushButton* go_button_ = nullptr;

	/*
	****************************************
	*/

	COVERAGE_PLOT_ELEMENTS plot_elements_;

	Location location_;

	QMap<QCustomPlot*, Location> plot_location_;

	QList<std::tuple<QString, int, int>> available_peak_locations_;

	bool get_location_by_gene_name_ = false;

	void set_layout();

	void set_property();

	void set_top_layout();

	void set_bottom_layout();

	void set_middle_layout();

	void set_figure_layout();

	void set_information_layout();

	bool get_location();

	bool get_location_by_gene_name();

	void expand_region();

	bool find_gene_in_region();

	bool find_peak_in_region();

	bool calculate_matrix();

	void smooth_matrix();

	void draw_plot();

	void update();

public slots:

	void s_set_graph_settings();

	void s_go();

	void s_left();

	void s_right();

	void s_save_png();
	void s_save_bmp();
	void s_save_jpg();
	void s_save_pdf();
	void s_save_pdf_and_png();

	void s_previous_plot();
	void s_next_plot();
	void s_clear_plot();
	void s_pop_plot();
	void s_new_plot();

	void s_search_gene_location();

	void s_check_data(void* data, soap::VariableType type, void* item);
};

