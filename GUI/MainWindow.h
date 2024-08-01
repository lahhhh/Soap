#pragma once

#include "Identifier.h"

#include <QMainWindow>
#include <QTreeWidgetItem>
#include <QCloseEvent>

#include "BulkRnaItem.h"
#include "SingleCellRnaItem.h"
#include "SingleCellAtacItem.h"
#include "SingleCellMultiomeItem.h"
#include "DataFrameItem.h"
#include "GeneNameItem.h"
#include "StringVectorItem.h"
#include "NumericMatrixItem.h"

#include "VariableTreeWidget.h"
#include "qCustomPlot.h"
#include "Switch.h"

class MainWindow : 
	public QMainWindow
{
	Q_OBJECT

public:

	MainWindow(QWidget* parent = nullptr);
	~MainWindow();

	SignalEmitter* signal_emitter_{ nullptr };

	QWidget* main_interface_{ nullptr };

	QMenuBar* menubar_{ nullptr };

	/*
	* **************************************
	*/

	QHBoxLayout* main_layout_{ nullptr };
	QHBoxLayout* right_top_layout_{ nullptr };

	QVBoxLayout* left_layout_{ nullptr };
	QVBoxLayout* right_layout_{ nullptr };

	/*
* **************************************
*/

	QWidget* left_panel_{ nullptr };
	QLabel* variable_label_{ nullptr };
	VariableTreeWidget* variable_tree_widget_{ nullptr };
	QLabel* message_label_{ nullptr };
	QPushButton* log_export_button_{ nullptr };
	InformationTextBrowser* information_area_{ nullptr };

	/*
* **************************************
*/

	QPushButton* graph_setting_button_{ nullptr };
	Switch* graph_setting_switch_{ nullptr };

	QPushButton* previous_picture_button_{ nullptr };
	QPushButton* next_picture_button_{ nullptr };
	QPushButton* clear_picture_button_{ nullptr };
	QPushButton* pop_picture_button_{ nullptr };
	QPushButton* save_picture_button_{ nullptr };
	QMenu* menu_save_picture_{ nullptr };
	QAction* action_save_png_{ nullptr };
	QAction* action_save_jpg_{ nullptr };
	QAction* action_save_bmp_{ nullptr };
	QAction* action_save_pdf_{ nullptr };
	QAction* action_save_pdf_and_png_{ nullptr };

	PlotsSuite* draw_suite_{ nullptr };

	QCustomPlot* draw_area_{ nullptr };

	template<typename VariableType>
	void set_normal_item(void* ptr, const QString& name) {

		VariableType* data = static_cast<VariableType*>(ptr);
		QString title = this->signal_emitter_->get_unique_name(name);

		auto VariableIdentifier = soap::type<VariableType>();
		this->signal_emitter_->update_information(title, VariableIdentifier, data, true);

		auto* item = new soap::item_type<VariableType>(
			title, 
			this->signal_emitter_->new_index_tree(), 
			data, 
			this->draw_suite_,
			this->information_area_,
			this->signal_emitter_
			);
		
		this->variable_tree_widget_->addTopLevelItem(item);
		
		G_LOG(title + " has been constructed");
	}

private slots:

	void __s_work_finished();

	void s_multiple_group_gsea_plot();

	void s_multiple_group_enrichment_plot_select_pathway();

	void s_multiple_group_enrichment_plot_top_n();

	void s_faq();

	void s_scRNAseq_pipeline();

	void closeEvent(QCloseEvent* e);

	void s_export_log();

	void s_set_graph_setting();

	void s_load_10X_scMultiome();

	void s_load_10X_scRNA();

	void s_load_fragments_scATAC();

	void s_load_scRNAseq_count_matrix();
	void s_batch_load_count_matrix();

	void s_test();

	void s_save_pdf_and_png();

	void s_save_png();

	void s_save_bmp();

	void s_save_jpg();

	void s_save_pdf();

	void s_previous_plot();

	void s_next_plot();

	void s_clear_plot();

	void s_pop_plot();

	void s_new_plot();

	void s_read_data_frame_default();
	void s_read_data_frame_fast();

	void s_load_item();

	void s_receive_data_frame(CustomMatrix* mat);

	void __s_receive_message(const QString& message, int mode);

	void s_create_data(void* data, soap::VariableType type, QString name);

	void s_create_string_vector_from_input();

	void s_generate_wix();

private:

	void prepare();

	void set_signal_emitter();

	void set_variable_tree_widget();

	void set_visualize_menu();

	void set_file_menu();

	void set_guide_menu();

	void set_utility_menu();

	void set_left_layout();

	void set_right_layout();

	void set_main_layout();

	void set_plot_suite();

	void set_property();
};
