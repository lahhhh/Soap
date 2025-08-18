#include <WinSock2.h>
#include <boost/process.hpp>

#include "MainWindow.h"
#include <iostream>

#include <QFileDialog>
#include <QTextEdit>

#include "Read10xRnaWorker.h"
#include "BatchLoading10XRnaWorker.h"
#include "Read10xRnaAndAtacWorker.h"
#include "PlotWindow.h"
#include "TableReadingWorker.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"
#include "CountMatrixReadingWorker.h"
#include "CountMatrixBatchLoadingWorker.h"
#include "LoadAtacFromFragmentsWorker.h"

#include "GraphSettingDialog.h"
#include "YesOrNoDialog.h"
#include "PageWindow.h"
#include "EnrichmentMultiGroupPlotDialog.h"
#include "CustomPlot.h"
#include "Poisson.h"

#include "TwoBitFileProcessor.h"

#include "SoapGUI.h"


MainWindow::MainWindow(QWidget* parent)
	: QMainWindow(parent)
{

	set_signal_emitter();

	this->menubar_ = new QMenuBar(this);
	this->setMenuBar(this->menubar_);

	set_file_menu();
	set_visualize_menu();
	set_guide_menu();
	set_utility_menu();
	set_left_layout();
	set_right_layout();
	set_main_layout();
	set_plot_suite();

	prepare();

	this->main_interface_ = new QWidget();
	this->main_interface_->setLayout(this->main_layout_);
	this->setCentralWidget(this->main_interface_);

	set_property();
}

void MainWindow::set_property() {
	this->setGeometry(100, 100, 1300, 800);
	this->setWindowTitle("SOAP");
	G_SET_ICON;
};

void MainWindow::prepare() {
	connect(this->signal_emitter_, &SignalEmitter::x_log, this->information_area_, &InformationTextBrowser::s_receive_log);
	connect(this->signal_emitter_, &SignalEmitter::x_notice, this->information_area_, &InformationTextBrowser::s_receive_notice);
	connect(this->signal_emitter_, &SignalEmitter::x_warn, this->information_area_, &InformationTextBrowser::s_receive_warning);
};

MainWindow::~MainWindow()
{
	delete this->draw_suite_;
}

void MainWindow::set_plot_suite() {
	this->draw_suite_ = new PlotsSuite();

	connect(this->graph_setting_switch_, &Switch::toggled, this->draw_suite_, &PlotsSuite::s_setting_activate);
	connect(this->draw_suite_, &PlotsSuite::x_plot_prepared, this, &MainWindow::s_new_plot);

	this->draw_suite_->prepare();
};

void MainWindow::set_main_layout() {
	this->main_layout_ = new QHBoxLayout();

	this->main_layout_->addWidget(this->left_panel_);
	this->main_layout_->addLayout(this->right_layout_);
};

void MainWindow::set_right_layout() {
	this->right_layout_ = new QVBoxLayout();

	this->graph_setting_button_ = new QPushButton("Graph Settings");
	this->graph_setting_button_->setFixedSize(150, 30);

	connect(this->graph_setting_button_, &QPushButton::clicked, this, &MainWindow::s_set_graph_setting);

	this->graph_setting_switch_ = new Switch(false, this->graph_setting_button_);
	this->graph_setting_switch_->setFixedSize(120, 30);


	G_SET_PLOTSUITE_BUTTON;
	this->menu_save_picture_ = new QMenu();

	this->action_save_png_ = new QAction("PNG");
	this->action_save_jpg_ = new QAction("JPG");
	this->action_save_bmp_ = new QAction("BMP");
	this->action_save_pdf_ = new QAction("PDF");
	this->action_save_pdf_and_png_ = new QAction("PDF and PNG");

	this->menu_save_picture_->addAction(this->action_save_png_);
	this->menu_save_picture_->addAction(this->action_save_jpg_);
	this->menu_save_picture_->addAction(this->action_save_bmp_);
	this->menu_save_picture_->addAction(this->action_save_pdf_);
	this->menu_save_picture_->addAction(this->action_save_pdf_and_png_);

	connect(this->action_save_jpg_, &QAction::triggered, this, &MainWindow::s_save_jpg);
	connect(this->action_save_png_, &QAction::triggered, this, &MainWindow::s_save_png);
	connect(this->action_save_bmp_, &QAction::triggered, this, &MainWindow::s_save_bmp);
	connect(this->action_save_pdf_, &QAction::triggered, this, &MainWindow::s_save_pdf);
	connect(this->action_save_pdf_and_png_, &QAction::triggered, this, &MainWindow::s_save_pdf_and_png);

	connect(this->previous_picture_button_, &QPushButton::clicked, this, &MainWindow::s_previous_plot);
	connect(this->next_picture_button_, &QPushButton::clicked, this, &MainWindow::s_next_plot);
	connect(this->clear_picture_button_, &QPushButton::clicked, this, &MainWindow::s_clear_plot);
	connect(this->pop_picture_button_, &QPushButton::clicked, this, &MainWindow::s_pop_plot);

	this->save_picture_button_->setMenu(this->menu_save_picture_);

	this->right_top_layout_ = new QHBoxLayout();

	this->right_top_layout_->addWidget(this->graph_setting_button_);
	this->right_top_layout_->addWidget(this->graph_setting_switch_);

	this->right_top_layout_->addStretch();

	this->right_top_layout_->addWidget(this->previous_picture_button_);
	this->right_top_layout_->addWidget(this->next_picture_button_);
	this->right_top_layout_->addWidget(this->clear_picture_button_);
	this->right_top_layout_->addWidget(this->pop_picture_button_);
	this->right_top_layout_->addWidget(this->save_picture_button_);

	this->right_layout_->addLayout(this->right_top_layout_);
};

void MainWindow::set_left_layout() {
	this->left_panel_ = new QWidget(this);

	this->left_panel_->setFixedWidth(500);
	this->left_panel_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);

	this->left_layout_ = new QVBoxLayout();

	this->left_layout_->setContentsMargins(0, 0, 0, 0);
	this->left_panel_->setLayout(this->left_layout_);

	G_SET_LABEL(this->variable_label_, "Variable", soap::LargeSize);
	G_SET_LARGE_LABEL(this->variable_label_);

	G_SET_LABEL(this->message_label_, "Message", soap::LargeSize);
	G_SET_LARGE_LABEL(this->message_label_);

	G_SET_BUTTON(this->log_export_button_, "Export", soap::MiddleSize);
	G_SET_BUTTON(this->log_clear_button_, "Clear", soap::MiddleSize);
	connect(this->log_export_button_, &QPushButton::clicked, this, &MainWindow::s_export_log);
	connect(this->log_clear_button_, &QPushButton::clicked, this, &MainWindow::s_clear_log);

	QHBoxLayout* row_layout = new QHBoxLayout;
	row_layout->addWidget(this->message_label_);
	row_layout->addStretch();
	row_layout->addWidget(this->log_export_button_);
	row_layout->addWidget(this->log_clear_button_);

	this->information_area_ = new InformationTextBrowser(this);

	set_variable_tree_widget();

	this->left_layout_->addWidget(this->variable_label_);
	this->left_layout_->addWidget(this->variable_tree_widget_);
	this->left_layout_->addLayout(row_layout);
	this->left_layout_->addWidget(this->information_area_);

	this->left_layout_->setStretchFactor(this->variable_tree_widget_, 5);
	this->left_layout_->setStretchFactor(this->information_area_, 5);
};

void MainWindow::set_utility_menu() {

	auto menu_utility = this->menubar_->addMenu("Utility");

	auto menu_create_variable = menu_utility->addMenu("Create Variable");

	auto menu_create_string_vector = menu_create_variable->addMenu("String Vector");

	menu_create_string_vector->addAction("From Input", this, &MainWindow::s_create_string_vector_from_input);
};

void MainWindow::s_create_string_vector_from_input() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Create String Vector",
		{ "Name", "Content" },
		{soap::InputStyle::StringLineEdit, soap::InputStyle::TextEdit}
	);

	if (settings.isEmpty()) {
		return;
	}

	QString var_name = settings[0];

	if (var_name.isEmpty()) {
		G_WARN("Variable Name can not be empty.");
		return;
	}

	StringVector* data = new StringVector(custom::split_lines(settings[1]));
	this->set_normal_item<StringVector>(data, var_name);
};

void MainWindow::set_guide_menu() {

	auto menu_guide = this->menubar_->addMenu("Guide");

	auto menu_pipeline = menu_guide->addMenu("Pipeline");

	menu_guide->addAction("FAQ", this, &MainWindow::s_faq);

	menu_pipeline->addAction("Single Cell RNA-seq Pipeline", this, &MainWindow::s_scRNAseq_pipeline);
};

void MainWindow::set_visualize_menu() {

	auto menu_visualize = this->menubar_->addMenu("Visualize");

	auto menu_multiple_group_plot = menu_visualize->addMenu("Multi Group Plot");

	auto menu_multiple_group_plot_enrichment = menu_multiple_group_plot->addMenu("Enrichment Barplot");

	menu_multiple_group_plot_enrichment->addAction("Select Pathway", this, &MainWindow::s_multiple_group_enrichment_plot_select_pathway);
	
	menu_multiple_group_plot_enrichment->addAction("Top N", this, &MainWindow::s_multiple_group_enrichment_plot_top_n);

	menu_multiple_group_plot->addAction("GSEA Heatmap", this, &MainWindow::s_multiple_group_gsea_plot);
};

void MainWindow::set_file_menu() {

	auto menu = this->menubar_->addMenu("File");

	auto read_dataframe_menu = menu->addMenu("Read Data Frame");

	read_dataframe_menu->addAction("Default", this, &MainWindow::s_read_data_frame_default);
	read_dataframe_menu->addAction("Fast(No Unicode)", this, &MainWindow::s_read_data_frame_fast);

	auto scrnaseq_menu = menu->addMenu("sc-RNA seq");

	auto rna_10X_menu = scrnaseq_menu->addMenu("Load 10X");

	rna_10X_menu->addAction("Default", this, &MainWindow::s_load_10X_scRNA);
	rna_10X_menu->addAction("Default(folder)", this, &MainWindow::s_load_10X_scRNA_folder);
	rna_10X_menu->addAction("Batch", this, &MainWindow::s_batch_load_10X_scRNA);

	auto count_menu = scrnaseq_menu->addMenu("Load Count Matrix");

	count_menu->addAction("Default", this, &MainWindow::s_load_scRNAseq_count_matrix);
	count_menu->addAction("Batch", this, &MainWindow::s_batch_load_count_matrix);

	auto scatac_menu = menu->addMenu("sc-ATAC seq");

	scatac_menu->addAction("Load From Fragments", this, &MainWindow::s_load_fragments_scATAC);

	auto scmultiome_menu = menu->addMenu("sc-Multiome seq");

	scmultiome_menu->addAction("Load 10X scMultiome", this, &MainWindow::s_load_10X_scMultiome);

	menu->addAction("Load SOAP Item", this, &MainWindow::s_load_item);
};

void MainWindow::set_signal_emitter() {
	this->signal_emitter_ = new SignalEmitter(this);
	connect(this->signal_emitter_, &SignalEmitter::x_data_create_soon, this, &MainWindow::s_create_data);
};

void MainWindow::s_create_data(void* data, soap::VariableType type, QString name) {
	if (data == nullptr) {
		return;
	}
	switch (type)
	{
	case soap::VariableType::AnyVariable:
		break;
	case soap::VariableType::SingleCellRna:
		set_normal_item<SingleCellRna>(data, name);
		break;
	case soap::VariableType::SingleCellAtac:
		set_normal_item<SingleCellAtac>(data, name);
		break;
	case soap::VariableType::DifferentialAnalysis:
		set_normal_item<DifferentialAnalysis>(data, name);
		break;
	case soap::VariableType::DenseInt:
		set_normal_item<DenseInt>(data, name);
		break; 
	case soap::VariableType::BulkRna:
		set_normal_item<BulkRna>(data, name);
		break;
	case soap::VariableType::DenseDouble:
		set_normal_item<DenseDouble>(data, name);
		break;
	case soap::VariableType::SparseDouble:
		set_normal_item<SparseDouble>(data, name);
		break;
	case soap::VariableType::SparseInt:
		set_normal_item<SparseInt>(data, name);
		break;
	case soap::VariableType::Metadata:
		set_normal_item<Metadata>(data, name);
		break;
	case soap::VariableType::Embedding:
		set_normal_item<Embedding>(data, name);
		break;
	case soap::VariableType::Enrichment:
		set_normal_item<Enrichment>(data, name);
		break;
	case soap::VariableType::DataFrame:
		set_normal_item<DataFrame>(data, name);
		break;
	case soap::VariableType::GSEA:
		set_normal_item<GSEA>(data, name);
		break;
	case soap::VariableType::CellChat:
		set_normal_item<CellChat>(data, name);
		break;
	case soap::VariableType::CNV:
		set_normal_item<CNV>(data, name);
		break;
	case soap::VariableType::SingleCellMultiome:
		set_normal_item<SingleCellMultiome>(data, name);
		break;
	case soap::VariableType::Pando:
		set_normal_item<Pando>(data, name);
		break;
	case soap::VariableType::GenomicRange:
		set_normal_item<GenomicRange>(data, name);
		break;
	case soap::VariableType::MotifPosition:
		set_normal_item<MotifPosition>(data, name);
		break;
	case soap::VariableType::CoverageTrack:
		set_normal_item<CoverageTrack>(data, name);
		break;
	case soap::VariableType::Footprint:
		set_normal_item<Footprint>(data, name);
		break;
	case soap::VariableType::StringVector:
		set_normal_item<StringVector>(data, name);
		break;
	case soap::VariableType::GeneName:
		set_normal_item<GeneName>(data, name);
		break;
	case soap::VariableType::NumericMatrix:
		set_normal_item<NumericMatrix>(data, name);
		break;
	case soap::VariableType::ChromVAR:
		set_normal_item<ChromVAR>(data, name);
		break;
	case soap::VariableType::Cicero:
		set_normal_item<Cicero>(data, name);
		break;
	default:
		break;
	}
};

void MainWindow::s_faq() {
	PageWindow::show_this("FAQ", FILE_FAQ);
};

void MainWindow::set_variable_tree_widget() {

	this->variable_tree_widget_ = new VariableTreeWidget(this, this->draw_suite_, this->information_area_, this->signal_emitter_);
	this->variable_tree_widget_->setFixedWidth(500);
	this->variable_tree_widget_->setColumnCount(3);
	this->variable_tree_widget_->setHeaderLabels(QStringList{ "Object" , "Type" , "Size" });
	this->variable_tree_widget_->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
	this->variable_tree_widget_->header()->setStretchLastSection(false);
	this->variable_tree_widget_->setTextElideMode(Qt::ElideNone);
	this->variable_tree_widget_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
};

void MainWindow::closeEvent(QCloseEvent* e) {
	bool close = YesOrNoDialog::get_response("Close soap", "Sure to Quit?");
	if (!close) {
		e->ignore();
	}
};

void MainWindow::s_batch_load_count_matrix() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Batch Loading Settings",
		{ "Choose Files", "Delimiter" },
		{ soap::InputStyle::MultiFile, soap::InputStyle::ComboBox},
		{ {}, { "auto-detect", "tab(\'\\t\')", "space(\' \')", "comma(\',\')" }}
	);

	if (settings.isEmpty())return;

	QString delimiter;
	if (settings[1] == "auto-detect") {
		delimiter = "auto-detect";
	}
	else {
		delimiter = settings[1];
		if (delimiter == "tab(\'\\t\')") {
			delimiter = "\t";
		}
		else if (delimiter == "space(\' \')") {
			delimiter = " ";
		}
		else {
			delimiter = ",";
		}
	}

	QStringList files = multiple_file_to_list(settings[0]);

	if (files.size() < 2) {
		G_WARN("Invalid File Input Number.");
		return;
	}
	CountMatrixBatchLoadingWorker* worker = new CountMatrixBatchLoadingWorker(
		files,
		delimiter
	);
	G_LINK_WORKER_THREAD(CountMatrixBatchLoadingWorker, x_data_create_soon, MainWindow, s_create_data);
};

void MainWindow::s_load_scRNAseq_count_matrix() {

	QString file_path = QFileDialog::getOpenFileName(this, tr("Choose Count Matrix File"));
	if (file_path.isEmpty())return;

	QStringList setting = CommonDialog::get_response(
		this->signal_emitter_, "Delimiter Setting",
		{ "Delimiter", "Custom" },
		{ soap::InputStyle::ComboBox, soap::InputStyle::StringLineEdit},
		{ { "auto-detect", "tab(\'\\t\')", "space(\' \')", "comma(\',\')", "custom"}}
	);

	if (setting.isEmpty())return;

	QString delimiter;
	if (setting[0] == "custom") {
		delimiter = setting[1];
	}
	else if (setting[0] == "auto-detect") {
		delimiter = "auto-detect";
	}
	else {
		delimiter = setting[0];
		if (delimiter == "tab(\'\\t\')") {
			delimiter = "\t";
		}
		else if (delimiter == "space(\' \')") {
			delimiter = " ";
		}
		else {
			delimiter = ",";
		}
	}
	G_LOG("Start loading count matrix from " + file_path);

	CountMatrixReadingWorker* worker = new CountMatrixReadingWorker(file_path, delimiter);
	G_LINK_WORKER_THREAD(CountMatrixReadingWorker, x_data_create_soon, MainWindow, s_create_data);
};

void MainWindow::s_set_graph_setting() {
	GraphSettingDialog::set_graph_setting(this->draw_suite_->graph_settings_);
	this->graph_setting_switch_->set_status(this->draw_suite_->graph_settings_.active());
};

void MainWindow::s_scRNAseq_pipeline() {
	PageWindow::show_this("Single-cell RNA-seq Pipeline", FILE_SINGLE_CELL_RNA_PIPELINE);
};

void MainWindow::s_load_10X_scMultiome() {

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose Files",
		{ "Barcodes File", "Feature File", "Matrix File" },
		{ soap::InputStyle::ChooseOpenFile, soap::InputStyle::ChooseOpenFile, soap::InputStyle::ChooseOpenFile }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString barcodes_file_name = settings[0];
	QString features_file_name = settings[1];
	QString matrix_file_name = settings[2];

	if (barcodes_file_name.isEmpty() || features_file_name.isEmpty() || matrix_file_name.isEmpty()) {
		return;
	}

	Read10XMultiomeWorker* worker = new Read10XMultiomeWorker(barcodes_file_name, features_file_name, matrix_file_name);
	G_LINK_WORKER_THREAD(Read10XMultiomeWorker, x_data_create_soon, MainWindow, s_create_data);

};

void MainWindow::s_load_fragments_scATAC() {

	QString file_name = QFileDialog::getOpenFileName(this, "Choose Fragments File", "", "Fragments File(*.tsv.gz)");

	if (file_name.isEmpty()) {
		return;
	}

	LoadAtacFromFragmentsWorker* worker = new LoadAtacFromFragmentsWorker(file_name);
	G_LINK_WORKER_THREAD(LoadAtacFromFragmentsWorker, x_data_create_soon, MainWindow, s_create_data);
}

void MainWindow::s_batch_load_10X_scRNA() {

	auto dir = QFileDialog::getExistingDirectory(nullptr, "Choose Root Directory");

	if (dir.isEmpty()) {
		return;
	}

	auto worker = new BatchLoading10XRnaWorker(dir);

	G_LINK_WORKER_THREAD(BatchLoading10XRnaWorker, x_data_create_soon, MainWindow, s_create_data);
};

void MainWindow::s_load_10X_scRNA_folder() {

	QString dir_name = QFileDialog::getExistingDirectory(this, tr("Choose Directory"));
	if (!dir_name.isEmpty()) {

		Read10XRnaWorker* worker = new Read10XRnaWorker(dir_name);
		G_LINK_WORKER_THREAD(Read10XRnaWorker, x_data_create_soon, MainWindow, s_create_data);
	}
};

void MainWindow::s_load_10X_scRNA()
{

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Choose Files",
		{ "Barcodes File", "Feature File", "Matrix File" },
		{ soap::InputStyle::ChooseOpenFile, soap::InputStyle::ChooseOpenFile, soap::InputStyle::ChooseOpenFile }
	);

	if (settings.isEmpty()) {
		return;
	}

	QString barcodes_file_name = settings[0];
	QString features_file_name = settings[1];
	QString matrix_file_name = settings[2];

	if (barcodes_file_name.isEmpty() || features_file_name.isEmpty() || matrix_file_name.isEmpty()) {
		return;
	}

	Read10XRnaWorker* worker = new Read10XRnaWorker(barcodes_file_name, features_file_name, matrix_file_name);
	G_LINK_WORKER_THREAD(Read10XRnaWorker, x_data_create_soon, MainWindow, s_create_data);
}

void MainWindow::s_load_item() {
	QString file_path = QFileDialog::getOpenFileName(this, "Select Item File", "", "soap item(*.sif)");
	if (file_path.isEmpty())return;

	ItemIOWorker* worker = new ItemIOWorker(file_path);
	G_LINK_WORKER_THREAD(ItemIOWorker, x_data_create_soon, MainWindow, s_create_data)
};

void MainWindow::s_clear_log() {

	this->information_area_->clear();
};

void MainWindow::s_export_log() {
	QString file_name = QFileDialog::getSaveFileName(this, "Set Log Name", "", "LOG(*.log)");
	if (file_name.isEmpty())return;

	QFile file(file_name);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		G_WARN("Can not create file " + file_name + " .");
		return;
	}

	QTextStream out(&file);
	out << this->information_area_->toPlainText();
	G_LOG("Log Saved.");
};

void MainWindow::s_save_png() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);
	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "PNG(*.png)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->savePng(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MainWindow::s_save_bmp() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);
	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "BMP(*.bmp)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->saveBmp(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MainWindow::s_save_jpg() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);
	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set Picture Name", "", "JPG(*.jpg)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->saveJpg(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MainWindow::s_save_pdf_and_png() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"Picture Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);
	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set PNG Name", "", "PNG(*.PNG)");
	if (picture_name.isEmpty())return;

	QString pdf_name = QFileDialog::getSaveFileName(this, "Set PDF Name", "", "PDF(*.pdf)");
	if (pdf_name.isEmpty())return;

	bool success = this->draw_area_->savePdf(pdf_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (!success) {
		G_LOG("PDF saving failed")
			return;
	}
	else {
		G_LOG("PDF saved.");
	}

	std::string pdf_file_name = pdf_name.toStdString();
	std::string picture_file_name = picture_name.toStdString();

	if (!custom::save_pdf_page_as_png(pdf_file_name, 0, picture_file_name)) {
		G_WARN("PNG Saving Failed.");
	}
	else {
		G_LOG("PNG saved.");
	}
}

void MainWindow::s_save_pdf() {
	QStringList picture_size = CommonDialog::get_response(
		this->signal_emitter_,
		"PDF Settings",
		{ "Height (pixel):" + QString::number(this->draw_area_->height()), "Width (pixel):" + QString::number(this->draw_area_->width()) },
		{ soap::InputStyle::IntegerLineEdit, soap::InputStyle::IntegerLineEdit}
	);
	if (picture_size.isEmpty())return;

	QString picture_name = QFileDialog::getSaveFileName(this, "Set PDF Name", "", "PDF(*.pdf)");
	if (picture_name.isEmpty())return;

	bool success = this->draw_area_->savePdf(picture_name, picture_size[1].toInt(), picture_size[0].toInt());
	if (success) {
		G_LOG("Picture saved.");
	}
	else {
		G_LOG("Picture saving failed");
	}
};

void MainWindow::s_previous_plot() {

	if (this->draw_suite_->current_plot_id_ <= 1) {
		return;
	}

	this->right_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[-- this->draw_suite_->current_plot_id_];
	this->right_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
	
};

void MainWindow::s_next_plot() {
	if (this->draw_suite_->current_plot_id_ == this->draw_suite_->maximum_plot_id_) {
		return;
	}

	this->right_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);

	this->draw_area_ = this->draw_suite_->plots_[++ this->draw_suite_->current_plot_id_];
	this->right_layout_->addWidget(this->draw_area_);
	this->draw_area_->setVisible(true);
	this->draw_area_->replot();
	
};

void MainWindow::s_clear_plot() {
	this->right_layout_->removeWidget(this->draw_area_);
	this->draw_area_->setVisible(false);
	this->draw_area_ = nullptr;
	this->draw_suite_->clear();
};

void MainWindow::s_pop_plot() {

	PlotWindow::show_plot(this->draw_suite_, "Figure", this->draw_area_->width(), this->draw_area_->height());
};

void MainWindow::s_new_plot() {

	if (this->draw_area_ == nullptr) {
		this->draw_area_ = this->draw_suite_->current_plot_; 
		this->right_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();	
	}
	else {

		this->right_layout_->removeWidget(this->draw_area_);
		this->draw_area_->setVisible(false);

		this->draw_area_ = this->draw_suite_->current_plot_;
		this->right_layout_->addWidget(this->draw_area_);
		this->draw_area_->setVisible(true);
		this->draw_area_->replot();	
	}
};

void MainWindow::s_read_data_frame_default() {

	QString file_name = QFileDialog::getOpenFileName();
	if (file_name.isEmpty()) {
		return;
	}

	G_LOG("Begin loading " + file_name);

	TableReadingWorker* worker = new TableReadingWorker(file_name, false);
	G_LINK_WORKER_THREAD(TableReadingWorker, x_data_frame_ready, MainWindow, s_receive_data_frame);
};

void MainWindow::s_read_data_frame_fast() {

	QString file_name = QFileDialog::getOpenFileName();
	if (file_name.isEmpty()) {
		return;
	}

	G_LOG("Begin loading " + file_name);

	TableReadingWorker* worker = new TableReadingWorker(file_name, true);
	G_LINK_WORKER_THREAD(TableReadingWorker, x_data_frame_ready, MainWindow, s_receive_data_frame);
};

void MainWindow::s_receive_data_frame(CustomMatrix* mat) {

	DataFrame* dataFrame = new DataFrame(*mat);
	delete mat;

	this->set_normal_item<DataFrame>(dataFrame, "DataFrame");
	G_LOG("DataFrame loading finished.");
};

void MainWindow::__s_receive_message(const QString& message, int mode) {
	if (mode == 0) {
		G_LOG(message);
	}
	else if (mode == 1) {
		G_WARN(message);
	}
	else {
		G_NOTICE(message);
	}
};

void MainWindow::__s_work_finished() {
};

void MainWindow::s_multiple_group_gsea_plot() {
	auto gseas = this->signal_emitter_->get_type_variable(soap::VariableType::GSEA);
	if (gseas.size() < 2) {
		this->signal_emitter_->x_warn("No enough gsea found.");
		return;
	}

	QStringList gsea_name, gsea_paths;
	for (auto& gsea : gseas.keys()) {
		auto ptr = (GSEA*)gseas[gsea];
		gsea_name << gsea;
		gsea_paths << ptr->mat_.rownames_;
	}
	gsea_paths = custom::unique(gsea_paths);

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"GSEA Plot Settings",
		{ "Group:Custom Label", "Pathway:Custom Label", "NES Filter", "NES Threshold (<):0.05" },
		{ soap::InputStyle::SimpleChoiceWithLineEditAndColorChoiceLayout,
		soap::InputStyle::MultipleDoubleLineEditWithCompleterLayout
		, soap::InputStyle::ComboBox, soap::InputStyle::NumericLineEdit},
		{ gsea_name, gsea_paths, { "No Filter", "P", "FDR", "FWER" }}
	);
	if (settings.isEmpty())return;

	auto [gsea_names, gsea_alias, gsea_colors] = simple_choice_with_line_edit_and_color_choice_layout_to_tuple(settings[0]);
	if (gsea_names.isEmpty())return;

	auto [path_names, path_alias] = multiple_double_line_edit_with_completer_layout_to_pair(settings[1]);
	if (path_names.isEmpty())return;

	QString filter_type = settings[2];
	double threshold = settings[3].toDouble();
	if (threshold <= 0) {
		G_WARN("Illegal threshold value : " + settings[3]);
		return;
	}
	int gsea_number = gsea_names.size(), path_number = path_names.size();

	QStringList gsea_labels = gsea_names, path_labels = path_names;
	for (int i = 0; i < gsea_number; ++i) {
		if (!gsea_alias[i].isEmpty()) {
			gsea_labels[i] = gsea_alias[i];
		}
	}
	for (int i = 0; i < path_number; ++i) {
		if (!path_alias[i].isEmpty()) {
			path_labels[i] = path_alias[i];
		}
	}
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(path_number, gsea_number);

	for (int i = 0; i < gsea_number; ++i) {
		GSEA* ptr = (GSEA*)gseas[gsea_names[i]];
		auto nes = ptr->get_nes(path_names, filter_type, threshold);
		for (int j = 0; j < path_number; ++j) {
			data(j, i) = nes[j];
		}
	}
	auto& gs = this->draw_suite_->graph_settings_;
	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);
	custom_plot::patch::set_range(axis_rect, QCPRange(0, 20 * gsea_number), QCPRange(0, 20 * (path_number + 1) - 5));
	QColor col;
	for (int i = 0; i < path_number; ++i) {
		for (int j = 0; j < gsea_number; ++j) {
			double val = data(i, j);
			if (val > 0) {
				val = std::min(1.0, val / 3.);
				col = custom_plot::utility::gradient_color(val, gs.get_gradient_middle_color(), gs.get_gradient_high_color());
			}
			else {
				val = std::max(-1.0, val / 3.);
				col = custom_plot::utility::gradient_color(val + 1, gs.get_gradient_low_color(), gs.get_gradient_middle_color());
			}
			custom_plot::patch::shape_borderless(draw_area, axis_rect, QVector<double>() << j * 20 + 1 << j * 20 + 1 << j * 20 + 19 << j * 20 + 19,
				QVector<double>() << (path_number - i - 1) * 20 + 1 << (path_number - i - 1) * 20 + 19 << (path_number - i - 1) * 20 + 19 << (path_number - i - 1) * 20 + 1,
				col);
		}
	}
	for (int i = 0; i < gsea_number; ++i) {
		custom_plot::patch::shape_borderless(draw_area, axis_rect, QVector<double>() << i * 20 + 1 << i * 20 + 1 << i * 20 + 19 << i * 20 + 19,
			QVector<double>() << path_number * 20 + 5 << path_number * 20 + 10 << path_number * 20 + 10 << path_number * 20 + 5, gsea_colors[i]);
	}
	custom_plot::patch::remove_left_bottom_axis(axis_rect);
	custom_plot::set_left_axis_label(axis_rect, Eigen::ArrayXd::LinSpaced(path_number, path_number * 20 - 10, 10), path_labels, 0, gs);
	custom_plot::add_square_legend(draw_area, legend_layout, gsea_labels, gsea_colors, gs.get_legend_label_font(), gs.get_legend_title("Group"), gs.get_legend_title_font(), gs.get_legend_column_width(), gs.get_legend_row_width());
	custom_plot::add_gradient_legend(draw_area, legend_layout, -3, 3, "NES", gs);
	custom_plot::add_title(draw_area, "GSEA Results", gs);

	this->draw_suite_->update(draw_area);

};

void MainWindow::s_multiple_group_enrichment_plot_top_n() {

	auto enrichments = this->signal_emitter_->get_type_variable_std(soap::VariableType::Enrichment);
	if (enrichments.size() < 2) {
		G_WARN("No enough enrichments found.");
		return;
	}

	auto enrichment_names = custom::keys(enrichments);

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"Plot Settings",
		{ "Enrichment:Label", "top n:5", "Capitalize First:yes"},
		{soap::InputStyle::SimpleChoiceWithLineEdit, soap::InputStyle::IntegerLineEdit,
		soap::InputStyle::SwitchButton},
		{enrichment_names}
	);

	if (settings.isEmpty()) {
		return;
	}

	auto [enrichment_use, label_use] = simple_choice_with_line_edit_layout_to_pair(settings[0]);
	if (enrichment_use.isEmpty()) {
		return;
	}

	int n_path = settings[1].toInt();
	if (n_path < 1) {
		return;
	}

	bool capitalize_first = switch_to_bool(settings[2]);

	QStringList groups;
	QList<QColor> group_colors, colors;
	QVector<double> counts;
	QStringList paths;

	QStringList final_labels = custom::ifelse<QStringList>(custom::sapply(label_use, [](auto&& l) {return l.isEmpty(); }), enrichment_use, label_use);

	auto&& gs = this->draw_suite_->graph_settings_;

	auto palette = gs.palette(final_labels);
	int n_enrichment_use = enrichment_use.size();
	for (int i = 0; i < n_enrichment_use; ++i) {
		auto enrichment = static_cast<Enrichment*>(enrichments[enrichment_use[i]]);
		auto&& enrichment_pathway_names = enrichment->mat_.get_const_qstring_reference(METADATA_ENRICHMENT_PATHWAY_NAMES);

		int enrich_size = enrichment_pathway_names.size();

		if (enrich_size == 0) {
			continue;
		}

		int n_sub_path = enrich_size >= n_path ? n_path : enrich_size;

		paths << enrichment_pathway_names.sliced(0, n_sub_path);
		groups << (label_use[i].isEmpty() ? enrichment_use[i] : label_use[i]);
		group_colors << palette[i];
		colors << QVector<QColor>(n_sub_path, palette[i]);
		counts << enrichment->mat_.get_double(METADATA_ENRICHMENT_COUNT).sliced(0, n_sub_path);
	}

	if (groups.isEmpty()) return;

	if (capitalize_first) {
		std::ranges::for_each(paths, [](auto&& path) {path = custom::capitalize_first(path); });
	}

	double maximum_length = std::ranges::max(counts);
	int bar_number = paths.size();
	QVector<double> x, y;
	QPen pen;
	pen.setWidth(20);

	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	axis_rect->axis(QCPAxis::atLeft)->setBasePen(QPen(Qt::NoPen));
	axis_rect->axis(QCPAxis::atLeft)->setTickLength(0);
	axis_rect->axis(QCPAxis::atLeft)->setSubTickPen(QPen(Qt::NoPen));
	axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);

	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);

	for (int i = 0; i < bar_number; ++i) {
		draw_area->addGraph(axis_rect->axis(QCPAxis::atLeft), axis_rect->axis(QCPAxis::atBottom));
		draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
		x.clear();
		y.clear();
		pen.setColor(colors[i]);
		x << bar_number - i;
		y << counts[i];
		draw_area->graph()->setPen(pen);
		draw_area->graph()->setData(x, y);
	}
	custom_plot::patch::set_range(axis_rect, QCPRange(0, maximum_length + 1), QCPRange(0, bar_number + 1));

	QSharedPointer<QCPAxisTickerText> y_ticker(new QCPAxisTickerText);
	x = custom::linspaced(bar_number, 1, bar_number);
	y_ticker->setTicks(QVector<double>(x.rbegin(), x.rend()), paths);
	axis_rect->axis(QCPAxis::atLeft)->setTicker(y_ticker);
	axis_rect->axis(QCPAxis::atLeft)->setTickLabelFont(gs.get_left_label_font());
	axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 0);
	axis_rect->axis(QCPAxis::atLeft)->setTickLabelRotation(gs.get_left_label_angle());

	pen.setColor(Qt::black);
	pen.setWidth(3);
	axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
	custom_plot::set_bottom_title(axis_rect, "Gene Number", gs, true);

	custom_plot::add_square_legend(draw_area, legend_layout, groups, group_colors, gs.get_legend_label_font(), gs.get_legend_title("Group"), gs.get_legend_title_font(), gs.get_legend_column_width(), gs.get_legend_row_width());

	custom_plot::add_title(draw_area, "Enrichment Results", gs);

	this->draw_suite_->update(draw_area);
};

void MainWindow::s_multiple_group_enrichment_plot_select_pathway() {
	auto enrichments = this->signal_emitter_->get_type_variable_std(soap::VariableType::Enrichment);
	if (enrichments.size() < 2) {
		G_WARN("No enough enrichments found.");
		return;
	}

	QMap<QString, QStringList> en2path;

	for (auto&& [name, data] : enrichments) {
		auto enrichment = static_cast<Enrichment*>(data);
		en2path[name] = custom::sorted(enrichment->mat_.get_qstring(METADATA_ENRICHMENT_PATHWAY_GENES));
	}
	auto [capitalize_first, settings] = EnrichmentMultiGroupPlotDialog::get_plot_settings(en2path);
	if (settings.isEmpty())return;

	auto& gs = this->draw_suite_->graph_settings_;
	bool palette_active = gs.active() && gs.bool_info_[IS_PALATTE_ACTIVE];
	auto palette = gs.palette_;

	QStringList groups;
	QList<QColor> group_colors, colors;
	QVector<double> counts;
	QStringList paths;
	for (const auto& group : settings) {
		auto [orig, alias] = factor_double_line_edit_with_completer_to_pair(group[0]);

		if (orig.second.isEmpty()) {
			continue;
		}

		QString group_name;
		if (alias.first.isEmpty()) {
			group_name = orig.first;
		}
		else {
			group_name = alias.first;
		}

		QColor group_color = QColor::fromString(group[1]);
		if (palette_active) {
			auto iter = palette.find(group_name);
			if (iter != palette.end()) {
				group_color = *iter;
			}
		}

		QStringList path_names = orig.second;
		auto ptr_enrichment = static_cast<Enrichment*>(enrichments[orig.first]);
		auto&& enrichment_pathway_names = ptr_enrichment->mat_.get_const_qstring_reference(METADATA_ENRICHMENT_PATHWAY_NAMES);

		auto filter = custom::in(path_names, enrichment_pathway_names);
		int n_valid_path = filter.count();
		if (n_valid_path == 0) {
			continue;
		}

		QStringList original_names = custom::sliced(path_names, filter);
		custom::assign(path_names, alias.second, custom::sapply(alias.second, [](auto&& name) {return !name.isEmpty(); }));

		path_names = custom::sliced(path_names, filter);
		paths << path_names;

		groups << group_name;
		group_colors << group_color;
		colors << QVector<QColor>(n_valid_path, group_color);
		counts << custom::reordered(ptr_enrichment->mat_.get_double(METADATA_ENRICHMENT_COUNT), custom::index_of(original_names, enrichment_pathway_names));
	}
	if (groups.isEmpty()) return;

	if (capitalize_first) {
		std::ranges::for_each(paths, [](auto&& path) {path = custom::capitalize_first(path); });
	}

	double maximum_length = std::ranges::max(counts);
	int bar_number = paths.size();
	QVector<double> x, y;
	QPen pen;
	pen.setWidth(20);

	auto [draw_area, axis_rect, legend_layout] = custom_plot::prepare(gs);

	axis_rect->axis(QCPAxis::atLeft)->setBasePen(QPen(Qt::NoPen));
	axis_rect->axis(QCPAxis::atLeft)->setTickLength(0);
	axis_rect->axis(QCPAxis::atLeft)->setSubTickPen(QPen(Qt::NoPen));
	axis_rect->axis(QCPAxis::atLeft)->grid()->setVisible(false);

	axis_rect->axis(QCPAxis::atBottom)->grid()->setVisible(false);

	for (int i = 0; i < bar_number; ++i) {
		draw_area->addGraph(axis_rect->axis(QCPAxis::atLeft), axis_rect->axis(QCPAxis::atBottom));
		draw_area->graph()->setLineStyle(QCPGraph::lsImpulse);
		x.clear();
		y.clear();
		pen.setColor(colors[i]);
		x << bar_number - i;
		y << counts[i];
		draw_area->graph()->setPen(pen);
		draw_area->graph()->setData(x, y);
	}
	custom_plot::patch::set_range(axis_rect, QCPRange(0, maximum_length + 1), QCPRange(0, bar_number + 1));

	QSharedPointer<QCPAxisTickerText> y_ticker(new QCPAxisTickerText);
	x = custom::linspaced(bar_number, 1, bar_number);
	y_ticker->setTicks(QVector<double>(x.rbegin(), x.rend()), paths);
	axis_rect->axis(QCPAxis::atLeft)->setTicker(y_ticker);
	axis_rect->axis(QCPAxis::atLeft)->setTickLabelFont(gs.get_left_label_font());
	axis_rect->axis(QCPAxis::atLeft)->setTickLength(0, 0);
	axis_rect->axis(QCPAxis::atLeft)->setTickLabelRotation(gs.get_left_label_angle());

	pen.setColor(Qt::black);
	pen.setWidth(3);
	axis_rect->axis(QCPAxis::atBottom)->setTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setSubTickPen(pen);
	axis_rect->axis(QCPAxis::atBottom)->setBasePen(pen);
	custom_plot::set_bottom_title(axis_rect, "Gene Number", gs, true);

	custom_plot::add_square_legend(draw_area, legend_layout, groups, group_colors, gs.get_legend_label_font(), gs.get_legend_title("Group"), gs.get_legend_title_font(), gs.get_legend_column_width(), gs.get_legend_row_width());

	custom_plot::add_title(draw_area, "Enrichment Results", gs);

	this->draw_suite_->update(draw_area);
};

std::pair<QString, int> listFilesRecursively(const QString& currentDir, int tab_count) {

	QString res;

	QDirIterator it(currentDir, QDir::NoDotAndDotDot | QDir::AllEntries, QDirIterator::NoIteratorFlags);

	static int dircount{ 0 };
	static int compcount{ 0 };

	while (it.hasNext()) {
		it.next();
		QFileInfo fileInfo = it.fileInfo();
		QString fn = fileInfo.fileName();

		if (fileInfo.isFile()) {

			QString compid = "Component_" + QString::number(compcount++);
			QString fileid = fn;
			QString filename = fn;
			QString Source = fileInfo.absoluteFilePath();
			Source.replace("/", "\\\\");

			for (int i = 0; i < tab_count; ++i) {
				res += "\t";
			}

			res += "<Component Id=\"" + compid + "\" Guid=\"*\" >\n";

			for (int i = 0; i < tab_count + 1; ++i) {
				res += "\t";
			}

			res += "<File Id=\"" + fileid + "\" Name=\"" + filename + "\" Source=\"" + Source + "\" />\n";
			
			for (int i = 0; i < tab_count; ++i) {
				res += "\t";
			}
				
			res += "</Component>\n\n";
		}

		if (fileInfo.isDir()) {

			QString dirid = "dir" + QString::number(dircount++);
			QString dirname = fn;

			for (int i = 0; i < tab_count; ++i) {
				res += "\t";
			}

			res += "<Directory Id=\"" + dirid + "\" Name=\"" + dirname + "\" >\n";

			res += listFilesRecursively(fileInfo.absoluteFilePath(), tab_count + 1).first;

			for (int i = 0; i < tab_count; ++i) {
				res += "\t";
			}

			res += "</Directory>\n\n";
		}
	}

	return { res, compcount };
}