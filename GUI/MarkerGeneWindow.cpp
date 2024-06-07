#include "MarkerGeneWindow.h"
#include <QTableWidgetItem>
#include <QHeaderView>

#include "CellMarkerTreeWidget.h"
#include "Identifier.h"

#include "SoapGUI.h"

MarkerGeneWindow::MarkerGeneWindow(
	const QString& cell_type, 
	const QStringList& marker_names, 
	CellMarkerTreeWidget::MarkerWidgetType type,
	QWidget* parent) :
	QMainWindow(parent),
	marker_names_(marker_names)
{
	this->main_interface_ = new QWidget(this);
	this->setCentralWidget(this->main_interface_);
	this->main_interface_->setLayout(&this->main_layout_);

	QIcon eye_icon(FILE_EYE_ICON_PNG);

	int row_num = this->marker_names_.size();
	if (type == CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget) {
		this->matrix_table_ = new QTableWidget(row_num, 4, this->main_interface_);
		set_matrix_table();

		for (int i = 0; i < row_num; ++i) {
			this->matrix_table_->setItem(i, 0, new QTableWidgetItem(this->marker_names_[i]));

			G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30))
			connect(button0, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked1);

			G_SET_NEW_BUTTON(button1, "►", QSize(30, 30))
			connect(button1, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked2);

			G_SET_NEW_BUTTON(button2, "▼", QSize(30, 30))
			connect(button2, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked3);

			this->matrix_table_->setCellWidget(i, 1, button0);
			this->matrix_table_->setCellWidget(i, 2, button1);
			this->matrix_table_->setCellWidget(i, 3, button2);
			this->matrix_table_->setRowHeight(i, 30);

		}
		this->matrix_table_->setColumnWidth(1, 30);
		this->matrix_table_->setColumnWidth(2, 30);
		this->matrix_table_->setColumnWidth(3, 30);
	}
	else if (type == CellMarkerTreeWidget::MarkerWidgetType::AlternativeWidget) {
		this->matrix_table_ = new QTableWidget(row_num, 3, this->main_interface_);
		set_matrix_table();
		for (int i = 0; i < row_num; ++i) {
			this->matrix_table_->setItem(i, 0, new QTableWidgetItem(this->marker_names_[i]));

			G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30))
			connect(button0, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked1);

			G_SET_NEW_BUTTON(button1, "►", QSize(30, 30))
			connect(button1, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked3);

			this->matrix_table_->setCellWidget(i, 1, button0);
			this->matrix_table_->setCellWidget(i, 2, button1);
			this->matrix_table_->setRowHeight(i, 30);

		}
		this->matrix_table_->setColumnWidth(1, 30);
		this->matrix_table_->setColumnWidth(2, 30);
	}
	else if (type == CellMarkerTreeWidget::MarkerWidgetType::ActiveWidget) {
		this->matrix_table_ = new QTableWidget(row_num, 3, this->main_interface_);
		this->set_matrix_table();
		for (int i = 0; i < row_num; ++i) {
			this->matrix_table_->setItem(i, 0, new QTableWidgetItem(this->marker_names_[i]));

			G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30))
			connect(button0, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked1);

			G_SET_NEW_BUTTON(button1, "◄", QSize(30, 30))
			connect(button1, &QPushButton::clicked, this, &MarkerGeneWindow::s_marker_clicked2);

			this->matrix_table_->setCellWidget(i, 1, button0);
			this->matrix_table_->setCellWidget(i, 2, button1);
			this->matrix_table_->setRowHeight(i, 30);
		}
		this->matrix_table_->setColumnWidth(1, 30);
		this->matrix_table_->setColumnWidth(2, 30);
	}

	this->setAttribute(Qt::WA_DeleteOnClose);
	this->setWindowIcon(QIcon(":/soap/Image/mainwindowicon.jpg"));
	this->resize(300, 600);
	this->setWindowTitle("Markers of " + cell_type);
	this->show();
}

void MarkerGeneWindow::set_matrix_table() {
	this->main_layout_.addWidget(this->matrix_table_);
	this->matrix_table_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	this->matrix_table_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	this->matrix_table_->horizontalHeader()->setHidden(true);
};

void MarkerGeneWindow::s_marker_clicked1() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());
	QModelIndex index = this->matrix_table_->indexAt(btn->pos());
	emit x_show_gene(this->marker_names_[index.row()]);
};

void MarkerGeneWindow::s_marker_clicked2() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());
	QModelIndex index = this->matrix_table_->indexAt(btn->pos());
	emit x_to_alternative(this->marker_names_[index.row()]);
};

void MarkerGeneWindow::s_marker_clicked3() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());
	QModelIndex index = this->matrix_table_->indexAt(btn->pos());
	emit x_to_active(this->marker_names_[index.row()]);
};

