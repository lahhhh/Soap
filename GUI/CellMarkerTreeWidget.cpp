#include "CellMarkerTreeWidget.h"

#include <QHeaderView>

#include "MarkerGeneWindow.h"
#include "Identifier.h"

#include "SoapGUI.h"


CellMarkerTreeWidget::CellMarkerTreeWidget(CellMarkerTreeWidget::MarkerWidgetType type, QWidget* parent) :
	QTreeWidget(parent),
	type_(type)
{
	connect(this, &QTreeWidget::itemDoubleClicked, this, &CellMarkerTreeWidget::QTreeWidgetItemDoubleClicked);

	this->setContextMenuPolicy(Qt::CustomContextMenu);
	this->headerItem()->setHidden(true);
	this->setColumnCount(4);

	QTreeWidgetItem* item = new QTreeWidgetItem();
	this->addTopLevelItem(item);
	this->header()->setSectionResizeMode(QHeaderView::ResizeToContents);
}

void CellMarkerTreeWidget::QTreeWidgetItemDoubleClicked(QTreeWidgetItem* tree, int) {

	if (this->currentItem() == nullptr)return;

	QModelIndex index = this->indexFromItem(tree);
	int row = index.row() - 1;
	if (row < 0)return;

	QString cell_type = tree->text(0);
	MarkerGeneWindow* window;

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget) {
		window = new MarkerGeneWindow(cell_type, this->type_to_markers_[cell_type], this->type_, this);
	}
	else {
		window = new MarkerGeneWindow(cell_type, this->markers_[row], this->type_, this);
	}

	connect(window, &MarkerGeneWindow::x_to_alternative, this, &CellMarkerTreeWidget::x_gene_to_alternative);
	connect(window, &MarkerGeneWindow::x_to_active, this, &CellMarkerTreeWidget::x_gene_to_active);
	connect(window, &MarkerGeneWindow::x_show_gene, this, &CellMarkerTreeWidget::x_show_gene);
}

void CellMarkerTreeWidget::safe_clear() {
	int count = this->topLevelItemCount();
	this->addTopLevelItem(new QTreeWidgetItem());

	for (int i = 0; i < count; ++i) {
		delete this->topLevelItem(0);
	}

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget) {
		this->type_to_markers_.clear();
	}
	else {
		this->types_.clear();
		this->markers_.clear();
	}
};

void CellMarkerTreeWidget::set_types(const QMap<QString, QStringList>& type_to_markers) {

	this->safe_clear();
	this->type_to_markers_ = type_to_markers;
	int type_size = type_to_markers.size();
	auto types = type_to_markers.keys();

	QIcon eye_icon(FILE_EYE_ICON_PNG);

	for (int i = 0; i < type_size; ++i) {

		QTreeWidgetItem* item = new QTreeWidgetItem({ types[i] });

		G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30));
		connect(button0, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_1);

		G_SET_NEW_BUTTON(button1, "►", QSize(30, 30));
		connect(button1, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_2);

		G_SET_NEW_BUTTON(button2, "▼", QSize(30, 30));
		connect(button2, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_3);

		this->addTopLevelItem(item);
		this->setItemWidget(item, 1, button0);
		this->setItemWidget(item, 2, button1);
		this->setItemWidget(item, 3, button2);
	}

	this->setColumnWidth(1, 30);
	this->setColumnWidth(2, 30);
	this->setColumnWidth(3, 30);
};

void CellMarkerTreeWidget::s_receive_type(const QString& type, const QStringList& markers) {
	this->types_ << type;
	this->markers_ << markers;

	QIcon eye_icon(FILE_EYE_ICON_PNG);
	QIcon delete_icon(FILE_DELETE_ICON_PNG);

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::AlternativeWidget) {
		QTreeWidgetItem* item = new QTreeWidgetItem({ type });

		G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30));
		connect(button0, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_1);

		G_SET_NEW_BUTTON(button1, "▼", QSize(30, 30));
		connect(button1, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_3);

		G_SET_NEW_BUTTON_ICON(button2, delete_icon, QSize(30, 30));
		connect(button2, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_4);

		this->addTopLevelItem(item);

		this->setItemWidget(item, 1, button0);
		this->setItemWidget(item, 2, button1);
		this->setItemWidget(item, 3, button2);
	}
	else if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::ActiveWidget) {
		QTreeWidgetItem* item = new QTreeWidgetItem({ type });

		G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30));
		connect(button0, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_1);

		G_SET_NEW_BUTTON(button1, "▲", QSize(30, 30));
		connect(button1, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_2);

		G_SET_NEW_BUTTON_ICON(button2, delete_icon, QSize(30, 30));
		connect(button2, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_4);

		this->addTopLevelItem(item);

		this->setItemWidget(item, 1, button0);
		this->setItemWidget(item, 2, button1);
		this->setItemWidget(item, 3, button2);
	}

	this->setColumnWidth(1, 30);
	this->setColumnWidth(2, 30);
	this->setColumnWidth(3, 30);
};

void CellMarkerTreeWidget::s_receive_marker(const QString& gene) {
	this->types_ << gene;
	this->markers_ << QStringList{ gene };

	QIcon eye_icon(FILE_EYE_ICON_PNG);

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::AlternativeWidget) {
		QTreeWidgetItem* item = new QTreeWidgetItem({ gene });

		G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30));
		connect(button0, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_1);

		G_SET_NEW_BUTTON(button1, "▼", QSize(30, 30));
		connect(button1, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_3);

		G_SET_NEW_BUTTON_ICON(button2, FILE_DELETE_ICON_PNG, QSize(30, 30));
		connect(button2, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_4);

		this->addTopLevelItem(item);

		this->setItemWidget(item, 1, button0);
		this->setItemWidget(item, 2, button1);
		this->setItemWidget(item, 3, button2);
	}
	else if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::ActiveWidget) {
		QTreeWidgetItem* item = new QTreeWidgetItem({ gene });

		G_SET_NEW_BUTTON_ICON(button0, eye_icon, QSize(30, 30));
		connect(button0, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_1);

		G_SET_NEW_BUTTON(button1, "▲", QSize(30, 30));
		connect(button1, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_2);

		G_SET_NEW_BUTTON_ICON(button2, FILE_DELETE_ICON_PNG, QSize(30, 30));
		connect(button2, &QPushButton::clicked, this, &CellMarkerTreeWidget::s_type_clicked_4);

		this->addTopLevelItem(item);

		this->setItemWidget(item, 1, button0);
		this->setItemWidget(item, 2, button1);
		this->setItemWidget(item, 3, button2);
	}
	this->setColumnWidth(1, 30);
	this->setColumnWidth(2, 30);
	this->setColumnWidth(3, 30);
};

void CellMarkerTreeWidget::s_type_clicked_1() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());

	QModelIndex index = this->indexAt(btn->pos());
	int row = index.row() - 1;

	QString cell_type;
	QStringList markers;

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget) {
		cell_type = this->type_to_markers_.keys()[index.row() - 1];
		markers = this->type_to_markers_[cell_type];
	}
	else {
		cell_type = this->types_[row];
		markers = this->markers_[row];
	}
	emit x_show_type(cell_type, markers);
};

void CellMarkerTreeWidget::s_type_clicked_2() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());

	QModelIndex index = this->indexAt(btn->pos());
	int row = index.row() - 1;

	QString cell_type;
	QStringList markers;

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget) {
		cell_type = this->type_to_markers_.keys()[index.row() - 1];
		markers = this->type_to_markers_[cell_type];
	}
	else {
		cell_type = this->types_[row];
		markers = this->markers_[row];
	}
	emit x_type_to_alternative(cell_type, markers);
};

void CellMarkerTreeWidget::s_type_clicked_3() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());

	QModelIndex index = this->indexAt(btn->pos());
	int row = index.row() - 1;

	QString cell_type;
	QStringList markers;

	if (this->type_ == CellMarkerTreeWidget::MarkerWidgetType::MarkerWidget) {
		cell_type = this->type_to_markers_.keys()[index.row() - 1];
		markers = this->type_to_markers_[cell_type];
	}
	else {
		cell_type = this->types_[row];
		markers = this->markers_[row];
	}
	emit x_type_to_active(cell_type, markers);
};

void CellMarkerTreeWidget::s_type_clicked_4() {
	QPushButton* btn = dynamic_cast<QPushButton*>(QObject::sender());

	QTreeWidgetItem* item = this->itemAt(btn->pos());

	int row = this->indexOfTopLevelItem(item) - 1;

	this->types_.removeAt(row);
	this->markers_.removeAt(row);
	delete item;
};
