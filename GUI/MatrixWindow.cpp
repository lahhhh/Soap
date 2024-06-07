#include "MatrixWindow.h"

#include <QTableWidgetItem>

#include "SoapGUI.h"

#include "MessageDialog.h"
#include "CommonDialog.h"
#include "FileIO.h"

#define CHECK_COPY \
if (!copy_data) { \
connect(signal_emitter, &SignalEmitter::x_data_delete_soon, window, &MatrixWindow::s_check_data); \
connect(signal_emitter, &SignalEmitter::x_data_edit_soon, window, &MatrixWindow::s_check_data); \
}

MatrixWindow::MatrixWindow(
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from
) :
	QMainWindow(signal_emitter->widget_),
	signal_emitter_(signal_emitter),
	matrix_table_(nullptr),
	copy_data_(copy_data),
	from_(from)
{
	this->main_layout_ = new QHBoxLayout;

	this->main_interface_ = new QWidget(this);
	this->setCentralWidget(this->main_interface_);
	this->main_interface_->setLayout(this->main_layout_);

	this->menubar_ = new QMenuBar(this);
	this->setMenuBar(this->menubar_);

	auto menu_save = this->menubar_->addMenu("Save as...");
	auto menu_tsv = menu_save->addMenu("TSV");
	auto menu_csv = menu_save->addMenu("CSV");

	menu_tsv->addAction("This Page", this, &MatrixWindow::s_save_as_tsv_page);
	menu_tsv->addAction("All", this, &MatrixWindow::s_save_as_tsv_all);

	menu_csv->addAction("This Page", this, &MatrixWindow::s_save_as_csv_page);
	menu_csv->addAction("All", this, &MatrixWindow::s_save_as_csv_all);
	//------------------------------------------//

	this->right_panel_ = new QWidget();

	this->up_button_ = new QPushButton("▲", this);
	this->down_button_ = new QPushButton("▼", this);
	this->left_button_ = new QPushButton("◄", this);
	this->right_button_ = new QPushButton("►", this);

	this->up_button_->setFixedHeight(30);
	this->down_button_->setFixedHeight(30);
	this->left_button_->setFixedHeight(30);
	this->right_button_->setFixedHeight(30);

	QGridLayout* direction_layout = new QGridLayout;
	direction_layout->addWidget(this->up_button_, 0, 0, 1, 2); 
	direction_layout->addWidget(this->left_button_, 1, 0, 1, 1); 
	direction_layout->addWidget(this->right_button_, 1, 1, 1, 1); 
	direction_layout->addWidget(this->down_button_, 2, 0, 1, 2); 

	direction_layout->setRowStretch(0, 1);
	direction_layout->setRowStretch(1, 1);
	direction_layout->setRowStretch(2, 1);
	direction_layout->setColumnStretch(0, 1);
	direction_layout->setColumnStretch(1, 1);

	G_SET_LABEL(this->row_label_, "Row", QSize(60, 30));
	G_SET_LABEL(this->column_label_, "Column", QSize(60, 30));

	G_SET_LINEEDIT(this->row_line_edit_, "1", QSize(50, 30));
	G_SET_LINEEDIT(this->column_line_edit_, "1", QSize(50, 30));
	this->row_line_edit_->setValidator(new QIntValidator());
	this->column_line_edit_->setValidator(new QIntValidator());

	G_SET_LABEL(this->row_number_label_, "", QSize(100, 30));
	G_SET_LABEL(this->column_number_label_, "", QSize(100, 30));

	this->row_layout_ = new QHBoxLayout;
	this->row_layout_->addWidget(this->row_label_);
	this->row_layout_->addWidget(this->row_line_edit_);
	this->row_layout_->addWidget(this->row_number_label_);

	this->column_layout_ = new QHBoxLayout;
	this->column_layout_->addWidget(this->column_label_);
	this->column_layout_->addWidget(this->column_line_edit_);
	this->column_layout_->addWidget(this->column_number_label_);

	G_SET_BUTTON(this->go_button_, "Go", QSize(180, 30));

	this->go_layout_ = new QVBoxLayout;

	this->go_layout_->addLayout(this->row_layout_);
	this->go_layout_->addLayout(this->column_layout_);
	this->go_layout_->addWidget(this->go_button_);

	this->information_area_ = new InformationTextBrowser(this);

	this->right_layout_ = new QVBoxLayout;
	this->right_layout_->addLayout(direction_layout);
	this->right_layout_->addLayout(this->go_layout_);
	this->right_layout_->addWidget(this->information_area_);

	this->right_panel_->setLayout(this->right_layout_);

	this->right_panel_->setFixedWidth(220);
	//-------------------------------//
	this->table_layout_ = new QVBoxLayout;
	this->set_matrix_table();
	this->main_layout_->addLayout(this->table_layout_);
	this->main_layout_->addWidget(this->right_panel_);

	connect(this->up_button_, &QPushButton::clicked, this, &MatrixWindow::s_up_page);
	connect(this->down_button_, &QPushButton::clicked, this, &MatrixWindow::s_down_page);
	connect(this->left_button_, &QPushButton::clicked, this, &MatrixWindow::s_left_page);
	connect(this->right_button_, &QPushButton::clicked, this, &MatrixWindow::s_right_page);
	connect(this->go_button_, &QPushButton::clicked, this, &MatrixWindow::s_go_page);

	this->setAttribute(Qt::WA_DeleteOnClose);
	G_SET_ICON;
	this->setWindowTitle(title);
	this->resize(800, 800);
};

MatrixWindow::~MatrixWindow() {
	delete this->matrix_table_;
	for (auto ptr : this->actions_) {
		delete ptr.second;
	}
	if (this->copy_data_) {
		if (this->data_type_ == soap::TableType::EigenArrayXXd) {
			delete static_cast<Eigen::ArrayXXd*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::EigenArrayXXi) {
			delete static_cast<Eigen::ArrayXXi*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::CustomMatrix) {
			delete static_cast<CustomMatrix*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::EigenMatrixXi) {
			delete static_cast<Eigen::MatrixXi*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::EigenMatrixXd) {
			delete static_cast<Eigen::MatrixXd*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::EigenSparseMatrixDouble) {
			delete static_cast<Eigen::SparseMatrix<double>*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::EigenSparseMatrixInt) {
			delete static_cast<Eigen::SparseMatrix<int>*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::StringList) {
			delete static_cast<QStringList*>(this->data_);
		}
		else if (this->data_type_ == soap::TableType::GenomicRange) {
			delete static_cast<GenomicRange*>(this->data_);
		}
	}
}

void MatrixWindow::s_up_page() {
	if (this->current_row_start_ == 1)return;
	int row_start_target = this->current_row_start_ - this->height_;
	row_start_target = row_start_target < 1 ? 1 : row_start_target;
	this->update_table(row_start_target, this->current_column_start_);
};

void MatrixWindow::s_down_page() {
	if (this->current_row_end_ == this->nrow_)return;
	int row_start_target = this->current_row_end_ + 1;
	this->update_table(row_start_target, this->current_column_start_);
};

void MatrixWindow::s_left_page() {
	if (this->current_column_start_ == 1)return;
	int col_start_target = this->current_column_start_ - this->width_;
	col_start_target = col_start_target < 1 ? 1 : col_start_target;
	this->update_table(this->current_row_start_, col_start_target);
};

void MatrixWindow::s_right_page() {
	if (this->current_column_end_ == this->ncol_)return;
	int col_start_target = this->current_column_end_ + 1;
	this->update_table(this->current_row_start_, col_start_target);
};

void MatrixWindow::s_go_page() {
	int row_start_target = this->row_line_edit_->text().toInt();
	int col_start_target = this->column_line_edit_->text().toInt();
	this->update_table(row_start_target, col_start_target);
}

void MatrixWindow::update_table(int row_start_target, int col_start_target) {
	if (row_start_target > this->nrow_
		|| col_start_target > this->ncol_
		|| row_start_target < 0
		|| col_start_target < 0)
	{
		return;
	}
	this->matrix_table_->clear();

	this->current_row_start_ = row_start_target;
	this->current_row_end_ = row_start_target + this->height_ - 1;
	this->current_row_end_ = this->current_row_end_ > this->nrow_ ? this->nrow_ : this->current_row_end_;

	this->current_column_start_ = col_start_target;
	this->current_column_end_ = col_start_target + this->width_ - 1;
	this->current_column_end_ = this->current_column_end_ > this->ncol_ ? this->ncol_ : this->current_column_end_;

	this->matrix_table_->setRowCount(this->current_row_end_ - this->current_row_start_ + 1);
	this->matrix_table_->setColumnCount(this->current_column_end_ - this->current_column_start_ + 1);
	this->matrix_table_->setHorizontalHeaderLabels(this->colnames_.sliced(this->current_column_start_ - 1, this->current_column_end_ - this->current_column_start_ + 1));
	this->matrix_table_->setVerticalHeaderLabels(this->rownames_.sliced(this->current_row_start_ - 1, this->current_row_end_ - this->current_row_start_ + 1));

	if (this->data_type_ == soap::TableType::EigenArrayXXd) {
		auto mat = static_cast<Eigen::ArrayXXd*>(this->data_);
		for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
			for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->operator()(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::EigenArrayXXi) {
		auto mat = static_cast<Eigen::ArrayXXi*>(this->data_);
		for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
			for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->operator()(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::CustomMatrix) {
		auto mat = static_cast<CustomMatrix*>(this->data_);
		for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
			for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(mat->get_qstring(row - 1, col - 1)));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::EigenMatrixXi) {
		auto mat = static_cast<Eigen::MatrixXi*>(this->data_);
		for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
			for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->operator()(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::EigenMatrixXd) {
		auto mat = static_cast<Eigen::MatrixXd*>(this->data_);
		for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
			for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->operator()(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::EigenSparseMatrixDouble) {
		auto mat = static_cast<Eigen::SparseMatrix<double> *>(this->data_);
		for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
			for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->coeff(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::EigenSparseMatrixInt) {
		auto mat = static_cast<Eigen::SparseMatrix<int> *>(this->data_);
		for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
			for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->coeff(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::GenomicRange) {
		auto mat = static_cast<GenomicRange*>(this->data_);
		for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
			for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(mat->get_qstring(row - 1, col - 1)));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::MotifPosition) {
		auto mat = static_cast<MotifPosition*>(this->data_);
		for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
			for (int col = this->current_column_start_; col <= this->current_column_end_; ++col) {
				this->matrix_table_->setItem(row - this->current_row_start_, col - this->current_column_start_, new QTableWidgetItem(QString::number(mat->get_match_count(row - 1, col - 1))));
			}
		}
	}
	else if (this->data_type_ == soap::TableType::StringList) {
		auto mat = static_cast<QStringList*>(this->data_);
		for (int row = this->current_row_start_; row <= this->current_row_end_; ++row) {
			this->matrix_table_->setItem(row - this->current_row_start_, 0, new QTableWidgetItem(mat->at(row - 1)));
		}
	}
	this->row_line_edit_->setText(QString::number(this->current_row_start_));
	this->column_line_edit_->setText(QString::number(this->current_column_start_));
};

void MatrixWindow::set_matrix_table() {
	this->matrix_table_ = new QTableWidget();
	this->table_layout_->addWidget(this->matrix_table_);
	this->matrix_table_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	this->matrix_table_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	this->matrix_table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
	this->matrix_table_->setFont(QFont("Times New Roman", 15, 500));
};

void MatrixWindow::set_table_range() {

	this->row_number_label_->setText("(1~" + QString::number(this->nrow_) + ")");
	this->column_number_label_->setText("(1~" + QString::number(this->ncol_) + ")");
}

void MatrixWindow::show_mat(
	Eigen::ArrayXXd* data,
	const QStringList& rownames,
	const QStringList& colnames)
{
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::ArrayXXd(*data);
	}
	else {
		this->data_ = (void*)data;
	}
	this->data_type_ = soap::TableType::EigenArrayXXd;
	this->rownames_ = rownames;
	this->colnames_ = colnames;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(
	Eigen::MatrixXd* data,
	const QStringList& rownames,
	const QStringList& colnames)
{
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::MatrixXd(*data);
	}
	else {
		this->data_ = (void*)data;
	}
	this->data_type_ = soap::TableType::EigenMatrixXd;
	this->rownames_ = rownames;
	this->colnames_ = colnames;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(
	Eigen::ArrayXXi* data,
	const QStringList& rownames,
	const QStringList& colnames)
{
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::ArrayXXi(*data);
	}
	else {
		this->data_ = (void*)data;
	}
	this->data_type_ = soap::TableType::EigenArrayXXi;
	this->rownames_ = rownames;
	this->colnames_ = colnames;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(
	Eigen::MatrixXi* data,
	const QStringList& rownames,
	const QStringList& colnames)
{
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::MatrixXi(*data);
	}
	else {
		this->data_ = (void*)data;
	}
	this->data_type_ = soap::TableType::EigenMatrixXi;
	this->rownames_ = rownames;
	this->colnames_ = colnames;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(QStringList* data) {

	if (this->copy_data_) {
		this->data_ = (void*)new QStringList(*data);
	}
	else {
		this->data_ = (void*)data;
	}

	this->data_type_ = soap::TableType::StringList;
	this->rownames_ = _Cs cast<QString>(_Cs seq_n(0, data->size()));
	this->colnames_ = { "0" };
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(CustomMatrix* mat) {
	if (this->copy_data_) {
		this->data_ = (void*)new CustomMatrix(*mat);
	}
	else {
		this->data_ = (void*)mat;
	}
	this->data_type_ = soap::TableType::CustomMatrix;
	this->rownames_ = mat->rownames_;
	this->colnames_ = mat->colnames_;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(
	MotifPosition* data,
	const QStringList& rownames,
	const QStringList& colnames)
{
	if (this->copy_data_) {
		this->data_ = (void*)new MotifPosition(*data);
	}
	else {
		this->data_ = (void*)data;
	}
	this->data_type_ = soap::TableType::MotifPosition;
	this->rownames_ = rownames;
	this->colnames_ = colnames;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();

};

void MatrixWindow::show_mat(
	GenomicRange* data,
	const QStringList& rownames,
	const QStringList& colnames)
{
	if (this->copy_data_) {
		this->data_ = (void*)new GenomicRange(*data);
	}
	else {
		this->data_ = (void*)data;
	}

	this->data_type_ = soap::TableType::GenomicRange;
	this->rownames_ = rownames;
	this->colnames_ = colnames;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(SparseInt* data) {
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::SparseMatrix<int>(data->mat_);
	}
	else {
		this->data_ = (void*)&data->mat_;
	}
	this->data_type_ = soap::TableType::EigenSparseMatrixInt;
	this->rownames_ = data->rownames_;
	this->colnames_ = data->colnames_;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(SparseDouble* data) {
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::SparseMatrix<double>(data->mat_);
	}
	else {
		this->data_ = (void*)&data->mat_;
	}
	this->data_type_ = soap::TableType::EigenSparseMatrixDouble;
	this->rownames_ = data->rownames_;
	this->colnames_ = data->colnames_;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(DenseDouble* data) {

	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::MatrixXd(data->mat_);
	}
	else {
		this->data_ = (void*)&data->mat_;
	}
	this->data_type_ = soap::TableType::EigenMatrixXd;
	this->rownames_ = data->rownames_;
	this->colnames_ = data->colnames_;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(DenseInt* data) {

	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::MatrixXi(data->mat_);
	}
	else {
		this->data_ = (void*)&data->mat_;
	}
	this->data_type_ = soap::TableType::EigenMatrixXi;
	this->rownames_ = data->rownames_;
	this->colnames_ = data->colnames_;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_mat(Embedding* data) {
	if (this->copy_data_) {
		this->data_ = (void*)new Eigen::MatrixXd(data->data_.mat_);
	}
	else {
		this->data_ = (void*)&data->data_.mat_;
	}
	this->data_type_ = soap::TableType::EigenMatrixXd;
	this->rownames_ = data->data_.rownames_;
	this->colnames_ = data->data_.colnames_;
	this->nrow_ = this->rownames_.size();
	this->ncol_ = this->colnames_.size();

	this->set_table_range();
	this->update_table(1, 1);
	this->show();
};

void MatrixWindow::show_matrix(
	SparseInt* data,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data);
};

void MatrixWindow::show_matrix(
	SparseDouble* data,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data);
};

void MatrixWindow::show_matrix(
	DenseInt* data,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data);
};

void MatrixWindow::show_matrix(
	DenseDouble* data,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data);
};

void MatrixWindow::show_matrix(
	Embedding* data,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data);
};


void MatrixWindow::show_matrix(
	CustomMatrix* data,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data);
};

void MatrixWindow::show_matrix(
	GenomicRange* data,
	const QStringList& rownames,
	const QStringList& colnames,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data, rownames, colnames);
};

void MatrixWindow::show_matrix(
	MotifPosition* data,
	const QStringList& rownames,
	const QStringList& colnames,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(data, rownames, colnames);
};

void MatrixWindow::show_matrix(
	Eigen::ArrayXXd* mat,
	const QStringList& rownames,
	const QStringList& colnames,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(mat, rownames, colnames);
};

void MatrixWindow::show_matrix(
	Eigen::MatrixXd* mat,
	const QStringList& rownames,
	const QStringList& colnames,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(mat, rownames, colnames);
};

void MatrixWindow::show_matrix(
	Eigen::MatrixXd* mat,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	int nrow = mat->rows();
	int ncol = mat->cols();

	show_matrix(mat, _Cs cast<QString>(_Cs seq_n(0, nrow)), _Cs cast<QString>(_Cs seq_n(0, ncol)), title, signal_emitter, copy_data, from);
};

void MatrixWindow::show_matrix(
	Eigen::ArrayXXi* mat,
	const QStringList& rownames,
	const QStringList& colnames,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(mat, rownames, colnames);
};

void MatrixWindow::show_matrix(
	Eigen::MatrixXi* mat,
	const QStringList& rownames,
	const QStringList& colnames,
	const QString& title,
	SignalEmitter* signal_emitter,
	bool copy_data,
	void* from)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, copy_data, from);
	CHECK_COPY;
	window->show_mat(mat, rownames, colnames);
};

void MatrixWindow::show_matrix(
	QStringList* data,
	const QString& title,
	SignalEmitter* signal_emitter
)
{
	MatrixWindow* window = new MatrixWindow(title, signal_emitter, true, nullptr);
	window->show_mat(data);
};

void MatrixWindow::s_save_as_csv_all() {
	QString file_path = QFileDialog::getSaveFileName(this, "Please Set Table Name", "", "CSV(*.csv)");
	if (file_path.isEmpty())return;

	QFile file(file_path);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		G_WARN("Cannot open create file.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"CSV Settings",
		{ "quote:yes", "keep rownames:yes", "keep colnames:yes" },
		{ soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton }
	);

	if (settings.isEmpty()) {
		return;
	}

	bool quote = switch_to_bool(settings[0]);
	bool keep_rownames = switch_to_bool(settings[1]);
	bool keep_colnames = switch_to_bool(settings[2]);

	QTextStream stream(&file);

	if (this->data_type_ == soap::TableType::EigenArrayXXd) {
		auto mat = static_cast<Eigen::ArrayXXd*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenArrayXXi) {
		auto mat = static_cast<Eigen::ArrayXXi*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::CustomMatrix) {
		auto mat = static_cast<CustomMatrix*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames,
			keep_colnames);
	}
	else if (this->data_type_ == soap::TableType::EigenMatrixXi) {
		auto mat = static_cast<Eigen::MatrixXi*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenMatrixXd) {
		auto mat = static_cast<Eigen::MatrixXd*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenSparseMatrixDouble) {
		auto mat = static_cast<Eigen::SparseMatrix<double> *>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenSparseMatrixInt) {
		auto mat = static_cast<Eigen::SparseMatrix<int> *>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::GenomicRange) {
		auto mat = static_cast<GenomicRange*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::MotifPosition) {
		auto mat = static_cast<MotifPosition*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::StringList) {
		auto mat = static_cast<QStringList*>(this->data_);
		write_sv(*mat, stream, ',', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}

	stream.flush();

	G_LOG("CSV Writing finished.");
}

void MatrixWindow::s_save_as_tsv_all() {

	QString file_path = QFileDialog::getSaveFileName(this, "Please Set Table Name", "", "TSV(*.tsv)");
	if (file_path.isEmpty())return;

	QFile file(file_path);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		G_WARN("Cannot open create file.");
		return;
	}

	auto settings = CommonDialog::get_response(
		this->signal_emitter_,
		"TSV Settings",
		{ "quote:yes", "keep rownames:yes", "keep colnames:yes" },
		{ soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton, soap::InputStyle::SwitchButton }
	);

	if (settings.isEmpty()) {
		return;
	}

	bool quote = switch_to_bool(settings[0]);
	bool keep_rownames = switch_to_bool(settings[1]);
	bool keep_colnames = switch_to_bool(settings[2]);

	QTextStream stream(&file);

	if (this->data_type_ == soap::TableType::EigenArrayXXd) {
		auto mat = static_cast<Eigen::ArrayXXd*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenArrayXXi) {
		auto mat = static_cast<Eigen::ArrayXXi*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::CustomMatrix) {
		auto mat = static_cast<CustomMatrix*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames,
			keep_colnames);
	}
	else if (this->data_type_ == soap::TableType::EigenMatrixXi) {
		auto mat = static_cast<Eigen::MatrixXi*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenMatrixXd) {
		auto mat = static_cast<Eigen::MatrixXd*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenSparseMatrixDouble) {
		auto mat = static_cast<Eigen::SparseMatrix<double> *>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::EigenSparseMatrixInt) {
		auto mat = static_cast<Eigen::SparseMatrix<int> *>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::GenomicRange) {
		auto mat = static_cast<GenomicRange*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::MotifPosition) {
		auto mat = static_cast<MotifPosition*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}
	else if (this->data_type_ == soap::TableType::StringList) {
		auto mat = static_cast<QStringList*>(this->data_);
		write_sv(*mat, stream, '\t', quote,
			keep_rownames ? this->rownames_ : QStringList{},
			keep_colnames ? this->colnames_ : QStringList{});
	}

	stream.flush();

	G_LOG("TSV Writing finished.");
};

void MatrixWindow::s_save_as_tsv_page() {
	QString file_path = QFileDialog::getSaveFileName(this, "Please Set Table Name", "", "TSV(*.tsv)");
	if (file_path.isEmpty())return;

	QFile file(file_path);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		G_WARN("Cannot open create file.");
		return;
	}

	QString ret;
	int row_number = this->matrix_table_->rowCount();
	int column_number = this->matrix_table_->columnCount();
	for (int i = 0; i < column_number; ++i) {
		ret.append("\t");
		ret.append("\"" + this->matrix_table_->horizontalHeaderItem(i)->text() + "\"");
	}
	ret.append("\n");
	for (int i = 0; i < row_number; ++i) {
		ret.append("\"" + this->matrix_table_->verticalHeaderItem(i)->text() + "\"");
		for (int j = 0; j < column_number; ++j) {
			ret.append("\t");
			ret.append("\"" + this->matrix_table_->item(i, j)->text() + "\"");
		}
		ret.append("\n");
	}
	QTextStream stream(&file);
	stream << ret;
	stream.flush();

	G_LOG("TSV Writing finished.");
};

void MatrixWindow::s_save_as_csv_page() {
	QString file_path = QFileDialog::getSaveFileName(this, "Please Set Table Name", "", "CSV(*.csv)");
	if (file_path.isEmpty())return;

	QFile file(file_path);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		G_WARN("Cannot open create file.");
		return;
	}
	QString ret;
	int row_number = this->matrix_table_->rowCount();
	int column_number = this->matrix_table_->columnCount();
	for (int i = 0; i < column_number; ++i) {
		ret.append(",");
		ret.append("\"" + this->matrix_table_->horizontalHeaderItem(i)->text() + "\"");
	}
	ret.append("\n");
	for (int i = 0; i < row_number; ++i) {
		ret.append("\"" + this->matrix_table_->verticalHeaderItem(i)->text() + "\"");
		for (int j = 0; j < column_number; ++j) {
			ret.append(",");
			ret.append("\"" + this->matrix_table_->item(i, j)->text() + "\"");
		}
		ret.append("\n");
	}
	QTextStream stream(&file);
	stream << ret;
	stream.flush();

	G_LOG("CSV Writing finished.");
};

void MatrixWindow::s_check_data(void* from, soap::VariableType, void* item) {

	if (this->from_ == from || this->from_ == item) {

		this->close();
	}
};