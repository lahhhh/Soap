#pragma once

#include "Identifier.h"

#include <QMainWindow>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QMenuBar>
#include <QHeaderView>
#include <QLabel>
#include <QFileDialog>
#include <QMessageBox>

#include "SignalEmitter.h"
#include "InformationTextBrowser.h"

#include "Embedding.h"
#include "Metadata.h"
#include "SparseInt.h"
#include "SparseDouble.h"
#include "DenseInt.h"
#include "DenseDouble.h"
#include "Custom.h"
#include "GenomicRange.h"
#include "MotifPosition.h"

class MatrixWindow : public QMainWindow
{
	Q_OBJECT
public:
	MatrixWindow(
		const QString& title, 
		SignalEmitter* signal_emitter,
		bool copy_data, 
		void* from = nullptr
	);

private:
	void set_matrix_table();
	void set_table_range();
	void update_table(int row_start_target, int col_start_target);

	void show_mat(SparseInt* data);
	void show_mat(SparseDouble* data);
	void show_mat(DenseInt* data);
	void show_mat(DenseDouble* data);
	void show_mat(Embedding* data);
	void show_mat(CustomMatrix* data);
	void show_mat(Eigen::ArrayXXd* data, const QStringList& rownames, const QStringList& colnames);
	void show_mat(Eigen::ArrayXXi* data, const QStringList& rownames, const QStringList& colnames);
	void show_mat(Eigen::MatrixXd* data, const QStringList& rownames, const QStringList& colnames);
	void show_mat(Eigen::MatrixXi* data, const QStringList& rownames, const QStringList& colnames);
	void show_mat(QStringList* data);
	void show_mat(GenomicRange* data, const QStringList& rownames, const QStringList& colnames);	
	void show_mat(MotifPosition* data, const QStringList& rownames, const QStringList& colnames);

	SignalEmitter* signal_emitter_;
	InformationTextBrowser* information_area_{ nullptr };

	QHBoxLayout* main_layout_;
	QVBoxLayout* table_layout_;
	QVBoxLayout* right_layout_;
	QVBoxLayout* go_layout_;
	QHBoxLayout* row_layout_;
	QHBoxLayout* column_layout_;
	QHBoxLayout* left_right_layout_;

	QMenuBar* menubar_;
	std::map<QString, QAction* > actions_;

	QTableWidget* matrix_table_;

	QPushButton* up_button_;
	QPushButton* down_button_;
	QPushButton* left_button_;
	QPushButton* right_button_;

	QLabel* row_label_;
	QLabel* column_label_;

	QLabel* row_number_label_;
	QLabel* column_number_label_;

	QLineEdit* row_line_edit_;
	QLineEdit* column_line_edit_;

	QPushButton* go_button_;

	QWidget* main_interface_{ nullptr };
	QWidget* right_panel_{ nullptr };

	soap::TableType data_type_;

	static constexpr int width_ = 20;
	static constexpr int height_ = 50;

	void* data_ = nullptr;

	void* from_ = nullptr;

	int nrow_ = 0;
	int ncol_ = 0;

	int current_row_start_ = 1;
	int current_column_start_ = 1;

	int current_row_end_;
	int current_column_end_;

	QStringList rownames_;
	QStringList colnames_;

	bool copy_data_ = true;

public:

	~MatrixWindow();

	static void show_matrix(
		SparseInt* data, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		SparseDouble* data, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		DenseInt* data,
		const QString& title,
		SignalEmitter* signal_emitter,
		bool copy_data = true,
		void* from = nullptr);

	static void show_matrix(
		DenseDouble* data,
		const QString& title,
		SignalEmitter* signal_emitter,
		bool copy_data = true,
		void* from = nullptr);

	static void show_matrix(
		Embedding* data, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		CustomMatrix* mat_item, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		Eigen::ArrayXXd* data,
		const QStringList& rownames, 
		const QStringList& colnames, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		Eigen::ArrayXXi* data, 
		const QStringList& rownames, 
		const QStringList& colnames, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		Eigen::MatrixXd* data,
		const QStringList& rownames,
		const QStringList& colnames,
		const QString& title,
		SignalEmitter* signal_emitter,
		bool copy_data = true,
		void* from = nullptr);

	static void show_matrix(
		Eigen::MatrixXd* data,
		const QString& title,
		SignalEmitter* signal_emitter,
		bool copy_data = true,
		void* from = nullptr);

	static void show_matrix(
		Eigen::MatrixXi* data,
		const QStringList& rownames,
		const QStringList& colnames,
		const QString& title,
		SignalEmitter* signal_emitter,
		bool copy_data = true,
		void* from = nullptr);

	static void show_matrix(
		QStringList* data, 
		const QString& title, 
		SignalEmitter* signal_emitter);

	static void show_matrix(
		GenomicRange* data, 
		const QStringList& rownames, 
		const QStringList& colnames, 
		const QString& title, 
		SignalEmitter* signal_emitter, 
		bool copy_data = true, 
		void* from = nullptr);

	static void show_matrix(
		MotifPosition* data,
		const QStringList& rownames,
		const QStringList& colnames,
		const QString& title,
		SignalEmitter* signal_emitter,
		bool copy_data = true,
		void* from = nullptr);

private slots:

	void s_save_as_tsv_page();
	void s_save_as_csv_page();

	void s_save_as_tsv_all();
	void s_save_as_csv_all();

	void s_up_page();
	void s_down_page();
	void s_left_page();
	void s_right_page();

	void s_go_page();

	void s_check_data(void* from, soap::VariableType, void* item);
};
