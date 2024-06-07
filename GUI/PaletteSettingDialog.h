#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QComboBox>
#include <QPushButton>
#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>

class PaletteSettingDialog
	: public QDialog
{
	Q_OBJECT
public:

	explicit PaletteSettingDialog(const QMap<QString, QColor>& palette);

	static void get_response(QMap<QString, QColor>& palette);

	bool is_accepted_ = false;

	qsizetype nrow_ = 0;
	qsizetype ncol_ = 0;
	qsizetype last_column_row_ = 0;

	QVBoxLayout* all_layout_ = nullptr;
	QHBoxLayout* main_layout_ = nullptr;

	QList<qsizetype> rows_per_column_;

	QList<QVBoxLayout*> column_layouts_;
	QList<QHBoxLayout*> row_layouts_;
	QList<QLineEdit*> factor_line_edits_;
	QList<QLabel*> color_view_labels_;
	QList<QPushButton*> choose_color_buttons_;
	QList<QPushButton*> delete_buttons_;

	QHBoxLayout* add_layout_ = nullptr;

	QPushButton* finish_button_ = nullptr;
	QPushButton* cancel_button_ = nullptr;

private slots:
	void accept();
	void reject();

	void s_delete_palette();

	void s_add_palette();

	void s_choose_color();
};

