
#include "PaletteSettingDialog.h"

#include <QColorDialog>

#include "Identifier.h"

#include "SoapGUI.h"

PaletteSettingDialog::PaletteSettingDialog(const QMap<QString, QColor>& palette) {

	this->all_layout_ = new QVBoxLayout;
	this->main_layout_ = new QHBoxLayout;
	const qsizetype size = palette.size();

	if (size > 0) {
		this->ncol_ = size / 10 + 1;
		this->nrow_ = ceil((double)size / this->ncol_);
		if (this->nrow_ < 10) {
			this->nrow_ = 10;
		}
		this->ncol_ = ceil((double)size / this->nrow_);

		this->rows_per_column_.resize(this->ncol_, this->nrow_);

		this->last_column_row_ = this->rows_per_column_[this->ncol_ - 1] = size - (this->ncol_ - 1) * this->nrow_;
	}
	else {
		this->nrow_ = 10;
	}

	this->row_layouts_.reserve(size);
	this->factor_line_edits_.reserve(size);
	this->color_view_labels_.reserve(size);
	this->choose_color_buttons_.reserve(size);
	this->delete_buttons_.reserve(size);

	this->column_layouts_.reserve(this->ncol_);

	auto iter = palette.constKeyValueBegin();

	for (qsizetype col = 0; col < this->ncol_; col++) {

		QVBoxLayout* column_layout = new QVBoxLayout;

		this->column_layouts_ << column_layout;

		for (qsizetype i = 0; i < this->rows_per_column_[col]; ++i) {

			G_SET_NEW_LINEEDIT(line_edit, iter->first, soap::MiddleSize);
			G_SET_NEW_LABEL(label, iter->second.name(), QSize(80, 30));
			G_SET_LABEL_COLOR(label, iter->second);
			G_SET_NEW_BUTTON(choose_color_button, "Choose Color", QSize(100, 30));
			G_SET_NEW_BUTTON(delete_button, "-", QSize(30, 30));
			G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(column_layout, row_layout, line_edit, label, choose_color_button, delete_button);
			connect(delete_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_delete_palette);
			connect(choose_color_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_choose_color);

			this->row_layouts_ << row_layout;
			this->factor_line_edits_ << line_edit;
			this->color_view_labels_ << label;
			this->choose_color_buttons_ << choose_color_button;
			this->delete_buttons_ << delete_button;

			++iter;
		}

		this->main_layout_->addLayout(column_layout);
	}

	G_SET_NEW_BUTTON(add_button, "+ palette", soap::MiddleSize);
	G_SINGLE_ITEM_ROWLAYOUT(this->add_layout_, add_button);

	connect(add_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_add_palette);

	if (this->last_column_row_ < this->nrow_) {

		QVBoxLayout* column_layout;
		if (this->column_layouts_.size() > 0) {
			column_layout = this->column_layouts_.back();
		}
		else {
			column_layout = new QVBoxLayout;
			this->column_layouts_ << column_layout;
			this->main_layout_->addLayout(column_layout);
			++this->ncol_;
			this->rows_per_column_.resize(1, 0);
		}
		column_layout->addLayout(this->add_layout_);

		++this->last_column_row_;
		++this->rows_per_column_[this->ncol_ - 1];
	}
	else {
		QVBoxLayout* column_layout = new QVBoxLayout;

		this->column_layouts_ << column_layout;
		column_layout->addLayout(this->add_layout_);
		this->main_layout_->addLayout(column_layout);

		++this->ncol_;
		this->last_column_row_ = 1;
		this->rows_per_column_ << 1;
	}

	this->all_layout_->addLayout(this->main_layout_);

	G_SET_FINISH_BUTTON;

	G_SET_CANCEL_BUTTON;

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT(this->all_layout_, row_layout, this->finish_button_, this->cancel_button_);


	setLayout(this->all_layout_);

	connect(this->finish_button_, &QPushButton::clicked, this, &PaletteSettingDialog::accept);
	connect(this->cancel_button_, &QPushButton::clicked, this, &PaletteSettingDialog::reject);

	G_SET_ICON;
	this->setWindowTitle("Palette Setting");
	this->exec();
}

void PaletteSettingDialog::s_choose_color() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());

	qsizetype index = this->choose_color_buttons_.indexOf(button);

	QColor color = QColorDialog::getColor(
		QColor(this->color_view_labels_[index]->text()), 
		nullptr, 
		"Choose Color");

	if (color.isValid()) {
		this->color_view_labels_[index]->setText(color.name());
		this->color_view_labels_[index]->setStyleSheet("QLabel{color:" + color.name() + "}");
	}
};

void PaletteSettingDialog::get_response(QMap<QString, QColor>& palette) {
	PaletteSettingDialog dlg(palette);

	if (!dlg.is_accepted_) {
		return ;
	}

	palette.clear();

	const qsizetype size = dlg.row_layouts_.size();

	for (qsizetype i = 0; i < size; ++i) {
		palette[dlg.factor_line_edits_[i]->text()] = QColor(dlg.color_view_labels_[i]->text());
	}

};

void PaletteSettingDialog::s_add_palette() {
	if (this->last_column_row_ < this->nrow_) {

		QVBoxLayout* column_layout = this->column_layouts_.back();
		column_layout->removeItem(this->add_layout_);

		G_SET_NEW_LINEEDIT(line_edit, "", soap::MiddleSize);
		G_SET_NEW_LABEL(label, "", QSize(80, 30));
		G_SET_NEW_BUTTON(choose_color_button, "Choose Color", QSize(100, 30));
		G_SET_NEW_BUTTON(delete_button, "-", QSize(30, 30));
		G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(column_layout, row_layout, line_edit, label, choose_color_button, delete_button);
		connect(delete_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_delete_palette);
		connect(choose_color_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_choose_color);

		this->row_layouts_ << row_layout;
		this->factor_line_edits_ << line_edit;
		this->color_view_labels_ << label;
		this->choose_color_buttons_ << choose_color_button;
		this->delete_buttons_ << delete_button;

		column_layout->addLayout(this->add_layout_);

		++this->last_column_row_;
		++this->rows_per_column_[this->ncol_ - 1];
	}
	else {
		QVBoxLayout* column_layout = this->column_layouts_.back();
		column_layout->removeItem(this->add_layout_);

		G_SET_NEW_LINEEDIT(line_edit, "", soap::MiddleSize);
		G_SET_NEW_LABEL(label, "", QSize(80, 30));
		G_SET_NEW_BUTTON(choose_color_button, "Choose Color", QSize(100, 30));
		G_SET_NEW_BUTTON(delete_button, "-", QSize(30, 30));
		G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(column_layout, row_layout, line_edit, label, choose_color_button, delete_button);
		connect(delete_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_delete_palette);
		connect(choose_color_button, &QPushButton::clicked, this, &PaletteSettingDialog::s_choose_color);

		this->row_layouts_ << row_layout;
		this->factor_line_edits_ << line_edit;
		this->color_view_labels_ << label;
		this->choose_color_buttons_ << choose_color_button;
		this->delete_buttons_ << delete_button;

		column_layout = new QVBoxLayout;
		this->column_layouts_ << column_layout;
		column_layout->addLayout(this->add_layout_);
		this->main_layout_->addLayout(column_layout);

		++this->ncol_;
		this->last_column_row_ = 1;
		this->rows_per_column_ << 1;
	}
}

void PaletteSettingDialog::s_delete_palette() {

	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());

	qsizetype index = this->delete_buttons_.indexOf(button);

	QHBoxLayout* row_layout = this->row_layouts_[index];

	QLayoutItem* child;
	while ((child = row_layout->takeAt(0)) != 0) {
		if (child->widget()) {
			child->widget()->setParent(nullptr);
		}
		delete child;
	}

	delete row_layout;

	this->row_layouts_.remove(index);
	this->factor_line_edits_.remove(index);
	this->color_view_labels_.remove(index);
	this->choose_color_buttons_.remove(index);
	this->delete_buttons_.remove(index);
};

void PaletteSettingDialog::accept() {
	QDialog::accept();
	this->is_accepted_ = true;
}

void PaletteSettingDialog::reject() {
	QDialog::reject();
}