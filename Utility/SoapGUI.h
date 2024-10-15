#pragma once

#include "CustomPlot.h"

#include <QComboBox>

#define G_SINGLE_ITEM_ROWLAYOUT(layout,item) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item); \
layout->addStretch(); 

#define G_NEW_SINGLE_ITEM_ROWLAYOUT(layout,item) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item); \
layout->addStretch(); 

#define G_ADD_SINGLE_ITEM_ROWLAYOUT(main_layout, layout,item) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT(main_layout, layout,item) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout,item) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item); \
main_layout->addLayout(layout);

#define G_SINGLE_ITEM_ROWLAYOUT_NO_STRETCH(layout,item) \
layout = new QHBoxLayout; \
layout->addWidget(item); 

#define G_NEW_SINGLE_ITEM_ROWLAYOUT_NO_STRETCH(layout,item) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item); 

#define G_ADD_SINGLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout,item) \
layout = new QHBoxLayout; \
layout->addWidget(item); \
main_layout->addLayout(layout);

#define G_ADD_NEW_SINGLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout,item) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item); \
main_layout->addLayout(layout);

#define G_DOUBLE_ITEM_ROWLAYOUT(layout,item1,item2) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch();

#define G_NEW_DOUBLE_ITEM_ROWLAYOUT(layout,item1,item2) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch();

#define G_ADD_DOUBLE_ITEM_ROWLAYOUT(main_layout, layout,item1,item2) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch();\
main_layout->addLayout(layout);

#define G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT(main_layout, layout,item1,item2) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch();\
main_layout->addLayout(layout);

#define G_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(layout,item1,item2) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); 

#define G_LEFT_ALIGN_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(layout,item1,item2) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addStretch();

#define G_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(layout,item1,item2) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); 

#define G_LEFT_ALIGN_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(layout,item1,item2) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addStretch();

#define G_ADD_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout,item1,item2) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
main_layout->addLayout(layout);

#define G_ADD_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout,item1,item2) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout,item1,item2) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
main_layout->addLayout(layout);

#define G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout,item1,item2) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_TRIPLE_ITEM_ROWLAYOUT(layout, item1, item2, item3) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch();

#define G_NEW_TRIPLE_ITEM_ROWLAYOUT(layout, item1, item2, item3) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch();

#define G_ADD_TRIPLE_ITEM_ROWLAYOUT(main_layout, layout,item1,item2, item3) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch();\
main_layout->addLayout(layout);

#define G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT(main_layout, layout,item1,item2, item3) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch();\
main_layout->addLayout(layout);

#define G_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(layout, item1, item2, item3) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3);

#define G_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(layout, item1, item2, item3) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3);

#define G_ADD_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout, item1, item2, item3) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
main_layout->addLayout(layout);

#define G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout, item1, item2, item3) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
main_layout->addLayout(layout);

#define G_ADD_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout, item1, item2, item3) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout, item1, item2, item3) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_TRIPLE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout, item1, item2, item3) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
main_layout->addLayout(layout);

#define G_FOUR_ITEM_ROWLAYOUT(layout, item1, item2, item3, item4) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch(); \
layout->addWidget(item4); \
layout->addStretch();

#define G_NEW_FOUR_ITEM_ROWLAYOUT(layout, item1, item2, item3, item4) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch(); \
layout->addWidget(item4); \
layout->addStretch();

#define G_ADD_FOUR_ITEM_ROWLAYOUT(main_layout, layout,item1,item2, item3, item4) \
layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch(); \
layout->addWidget(item4); \
layout->addStretch();\
main_layout->addLayout(layout);

#define G_ADD_NEW_FOUR_ITEM_ROWLAYOUT(main_layout, layout,item1,item2, item3, item4) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addStretch(); \
layout->addWidget(item1); \
layout->addStretch(); \
layout->addWidget(item2); \
layout->addStretch(); \
layout->addWidget(item3); \
layout->addStretch(); \
layout->addWidget(item4); \
layout->addStretch();\
main_layout->addLayout(layout);

#define G_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(layout, item1, item2, item3, item4) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4);

#define G_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(layout, item1, item2, item3, item4) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4);

#define G_ADD_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout, item1, item2, item3, item4) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
main_layout->addLayout(layout);

#define G_ADD_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout, item1, item2, item3, item4) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout, item1, item2, item3, item4) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_FOUR_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout, item1, item2, item3, item4) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
main_layout->addLayout(layout);

#define G_ADD_FIVE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout, item1, item2, item3, item4, item5) \
layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
layout->addWidget(item5); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_ADD_NEW_FIVE_ITEM_ROWLAYOUT_NO_STRETCH(main_layout, layout, item1, item2, item3, item4, item5) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
layout->addWidget(item5); \
main_layout->addLayout(layout);

#define G_ADD_NEW_FIVE_ITEM_ROWLAYOUT_LEFT_ALIGN_NO_STRETCH(main_layout, layout, item1, item2, item3, item4, item5) \
QHBoxLayout* layout = new QHBoxLayout; \
layout->addWidget(item1); \
layout->addWidget(item2); \
layout->addWidget(item3); \
layout->addWidget(item4); \
layout->addWidget(item5); \
layout->addStretch(); \
main_layout->addLayout(layout);

#define G_SET_ICON setWindowIcon(QIcon(FILE_SOAP_ICON_JPG));

#define G_SET_BUTTON_FINISH_STYLE(button,text) \
button = new QPushButton(text, this);\
button->setFixedSize(120, 30);\
button->setFont(QFont("Arial", 10, 700));\
button->setStyleSheet("QPushButton{color:#000; background-color:#e8ffe8; border-radius:15px;}QPushButton:hover{background-color:#d8ffd8;}");

#define G_SET_BUTTON_CANCEL_STYLE(button,text) \
button = new QPushButton(text, this);\
button->setFixedSize(120, 30);\
button->setFont(QFont("Arial", 10, 700));\
button->setStyleSheet("QPushButton{color:#000; background-color:#ffe8e8; border-radius:15px;}QPushButton:hover{background-color:#ffd8d8;}");


#define G_SET_CANCEL_BUTTON \
this->cancel_button_ = new QPushButton("Cancel", this);\
this->cancel_button_->setFixedSize(120, 30);\
this->cancel_button_->setFont(QFont("Arial", 10, 700));\
this->cancel_button_->setStyleSheet("QPushButton{color:#000; background-color:#ffe8e8; border-radius:15px;}QPushButton:hover{background-color:#ffd8d8;}");

#define G_SET_FINISH_BUTTON \
this->finish_button_ = new QPushButton("Finish", this);\
this->finish_button_->setFixedSize(120, 30);\
this->finish_button_->setFont(QFont("Arial", 10, 700));\
this->finish_button_->setStyleSheet("QPushButton{color:#000; background-color:#e8ffe8; border-radius:15px;}QPushButton:hover{background-color:#d8ffd8;}");


#define G_SET_PLOTSUITE_BUTTON \
this->previous_picture_button_ = new QPushButton("◄", this);\
this->previous_picture_button_->setFixedSize(120, 30);\
this->previous_picture_button_->setFont(QFont("Arial", 10, 700));\
this->previous_picture_button_->setStyleSheet("QPushButton{color:#000; background-color:#e8e8ff; border-radius:15px;}QPushButton:hover{background-color:#d8d8ff;}");\
\
this->next_picture_button_ = new QPushButton("►", this);\
this->next_picture_button_->setFixedSize(120, 30);\
this->next_picture_button_->setFont(QFont("Arial", 10, 700));\
this->next_picture_button_->setStyleSheet("QPushButton{color:#000; background-color:#e8e8ff; border-radius:15px;}QPushButton:hover{background-color:#d8d8ff;}");\
\
this->clear_picture_button_ = new QPushButton("Clear", this);\
this->clear_picture_button_->setFixedSize(120, 30);\
this->clear_picture_button_->setFont(QFont("Arial", 10, 700));\
this->clear_picture_button_->setStyleSheet("QPushButton{color:#000; background-color:#e8e8ff; border-radius:15px;}QPushButton:hover{background-color:#d8d8ff;}");\
\
this->pop_picture_button_ = new QPushButton("Pop", this);\
this->pop_picture_button_->setFixedSize(120, 30);\
this->pop_picture_button_->setFont(QFont("Arial", 10, 700));\
this->pop_picture_button_->setStyleSheet("QPushButton{color:#000; background-color:#e8e8ff; border-radius:15px;}QPushButton:hover{background-color:#d8d8ff;}");\
\
this->save_picture_button_ = new QPushButton("Save Picture", this);\
this->save_picture_button_->setFixedSize(120, 30);\
this->save_picture_button_->setFont(QFont("Arial", 10, 700));\
this->save_picture_button_->setStyleSheet("QPushButton{color:#000; background-color:#e8e8ff; border-radius:15px;}QPushButton:hover{background-color:#d8d8ff;}");\

#define G_SET_BUTTON(button, text, size) \
button = new QPushButton(text, this);\
button->setFixedSize(size);

#define G_SET_BUTTON_ICON(button, icon_path, size) \
button = new QPushButton(this);\
button->setIcon(QIcon(icon_path)); \
button->setFixedSize(size);

#define G_SET_NEW_BUTTON_ICON(button, icon_path, size) \
QPushButton* button = new QPushButton(this);\
button->setIcon(QIcon(icon_path)); \
button->setFixedSize(size);

#define G_SET_NEW_BUTTON(button, text, size) \
QPushButton* button = new QPushButton(text, this);\
button->setFixedSize(size);

#define G_SET_LABEL(label, text, size) \
label = new QLabel(text, this);\
label->setFixedSize(size);

#define G_SET_EMPTY_LABEL(label, size) \
label = new QLabel(this);\
label->setFixedSize(size);

#define G_SET_LABEL_PRECISE(label, text, size, color, font)\
label = new QLabel(text, this);\
label->setFixedSize(size);\
label->setStyleSheet("QLabel{color: "+ color.name() + "}");\
label->setFont(font);

#define G_SET_NEW_LABEL(label, text, size) \
QLabel* label = new QLabel(text, this);\
label->setFixedSize(size);

#define G_SET_NEW_LABEL_ADJUST_SIZE(label, text) \
QLabel* label = new QLabel(text, this);\
label->adjustSize();

#define G_SET_LABEL_FIXED_HEIGHT(label, text, height) \
label = new QLabel(text, this);\
label->adjustSize();\
label->setFixedHeight(height);

#define G_SET_NEW_LABEL_FIXED_HEIGHT(label, text, height) \
QLabel* label = new QLabel(text, this);\
label->adjustSize();\
label->setFixedHeight(height);

#define G_SET_LINEEDIT(line_edit, text, size) \
line_edit = new QLineEdit(text, this);\
line_edit->setFixedSize(size);

#define G_SET_EMPTY_LINEEDIT(line_edit, size) \
line_edit = new QLineEdit(this);\
line_edit->setFixedSize(size);

#define G_SET_LINEEDIT_WITH_COMPLETER(line_edit, text, completer_content, size) \
line_edit = new QLineEdit(text, this);\
{\
	QCompleter* __completer = new QCompleter(completer_content, this);\
	line_edit->setCompleter(__completer);\
}\
line_edit->setFixedSize(size);

#define G_SET_NEW_LINEEDIT_WITH_COMPLETER(line_edit, text, completer_content, height)\
QLineEdit* line_edit = new QLineEdit(text, this);\
{\
	QCompleter* __completer = new QCompleter(completer_content, this);\
	line_edit->setCompleter(__completer);\
}\
line_edit->adjustSize();\
line_edit->setFixedHeight(height);

#define G_SET_NEW_LINEEDIT(line_edit, text, size) \
QLineEdit* line_edit = new QLineEdit(text, this);\
line_edit->setFixedSize(size);

#define G_SET_NEW_COMBOBOX(combo_box, items, height)\
QComboBox* combo_box = new QComboBox(this); \
combo_box->addItems(items); \
{\
	int __combo_box_width = custom_plot::utility::get_max_text_width(items, QFont("Arial", 15)) + 20;\
	if (__combo_box_width < 150) {\
		__combo_box_width = 150;\
	}\
	combo_box->setFixedWidth(__combo_box_width);\
}\
combo_box->setFixedHeight(height);

#define G_SET_COMBOBOX(combo_box, items, height)\
combo_box = new QComboBox(this); \
combo_box->addItems(items); \
{\
	int __width = custom_plot::utility::get_max_text_width(items, QFont("Arial", 15)) + 20;\
	if (__width < 150) {\
		__width = 150;\
	}\
	combo_box->view()->setMinimumWidth(__width);\
}\
combo_box->setFixedWidth(150);\
combo_box->setFixedHeight(height);

#define G_SET_COMBOBOX_FIXED_WIDTH(combo_box, items, width, height)\
combo_box = new QComboBox(this); \
combo_box->addItems(items); \
{\
	int __width = custom_plot::utility::get_max_text_width(items, QFont("Arial", 15)) + 20;\
	if (__width < 150) {\
		__width = 150;\
	}\
	combo_box->view()->setMinimumWidth(__width);\
}\
combo_box->setFixedWidth(width);\
combo_box->setFixedHeight(height);

#define G_SET_SWITCH(sw, status, link_item, size)\
sw = new Switch(status, link_item, this);\
sw->setFixedSize(size);

#define G_SET_NEW_SWITCH(sw, status, link_item, size)\
Switch* sw = new Switch(status, link_item, this);\
sw->setFixedSize(size);

#define G_SET_SWITCH_NO_LINK(sw, status, size)\
sw = new Switch(status, this);\
sw->setFixedSize(size);

#define G_SET_LABEL_COLOR(label, color)\
label->setStyleSheet("QLabel{color: "+ color.name() + "}");

#define G_SET_LARGE_LABEL(label)\
label->setFont(QFont("Arial", 20, 700));