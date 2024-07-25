#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QScrollArea>

#include "SignalEmitter.h"

#include "Switch.h"
#include "SimpleChoiceLayout.h"
#include "SimpleChoiceWithLineEditLayout.h"
#include "SimpleChoiceWithLineEditAndColorChoiceLayout.h"
#include "FactorChoiceLayout.h"
#include "FactorChoiceWithLineEditLayout.h"
#include "FactorDoubleLineEditWithCompleterLayout.h"
#include "MultipleDoubleLineEditWithCompleterLayout.h"
#include "MultipleLineEditLayout.h"
#include "MultipleLineEditWithCompleterLayout.h"
#include "ChooseMultiFileLayout.h"
#include "MultiLineEditWithColorChoiceLayout.h"
#include "MultiCheckBoxLayout.h"
#include "QTextEdit.h"
#include "CompareLayout.h"
#include "ChooseColorLayout.h"
#include "ChooseOpenFileLayout.h"
#include "ChooseSaveFileLayout.h"
#include "LogicLayout.h"

#include "LogicHandler.h"

class CommonDialog 
	: public QDialog
{
	Q_OBJECT
public:
	CommonDialog(
		SignalEmitter* signal_emitter,
		const QString& title, 
		const QStringList& labels, 
		const QList<soap::InputStyle>& input_style,
		const QList<QStringList>& qstring_list_items, 
		const QList<QMap<QString, QStringList>>& qstring_map_items,
		const QList<void*>& ptrs
	);

	static QStringList get_response(
		SignalEmitter* signal_emitter,
		const QString& title, 
		const QStringList& labels, 
		const QList<soap::InputStyle>& input_style,
		const QList<QStringList>& qstring_list_items = {},
		const QList<QMap<QString, QStringList>>& qstring_map_items = {},
		const QList<void*>& ptrs = {}
	);

private:

	SignalEmitter* signal_emitter_{ nullptr };

	QVBoxLayout* all_layout_{ nullptr };
	QHBoxLayout* main_layout_{ nullptr };

	QList<QLineEdit*> line_edits_integer_;
	QList<QLineEdit*> line_edits_numeric_;
	QList<QLineEdit*> line_edits_string_;
	QList<QLineEdit*> line_edits_with_completer_;
	QList<QLineEdit*> line_edits_with_completer_2_;
	QList<QComboBox*> combo_boxes_;
	QList<Switch*> switches_;
	QList<SimpleChoiceLayout*> simple_choice_layouts_;
	QList<FactorChoiceLayout*> factor_choice_layouts_;
	QList<MultipleLineEditLayout*> multiple_line_edit_layouts_;
	QList<MultipleLineEditWithCompleterLayout*> multiple_line_edit_with_completer_layouts_;
	QList<MultipleDoubleLineEditWithCompleterLayout*> multiple_double_line_edit_with_completer_layouts_;
	QList<SimpleChoiceWithLineEditLayout*> simple_choice_with_line_edit_layouts_;
	QList<SimpleChoiceWithLineEditAndColorChoiceLayout*> simple_choice_with_line_edit_and_color_choice_layouts_;
	QList<FactorChoiceWithLineEditLayout*> factor_choice_with_line_edit_layouts_;
	QList<FactorDoubleLineEditWithCompleterLayout*> factor_double_line_edit_with_completer_layouts_;
	QList<ChooseMultiFileLayout*> choose_multiple_file_layouts_;
	QList<MultiLineEditWithColorChoiceLayout*> multi_line_edit_with_color_choice_layouts_;
	QList<MultiCheckBoxLayout*> multi_checkbox_layouts_;
	QList<QTextEdit*> text_edits_;
	QList<CompareLayout*> compare_layouts_;
	QList<ChooseColorLayout*> choose_color_layouts_;
	QList<ChooseOpenFileLayout*> choose_open_file_layouts_;
	QList<ChooseSaveFileLayout*> choose_save_file_layouts_;
	QList<LogicLayout*> logic_layouts_;

	QPushButton* finish_button_{ nullptr };
	QPushButton* cancel_button_{ nullptr };

	bool is_accepted_ = false;

private slots:

	void accept();	
	void reject();
};

QStringList compare_layouts_to_list(const QString& cl);

bool switch_to_bool(const QString& sw);

QStringList multiple_file_to_list(const QString& mf);

QStringList multiple_line_edit_to_list(const QString& al);

QStringList multiple_line_edit_with_completer_to_list(const QString& al);

QPair<QStringList, QStringList> multiple_double_line_edit_with_completer_layout_to_pair(const QString& mdl);

QStringList simple_choice_to_list(const QString& sc);

QPair<QString, QStringList> factor_choice_to_pair(const QString& fc);

QPair<QPair<QString, QStringList>, QPair<QString, QStringList>> factor_choice_with_line_edit_to_pair(const QString& fcw);

QPair<QPair<QString, QStringList>, QPair<QString, QStringList>> factor_double_line_edit_with_completer_to_pair(const QString& fcw);

QPair<QStringList, QStringList> simple_choice_with_line_edit_layout_to_pair(const QString& scw);

std::tuple<QStringList, QStringList, QList<QColor>> simple_choice_with_line_edit_and_color_choice_layout_to_tuple(const QString& slc);

std::pair<QStringList, QList<QColor>> multi_line_edit_with_color_choice_layout_to_pair(const QString& mlc);

QStringList multi_check_box_to_list(const QString& mcb);