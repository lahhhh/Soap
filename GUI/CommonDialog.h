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
		const QList<LogicHandler*> & logic_handlers
	);

	static QStringList get_response(
		SignalEmitter* signal_emitter,
		const QString& title, 
		const QStringList& labels, 
		const QList<soap::InputStyle>& input_style,
		const QList<QStringList>& qstring_list_items = {},
		const QList<QMap<QString, QStringList>>& qstring_map_items = {},
		const QList<LogicHandler*>& logic_handlers = {}
	);

private:

	QList<LogicHandler*> logic_handlers;

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

inline QStringList compare_layouts_to_list(const QString& cl) {
	return cl.split(SOAP_DELIMITER);
}

inline bool switch_to_bool(const QString& sw) {
	return sw == SWITCH_ACCEPT;
};

inline QStringList multiple_file_to_list(const QString& mf) {
	QStringList ret = mf.split(SOAP_DELIMITER);
	const qsizetype size = ret.size();

	if (size > 1) {
		return ret.sliced(1, size - 1);
	}
	else {
		return {};
	}
}

inline QStringList multiple_line_edit_to_list(const QString& al) {
	QStringList ret = al.split(SOAP_DELIMITER);
	const qsizetype size = ret.size();

	if (size > 1) {
		return ret.sliced(1, size - 1);
	}
	else {
		return {};
	}
}

inline QStringList multiple_line_edit_with_completer_to_list(const QString& al) {
	QStringList ret = al.split(SOAP_DELIMITER);
	const qsizetype size = ret.size();

	if (size > 1) {
		return ret.sliced(1, size - 1);
	}
	else {
		return {};
	}
}

inline QPair<QStringList, QStringList> multiple_double_line_edit_with_completer_layout_to_pair(const QString& mdl) {
	if (mdl.isEmpty()) {
		return QPair<QStringList, QStringList>();
	}

	QStringList tmp = mdl.split(SOAP_DELIMITER), ret1, ret2;
	const qsizetype size = (tmp.size() - 1) / 2;

	for (qsizetype i = 0; i < size; ++i) {
		ret1 << tmp[i * 2 + 1];
		ret2 << tmp[i * 2 + 2];
	}
	return qMakePair(ret1, ret2);
}

inline QStringList simple_choice_to_list(const QString& sc) {
	QStringList ret = sc.split(SOAP_DELIMITER);

	const qsizetype size = ret.size();
	if (size > 1) {
		return ret.sliced(1, ret.size() - 1);
	}
	else {
		return {};
	}
}

inline QPair<QString, QStringList> factor_choice_to_pair(const QString& fc) {
	QStringList ret = fc.split(SOAP_DELIMITER);

	if (ret.size() == 1) {
		return { ret[0], {} };
	}
	else {
		return { ret[0], ret.sliced(1, ret.size() - 1) };
	}
}

inline QPair<QPair<QString, QStringList>, QPair<QString, QStringList>> factor_choice_with_line_edit_to_pair(const QString& fcw) {
	QStringList ret = fcw.split(SOAP_DELIMITER);

	const qsizetype size = ret.size() / 2 - 1;
	QPair<QString, QStringList> ret1, ret2;
	ret1.first = ret[0];
	ret2.first = ret[1];

	for (qsizetype i = 0; i < size; ++i) {
		ret1.second << ret[i * 2 + 2];
		ret2.second << ret[i * 2 + 3];
	}
	return qMakePair(ret1, ret2);
}

inline QPair<QPair<QString, QStringList>, QPair<QString, QStringList>> factor_double_line_edit_with_completer_to_pair(const QString& fcw) {
	QStringList ret = fcw.split(SOAP_DELIMITER);

	const qsizetype size = ret.size() / 2 - 1;

	QPair<QString, QStringList> ret1, ret2;
	ret1.first = ret[0];
	ret2.first = ret[1];

	for (qsizetype i = 0; i < size; ++i) {
		ret1.second << ret[i * 2 + 2];
		ret2.second << ret[i * 2 + 3];
	}
	return qMakePair(ret1, ret2);
}

inline QPair<QStringList, QStringList> simple_choice_with_line_edit_layout_to_pair(const QString& scw) {
	
	QPair<QStringList, QStringList> ret;

	if (scw.isEmpty()) {
		return ret;
	}

	QStringList tmp = scw.split(SOAP_DELIMITER);

	const qsizetype size = (tmp.size() - 1) / 2;

	for (qsizetype i = 0; i < size; ++i) {
		ret.first << tmp[i * 2 + 1];
		ret.second << tmp[i * 2 + 2];
	}

	return ret;
}

inline std::tuple<QStringList, QStringList, QList<QColor>> simple_choice_with_line_edit_and_color_choice_layout_to_tuple(const QString& slc) {
	if (slc.isEmpty()) {
		return std::tuple<QStringList, QStringList, QList<QColor>>();
	}
	QStringList tmp = slc.split(SOAP_DELIMITER);

	const qsizetype size = (tmp.size() - 1) / 3;

	QStringList ret1, ret2;
	ret1.reserve(size);
	ret2.reserve(size);

	QList<QColor> ret3;
	ret3.reserve(size);

	for (qsizetype i = 0; i < size; ++i) {
		ret1 << tmp[i * 3 + 1];
		ret2 << tmp[i * 3 + 2];
		ret3 << QColor(tmp[i * 3 + 3]);
	}

	return std::make_tuple(ret1, ret2, ret3);
}

inline std::pair<QStringList, QList<QColor>> multi_line_edit_with_color_choice_layout_to_pair(const QString& mlc) {

	if (mlc.isEmpty()) {
		return std::pair<QStringList, QList<QColor>>();
	}
	QStringList tmp = mlc.split(SOAP_DELIMITER);

	const qsizetype size = (tmp.size() - 1) / 2;

	QStringList factors;
	factors.reserve(size);

	QList<QColor> colors;
	colors.reserve(size);
	
	for (qsizetype i = 0; i < size; ++i) {
		factors << tmp[i * 2 + 1];
		colors << QColor(tmp[i * 2 + 2]);
	}

	return { factors, colors };
}

inline QStringList multi_check_box_to_list(const QString& mcb) {

	if (mcb.isEmpty()) {
		return {};
	}

	QStringList tmp = mcb.split(SOAP_DELIMITER);

	if (tmp.size() > 1) {
		return tmp.sliced(1);
	}
	else {
		return {};
	}
}