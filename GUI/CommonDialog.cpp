#include "CommonDialog.h"

#include <QGridLayout>

#include "SoapGUI.h"

// ban mouse wheel
void QComboBox::wheelEvent(QWheelEvent* e){}

CommonDialog::CommonDialog(
	SignalEmitter* signal_emitter,
	const QString& title,
	const QStringList& labels,
	const QList<soap::InputStyle>& input_style,
	const QList<QStringList>& qstring_list_items,
	const QList<QMap<QString, QStringList> >& qstring_map_items,
	const QList<void*>& ptrs
)
	:
	signal_emitter_(signal_emitter)
{
	this->main_layout_ = new QHBoxLayout;
	this->all_layout_ = new QVBoxLayout;

	int list_count{ 0 }, map_count{ 0 }, ptr_count{ 0 };
	const int size = labels.size();
	int ncol = size / 10 + 1;
	const int row_per_column = ceil((double)size / ncol);

	ncol = ceil((double)size / row_per_column);

	QVector<int> rows_per_column(ncol, row_per_column);
	rows_per_column[ncol - 1] = size - (ncol - 1) * row_per_column;

	for (int col = 0; col < ncol; col++) {

		QGridLayout* grid_layout = new QGridLayout;

		for (int i = 0; i < rows_per_column[col]; ++i)
		{
			int index = col * row_per_column + i;

			QString label, defaults;
			QStringList gui_list = labels[index].split(':');

			if (gui_list.size() > 1) {
				defaults = gui_list.last();
				label = gui_list.sliced(0, gui_list.size() - 1).join("");
			}
			else {
				label = gui_list.front();
			}

			G_SET_NEW_LABEL_FIXED_HEIGHT(row_label, label, 30);

			grid_layout->addWidget(row_label, i, 0);

			soap::InputStyle this_input_style = input_style[index];

			if (this_input_style == soap::InputStyle::IntegerLineEdit) {

				G_SET_NEW_LINEEDIT(line_edit, defaults, soap::MiddleSize);

				line_edit->setValidator(new QIntValidator());

				grid_layout->addWidget(line_edit, i, 1);

				this->line_edits_integer_ << line_edit;
			}
			else if (this_input_style == soap::InputStyle::NumericLineEdit) {
				G_SET_NEW_LINEEDIT(line_edit, defaults, soap::MiddleSize);

				line_edit->setValidator(new QDoubleValidator());

				grid_layout->addWidget(line_edit, i, 1);

				this->line_edits_numeric_ << line_edit;
			}
			else if (this_input_style == soap::InputStyle::StringLineEdit) {

				G_SET_NEW_LINEEDIT(line_edit, defaults, soap::MiddleSize);

				grid_layout->addWidget(line_edit, i, 1);

				this->line_edits_string_ << line_edit;
			}
			else if (this_input_style == soap::InputStyle::ComboBox) {

				QStringList items = qstring_list_items[list_count++];

				G_SET_NEW_COMBOBOX(combo_box, items, 30);
				grid_layout->addWidget(combo_box, i, 1);
				this->combo_boxes_ << combo_box;
			}
			else if (this_input_style == soap::InputStyle::LineEditWithCompleter) {

				G_SET_NEW_LINEEDIT_WITH_COMPLETER(line_edit, defaults, qstring_list_items[list_count++], 30);

				grid_layout->addWidget(line_edit, i, 1);

				this->line_edits_with_completer_ << line_edit;
			}
			else if (this_input_style == soap::InputStyle::LineEditWithCompleter2) {

				auto ele = new QLineEdit(defaults, this); 
				
				auto completer = new QCompleter(qstring_list_items[list_count++], this);
				completer->setFilterMode(Qt::MatchContains);

				ele->setCompleter(completer);
				ele->adjustSize();
				ele->setFixedHeight(30);

				grid_layout->addWidget(ele, i, 1);

				this->line_edits_with_completer_2_ << ele;
			}
			else if (this_input_style == soap::InputStyle::SwitchButton) {

				bool initial_value = defaults == "yes";

				G_SET_NEW_SWITCH(sw, initial_value, row_label, soap::MiddleSize);

				grid_layout->addWidget(sw, i, 1);

				this->switches_ << sw;
			}
			else if (this_input_style == soap::InputStyle::SimpleChoice) {
				SimpleChoiceLayout* scl = new SimpleChoiceLayout(label, qstring_list_items[list_count++], this);
				grid_layout->addWidget(scl, i, 1);
				this->simple_choice_layouts_ << scl;
			}
			else if (this_input_style == soap::InputStyle::FactorChoice) {
				FactorChoiceLayout* fcl = new FactorChoiceLayout(defaults, qstring_map_items[map_count++], this);
				grid_layout->addWidget(fcl, i, 1);
				this->factor_choice_layouts_ << fcl;
			}
			else if (this_input_style == soap::InputStyle::MultipleLineEdit) {
				MultipleLineEditLayout* all = new MultipleLineEditLayout(label, this->signal_emitter_, {}, this);
				grid_layout->addWidget(all, i, 1);
				this->multiple_line_edit_layouts_ << all;
			}
			else if (this_input_style == soap::InputStyle::MultipleLineEditWithCompleter) {
				MultipleLineEditWithCompleterLayout* mew = new MultipleLineEditWithCompleterLayout(label, qstring_list_items[list_count++], this->signal_emitter_, this);
				grid_layout->addWidget(mew, i, 1);
				this->multiple_line_edit_with_completer_layouts_ << mew;
			}
			else if (this_input_style == soap::InputStyle::SimpleChoiceWithLineEdit) {
				SimpleChoiceWithLineEditLayout* scw = new SimpleChoiceWithLineEditLayout(label, defaults, qstring_list_items[list_count++], this);
				grid_layout->addWidget(scw, i, 1);
				this->simple_choice_with_line_edit_layouts_ << scw;
			}
			else if (this_input_style == soap::InputStyle::SimpleChoiceWithLineEditAndColorChoiceLayout) {
				SimpleChoiceWithLineEditAndColorChoiceLayout* swc = new SimpleChoiceWithLineEditAndColorChoiceLayout(label, defaults, qstring_list_items[list_count++], this);
				grid_layout->addWidget(swc, i, 1);
				this->simple_choice_with_line_edit_and_color_choice_layouts_ << swc;
			}
			else if (this_input_style == soap::InputStyle::FactorChoiceWithLineEdit) {
				FactorChoiceWithLineEditLayout* fcw = new FactorChoiceWithLineEditLayout(label, defaults, qstring_map_items[map_count++], this);
				grid_layout->addWidget(fcw, i, 1);
				this->factor_choice_with_line_edit_layouts_ << fcw;
			}
			else if (this_input_style == soap::InputStyle::MultipleDoubleLineEditWithCompleterLayout) {
				MultipleDoubleLineEditWithCompleterLayout* mdl = new MultipleDoubleLineEditWithCompleterLayout(label, defaults, qstring_list_items[list_count++], this);
				grid_layout->addWidget(mdl, i, 1);
				this->multiple_double_line_edit_with_completer_layouts_ << mdl;
			}
			else if (this_input_style == soap::InputStyle::FactorDoubleLineEditWithCompleterLayout) {
				FactorDoubleLineEditWithCompleterLayout* fdl = new FactorDoubleLineEditWithCompleterLayout(label, defaults, qstring_map_items[map_count++], this);
				grid_layout->addWidget(fdl, i, 1);
				this->factor_double_line_edit_with_completer_layouts_ << fdl;
			}
			else if (this_input_style == soap::InputStyle::MultiFile) {
				ChooseMultiFileLayout* cmf = new ChooseMultiFileLayout(defaults, qstring_list_items[list_count++], this);
				grid_layout->addWidget(cmf, i, 1);
				this->choose_multiple_file_layouts_ << cmf;
			}
			else if (this_input_style == soap::InputStyle::MultiLineEditWithColorChoiceLayout) {
				MultiLineEditWithColorChoiceLayout* mlc = new MultiLineEditWithColorChoiceLayout(defaults, this);
				grid_layout->addWidget(mlc, i, 1);
				this->multi_line_edit_with_color_choice_layouts_ << mlc;
			}
			else if (this_input_style == soap::InputStyle::MultiCheckBox) {
				MultiCheckBoxLayout* mcb = new MultiCheckBoxLayout(qstring_list_items[list_count++], this);
				grid_layout->addWidget(mcb, i, 1);
				this->multi_checkbox_layouts_ << mcb;
			}
			else if (this_input_style == soap::InputStyle::TextEdit) {
				QTextEdit* te = new QTextEdit(defaults, this);
				grid_layout->addWidget(te, i, 1);
				this->text_edits_ << te;
			}
			else if (this_input_style == soap::InputStyle::CompareLayout) {

				int type{ 0 };
				if (!defaults.isEmpty()) {
					type = defaults.toInt();
				}

				CompareLayout* cl = new CompareLayout(
					label,
					qstring_map_items[map_count++],
					type,
					this
				);

				grid_layout->addWidget(cl, i, 1);
				this->compare_layouts_ << cl;
			}
			else if (this_input_style == soap::InputStyle::ColorChoice) {
				ChooseColorLayout* ccl = new ChooseColorLayout(
					defaults, this
				);

				grid_layout->addWidget(ccl, i, 1);
				this->choose_color_layouts_ << ccl;
			}
			else if (this_input_style == soap::InputStyle::ChooseOpenFile) {
				ChooseOpenFileLayout* cfl = new ChooseOpenFileLayout(
					defaults, this
				);

				grid_layout->addWidget(cfl, i, 1);
				this->choose_open_file_layouts_ << cfl;
			}
			else if (this_input_style == soap::InputStyle::ChooseSaveFile) {
				auto* cfl = new ChooseSaveFileLayout(
					defaults, this
				);

				grid_layout->addWidget(cfl, i, 1);
				this->choose_save_file_layouts_ << cfl;
			}
			else if (this_input_style == soap::InputStyle::LogicLayout) {
				LogicLayout* ll = new LogicLayout(static_cast<LogicHandler*>(ptrs[ptr_count++]), this);
				grid_layout->addWidget(ll, i, 1);
				this->logic_layouts_ << ll;
			}
		}
		this->main_layout_->addLayout(grid_layout);
	}
	this->main_layout_->setSizeConstraint(QLayout::SetFixedSize);	

	QWidget* main_widget = new QWidget(this);
	main_widget->setLayout(this->main_layout_);
	main_widget->setObjectName("MainWidget");
	main_widget->setStyleSheet("#MainWidget{background-color:#f4f4ff;}");

	QScrollArea* scroll_area = new QScrollArea(this);
	scroll_area->setWidget(main_widget);
	scroll_area->setObjectName("ScrollArea");
	scroll_area->setStyleSheet("#ScrollArea{background-color:#f4f4ff;}");

	this->all_layout_->addWidget(scroll_area);

	G_SET_FINISH_BUTTON;
	G_SET_CANCEL_BUTTON;

	G_ADD_NEW_DOUBLE_ITEM_ROWLAYOUT(this->all_layout_, row_layout, this->finish_button_, this->cancel_button_);

	this->setLayout(this->all_layout_);

	connect(this->finish_button_, &QPushButton::clicked, this, &CommonDialog::accept);
	connect(this->cancel_button_, &QPushButton::clicked, this, &CommonDialog::reject);

	G_SET_ICON;
	this->setWindowTitle(title);
	this->exec();
}

QStringList CommonDialog::get_response(
	SignalEmitter* signal_emitter,
	const QString& title,
	const QStringList& labels,
	const QList<soap::InputStyle>& input_style,
	const QList<QStringList>& qstring_list_items,
	const QList<QMap<QString, QStringList> >& qstring_map_items,
	const QList<void*>& ptrs
) {
	CommonDialog dlg(
		signal_emitter,
		title,
		labels,
		input_style,
		qstring_list_items,
		qstring_map_items,
		ptrs
	);
	if (!dlg.is_accepted_)return QStringList();

	QStringList ret;

	QMap<soap::InputStyle, int> counts;

	const int size = labels.size();

	for (int i = 0; i < size; ++i) {

		auto this_input_style = input_style[i];

		if (this_input_style == soap::InputStyle::IntegerLineEdit) {
			ret << dlg.line_edits_integer_[counts[this_input_style]++]->text();
		}
		else if (this_input_style == soap::InputStyle::NumericLineEdit) {
			ret << dlg.line_edits_numeric_[counts[this_input_style]++]->text();
		}
		else if (this_input_style == soap::InputStyle::StringLineEdit) {
			ret << dlg.line_edits_string_[counts[this_input_style]++]->text();
		}
		else if (this_input_style == soap::InputStyle::LineEditWithCompleter) {
			ret << dlg.line_edits_with_completer_[counts[this_input_style]++]->text();
		}
		else if (this_input_style == soap::InputStyle::LineEditWithCompleter2) {
			ret << dlg.line_edits_with_completer_2_[counts[this_input_style]++]->text();
		}
		else if (this_input_style == soap::InputStyle::MultipleLineEdit) {
			ret << dlg.multiple_line_edit_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::MultipleLineEditWithCompleter) {
			ret << dlg.multiple_line_edit_with_completer_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::ComboBox) {
			ret << dlg.combo_boxes_[counts[this_input_style]++]->currentText();
		}
		else if (this_input_style == soap::InputStyle::SwitchButton) {
			ret << dlg.switches_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::SimpleChoice) {
			ret << dlg.simple_choice_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::FactorChoice) {
			ret << dlg.factor_choice_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::SimpleChoiceWithLineEdit) {
			ret << dlg.simple_choice_with_line_edit_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::SimpleChoiceWithLineEditAndColorChoiceLayout) {
			ret << dlg.simple_choice_with_line_edit_and_color_choice_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::FactorChoiceWithLineEdit) {
			ret << dlg.factor_choice_with_line_edit_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::MultipleDoubleLineEditWithCompleterLayout) {
			ret << dlg.multiple_double_line_edit_with_completer_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::FactorDoubleLineEditWithCompleterLayout) {
			ret << dlg.factor_double_line_edit_with_completer_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::MultiFile) {
			ret << dlg.choose_multiple_file_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::MultiLineEditWithColorChoiceLayout) {
			ret << dlg.multi_line_edit_with_color_choice_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::MultiCheckBox) {
			ret << dlg.multi_checkbox_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::TextEdit) {
			ret << dlg.text_edits_[counts[this_input_style]++]->toPlainText();
		}
		else if (this_input_style == soap::InputStyle::CompareLayout) {
			ret << dlg.compare_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::ColorChoice) {
			ret << dlg.choose_color_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::ChooseOpenFile) {
			ret << dlg.choose_open_file_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::ChooseSaveFile) {
			ret << dlg.choose_save_file_layouts_[counts[this_input_style]++]->current_value();
		}
		else if (this_input_style == soap::InputStyle::LogicLayout) {
			ret << dlg.logic_layouts_[counts[this_input_style]++]->current_value();
		}
	}

	return ret;
}

void CommonDialog::accept() {
	QDialog::accept();
	this->is_accepted_ = true;
}

void CommonDialog::reject() {
	QDialog::reject();
}

QStringList compare_layouts_to_list(const QString& cl) {
	return cl.split(SOAP_DELIMITER);
}

bool switch_to_bool(const QString& sw) {
	return sw == SWITCH_ACCEPT;
};

QStringList multiple_file_to_list(const QString& mf) {
	QStringList ret = mf.split(SOAP_DELIMITER);
	const qsizetype size = ret.size();

	if (size > 1) {
		return ret.sliced(1, size - 1);
	}
	else {
		return {};
	}
}

QStringList multiple_line_edit_to_list(const QString& al) {
	QStringList ret = al.split(SOAP_DELIMITER);
	const qsizetype size = ret.size();

	if (size > 1) {
		return ret.sliced(1, size - 1);
	}
	else {
		return {};
	}
}

QStringList multiple_line_edit_with_completer_to_list(const QString& al) {
	QStringList ret = al.split(SOAP_DELIMITER);
	const qsizetype size = ret.size();

	if (size > 1) {
		return ret.sliced(1, size - 1);
	}
	else {
		return {};
	}
}

QPair<QStringList, QStringList> multiple_double_line_edit_with_completer_layout_to_pair(const QString& mdl) {
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

QStringList simple_choice_to_list(const QString& sc) {
	QStringList ret = sc.split(SOAP_DELIMITER);

	const qsizetype size = ret.size();
	if (size > 1) {
		return ret.sliced(1, ret.size() - 1);
	}
	else {
		return {};
	}
}

QPair<QString, QStringList> factor_choice_to_pair(const QString& fc) {
	QStringList ret = fc.split(SOAP_DELIMITER);

	if (ret.size() == 1) {
		return { ret[0], {} };
	}
	else {
		return { ret[0], ret.sliced(1, ret.size() - 1) };
	}
}

QPair<QPair<QString, QStringList>, QPair<QString, QStringList>> factor_choice_with_line_edit_to_pair(const QString& fcw) {
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

QPair<QPair<QString, QStringList>, QPair<QString, QStringList>> factor_double_line_edit_with_completer_to_pair(const QString& fcw) {
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

QPair<QStringList, QStringList> simple_choice_with_line_edit_layout_to_pair(const QString& scw) {

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

std::tuple<QStringList, QStringList, QList<QColor>> simple_choice_with_line_edit_and_color_choice_layout_to_tuple(const QString& slc) {
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

std::pair<QStringList, QList<QColor>> multi_line_edit_with_color_choice_layout_to_pair(const QString& mlc) {

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

QStringList multi_check_box_to_list(const QString& mcb) {

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
