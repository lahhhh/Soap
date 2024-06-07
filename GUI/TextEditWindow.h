#pragma once

#include "Identifier.h"

#include <QMainWindow>
#include <QHBoxLayout>
#include <QTextEdit>

#include "SignalEmitter.h"

class TextEditWindow :
	public QMainWindow
{
public:
	TextEditWindow(
		QString* note,
		void* from,
		SignalEmitter* signal_emitter);

	TextEditWindow(
		QStringList* string_vector_data,
		void* from,
		SignalEmitter* signal_emitter);


	enum class WorkMode {
		Note,
		StringVectorEdit
	};

	void set_layout();

	void set_property();

	WorkMode mode_{ WorkMode::Note };

	QString* note_{ nullptr };
	QStringList* string_vector_data_{ nullptr };

	void* from_{ nullptr };

	SignalEmitter* signal_emitter_{ nullptr };

	QHBoxLayout* main_layout_{ nullptr };

	QWidget* main_interface_{ nullptr };

	QTextEdit* text_edit_{ nullptr };

	bool edited_{ false };
	bool valid_{ true };

	static void view(
		QString* note,
		void* from,
		SignalEmitter* signal_emitter);

	static void view(
		QStringList* string_vector_data,
		void* from,
		SignalEmitter* signal_emitter);

private slots:

	void s_check_data(void* data, soap::VariableType type, void* item);

	void s_edit();

	void closeEvent(QCloseEvent* e);
};

