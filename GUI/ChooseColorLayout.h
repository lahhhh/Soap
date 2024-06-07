#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>

class ChooseColorLayout
	: public QWidget
{
public:

	ChooseColorLayout(
		const QString& default_color_name,
		QWidget* parent
	);

	QPushButton* button_{ nullptr };

	QString color_name_;

	void s_choose_color();

	QString current_value();
};

