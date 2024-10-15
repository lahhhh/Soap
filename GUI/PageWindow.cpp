#include "PageWindow.h"

#include <QFile>

#include "Custom.h"
#include "Identifier.h"

#include "SoapGUI.h"

PageWindow::PageWindow(
	const QString& title, 
	const QString& file_name, 
	QWidget* parent
):
	QMainWindow(parent),
	file_name_(file_name)
{
	this->main_interface_ = new QWidget(this);
	this->setCentralWidget(this->main_interface_);

	this->main_layout_ = new QHBoxLayout;
	this->main_interface_->setLayout(this->main_layout_);

	this->scroll_area_ = new QScrollArea(this);
	this->scroll_area_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
	this->main_layout_->addWidget(this->scroll_area_);
	this->scroll_area_->setFixedWidth(this->scroll_area_width_);

	this->text_browser_ = new QTextBrowser(this);
	this->main_layout_->addWidget(this->text_browser_);
	this->text_browser_->setFont(QFont("Arial", 10));

	this->scroll_area_widget_ = new QWidget(this->scroll_area_);

	this->scroll_area_layout_ = new QVBoxLayout;
	this->scroll_area_widget_->setLayout(this->scroll_area_layout_);	

	this->get_content();

	this->scroll_area_->setWidget(this->scroll_area_widget_);

	this->setAttribute(Qt::WA_DeleteOnClose);
	G_SET_ICON;
	this->resize(1000, 800);
	this->setWindowTitle(title);
	this->show();
};

void PageWindow::get_content() {
	QFile file(this->file_name_);

	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return;

	QTextStream in(&file);

	QString question, answer;
	QFont question_font("Arial", 12);
	QFontMetrics fm(question_font);

	const int width = this->scroll_area_width_ - 10;

	while (true) {
		question = in.readLine();

		if (question.isNull()) {
			break;
		}

		answer = in.readLine();

		if (answer.isNull()) {
			break;
		}

		QPushButton* button = new QPushButton(custom::string_next_line(fm, question, width - 30), this->scroll_area_);
		button->setFixedWidth(width);
		button->setFont(question_font);
		button->adjustSize();

		this->questions_ << button;
		this->scroll_area_layout_->addWidget(button);
		this->answers_ << answer;
		connect(button, &QPushButton::clicked, this, &PageWindow::s_question_asked);
	}

};

void PageWindow::s_question_asked() {
	QPushButton* button = dynamic_cast<QPushButton*>(QObject::sender());
	int index = this->questions_.indexOf(button);
	this->text_browser_->setText(this->answers_[index]);
};

void PageWindow::show_this(
	const QString& title,
	const QString& file_name,
	QWidget* parent
) {
	PageWindow* fq = new PageWindow(title, file_name, parent);
};