#pragma once

#include "Identifier.h"

#include <QMainWindow>
#include <QScrollArea>
#include <QPushButton>
#include <QTextBrowser>
#include <QHBoxLayout>

class PageWindow :
    public QMainWindow
{
public:
    PageWindow(
        const QString& title,
        const QString& file_name, 
        QWidget* parent = nullptr
    );

    QString file_name_;

    QWidget* main_interface_;

    QHBoxLayout* main_layout_;

    QScrollArea* scroll_area_;

    QWidget* scroll_area_widget_;

    QVBoxLayout* scroll_area_layout_;

    QTextBrowser* text_browser_;

    QList<QPushButton*> questions_;

    QStringList answers_;

    int scroll_area_width_{ 600 };

    void get_content();

    static void show_this(
        const QString& title,
        const QString& file_name,
        QWidget* parent = nullptr
    );

private slots:

    void s_question_asked();

};

