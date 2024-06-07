#pragma once
#include "qcustomplot.h"

#include "PlotUtility.h"
class SoapTextElement :
    public QCPTextElement
{
    Q_OBJECT

public:
    explicit SoapTextElement(QCustomPlot* parent_plot): QCPTextElement(parent_plot) {};

    SoapTextElement(QCustomPlot* parent_plot, const QString& text) : QCPTextElement(parent_plot, text) {};

    SoapTextElement(
        QCustomPlot* parent_plot, 
        const QString& text, 
        double pointSize
    ) :
        QCPTextElement(
            parent_plot, 
            text, 
            pointSize
        ) {};

    SoapTextElement(
        QCustomPlot* parent_plot, 
        const QString& text, 
        const QString& fontFamily, 
        double pointSize
    ) :
        QCPTextElement(
            parent_plot, 
            text, 
            fontFamily, 
            pointSize
        ) {};

    SoapTextElement(
        QCustomPlot* parent_plot,
        const QString& text,
        const QFont& font
    ) : 
        QCPTextElement(
            parent_plot,
            text, 
            font
        ){};

    virtual QSize minimumOuterSizeHint() const Q_DECL_OVERRIDE {
        QTextDocument td;
        td.setDefaultFont(mFont);
        td.setHtml(mText);
        QSize result(td.size().toSize());
        result.rwidth() += mMargins.left() + mMargins.right();
        result.rheight() += mMargins.top() + mMargins.bottom();
        return result;
    };

    virtual QSize maximumOuterSizeHint() const Q_DECL_OVERRIDE {
        QTextDocument td;
        td.setDefaultFont(mFont);
        td.setHtml(mText);
        QSize result(td.size().toSize());
        result.setWidth(QWIDGETSIZE_MAX);
        result.rheight() += mMargins.top() + mMargins.bottom();
        return result;
    };

protected:

    virtual void draw(QCPPainter* painter) Q_DECL_OVERRIDE {

        auto text_font = mainFont();

        QTextDocument td;
        
        td.setDefaultFont(mainFont());
        td.setHtml(mText);        

        auto [width, height] = td.size().toSize();
        painter->translate(mRect.x() + mRect.width() / 2 - width / 2, mRect.y() + mRect.height() / 2 - height / 2);

        td.drawContents(painter);
    };
};

