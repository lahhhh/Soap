#pragma once

#include "Identifier.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

#include "SignalEmitter.h"

class MultipleLineEditWithCompleterLayout : public QWidget
{
    Q_OBJECT
public:
    MultipleLineEditWithCompleterLayout(
        const QString& name, 
        const QStringList& completer,
        SignalEmitter* signal_emitter,
        QWidget* parent = nullptr
    );

    QString name_;

    QStringList completer_;

    SignalEmitter* signal_emitter_{ nullptr };

    QVBoxLayout* main_layout_;
    QVBoxLayout* item_layout_;

    QPushButton* add_button_;
    QPushButton* import_button_{ nullptr };

    QList<QHBoxLayout*> row_layouts_;

    QList<QLineEdit*> line_edits_;

    QList<QPushButton*> buttons_;

    QString current_value();

private slots:
    void s_delete_line();

    void s_add_line();

    void s_import();
};
