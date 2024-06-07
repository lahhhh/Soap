#pragma once

#include "VariableItem.h"

class NoteItem : 
    public VariableItem
{

public:

    NoteItem(QString* data, void* from, SignalEmitter* signal_emitter) :
        from_(from)
    {
        this->data_ = data;

        this->setText(0, "Note");

        this->signal_emitter_ = signal_emitter;
    };

    QString* data() {
        return static_cast<QString*>(this->data_);
    }

    void* from_{ nullptr };

    void __show_this() override;
};

