#include "NoteItem.h"

#include "TextEditWindow.h"

void NoteItem::__show_this() {

	TextEditWindow::view(static_cast<QString*>(this->data_), this->from_, this->signal_emitter_);
};