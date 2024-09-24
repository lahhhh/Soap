#include "CellChatItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"

void CellChatItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

