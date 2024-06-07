#pragma once

#define COUNT_ARGS_HELPER(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, N, ...) N
#define COUNT_ARGS(...) COUNT_ARGS_HELPER(dummy, ##__VA_ARGS__, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define NARGS(...) COUNT_ARGS(__VA_ARGS__)

#define CONCAT_DIRECT(a, b) a ## b
#define CONCAT_INDIRECT(a, b) CONCAT_DIRECT(a, b)

#define G_CLASS_FUNCTION_DEFAULT(class_type) \
class_type() = default; \
class_type(const class_type&) = default; \
class_type(class_type&&) = default; \
class_type& operator=(const class_type&) = default; \
class_type& operator=(class_type&&) = default; \
~class_type() = default;

#define G_QUICK_ACCESS(type, name) \
type* name(){ \
	for (auto&& [_n, _d] : SUBMODULES(*this, type)) { \
		return &_d; \
	} \
	return nullptr; \
}\
const type* name() const { \
	for (auto&& [_n, _d] : SUBMODULES(*this, type)) { \
		return &_d; \
	} \
	return nullptr; \
}

#define G_QUICK_ACCESS_ITEM(type, name) \
type##Item* name(){ \
	const int child_count = this->childCount();\
	for (int i = 0; i < child_count; ++i) {\
		VariableItem* item = static_cast<VariableItem*>(this->child(i));\
		if (item->data_type_ == soap::VariableType::type) {\
			return static_cast<type##Item*>(item);\
		}\
	}\
	return nullptr;\
}

#define G_QUICK_ACCESS_TYPE(type, name, dtype) \
type* name(){ \
	for (auto&& [_n, _d] : SUBMODULES(*this, type)) { \
		if(_d.data_type_ == type::DataType::dtype){ \
			return &_d; \
		} \
	} \
	return nullptr; \
}\
const type* name() const { \
	for (auto&& [_n, _d] : SUBMODULES(*this, type)) { \
		if(_d.data_type_ == type::DataType::dtype){ \
			return &_d; \
		} \
	} \
	return nullptr; \
}

#define G_QUICK_ACCESS_ITEM_TYPE(type, name, dtype) \
type##Item* name(){ \
	const int child_count = this->childCount();\
	for (int i = 0; i < child_count; ++i) {\
		VariableItem* item = static_cast<VariableItem*>(this->child(i));\
		if (item->data_type_ == soap::VariableType::type) {\
			type* _d = static_cast<type*>(item->data_);\
			if (_d->data_type_ == type::DataType::dtype) {\
				return static_cast<type##Item*>(item);\
			}\
		}\
	}\
	return nullptr;\
}

#define CREATE_ROOT_MENU \
this->menus_["Root"] = new QMenu();\
using ItemType = std::remove_pointer_t<decltype(this)>;

#define ADD_MAIN_MENU(menu_id) this->menus_[menu_id] = this->menus_["Root"]->addMenu(menu_id); 

#define ADD_MENU(menu_id, menu_name, parent_menu_id) this->menus_[menu_id] = this->menus_[parent_menu_id]->addMenu(menu_name);

#define ADD_MAIN_ACTION(action_name, slot_name) this->menus_["Root"]->addAction(action_name, this, &ItemType::slot_name);

#define ADD_ACTION(action_name, parent_menu_id,  slot_name) this->menus_[parent_menu_id]->addAction(action_name, this, &ItemType::slot_name);

#define G_SET_IDENTIFIER(x) static inline QString g_identifier(){ return x;}

#define G_UNLOCK this->signal_emitter_->unlock(this->index_tree_); 

#define G_GETLOCK if(!this->__check_lock()) return; 


#define G_LOG(X) this->information_area_->log(X);
#define G_NOTICE(X) this->information_area_->notice(X);
#define G_WARN(X) this->information_area_->warn(X);

#define G_LINK_WORKER_THREAD(worker_type, worker_signal, receiver_type, receiver_slot) \
QThread* __t = new QThread(); \
worker->moveToThread(__t); \
connect(__t, &QThread::started, worker, &worker_type::run); \
connect(worker, &worker_type::worker_signal, this, &receiver_type::receiver_slot); \
connect(worker, &worker_type::x_message, this, &receiver_type::__s_receive_message); \
connect(worker, &worker_type::x_results_ready, this, &receiver_type::__s_work_finished); \
connect(worker, &worker_type::x_results_ready, __t, &QThread::quit); \
connect(worker, &worker_type::x_results_ready, worker, &QObject::deleteLater); \
connect(__t, &QThread::finished, __t, &QObject::deleteLater); \
__t->start(); 

#define G_LINK_WORKER_THREAD_NO_RESPONSE(worker_type, receiver_type) \
QThread* __t = new QThread(); \
worker->moveToThread(__t); \
connect(__t, &QThread::started, worker, &worker_type::run); \
connect(worker, &worker_type::x_message, this, &receiver_type::__s_receive_message); \
connect(worker, &worker_type::x_results_ready, this, &receiver_type::__s_work_finished); \
connect(worker, &worker_type::x_results_ready, __t, &QThread::quit); \
connect(worker, &worker_type::x_results_ready, worker, &QObject::deleteLater); \
connect(__t, &QThread::finished, __t, &QObject::deleteLater); \
__t->start(); 

#define G_TASK_LOG(message) \
emit x_message(message, 0);

#define G_TASK_WARN(message) \
emit x_message(message, 1);

#define G_TASK_NOTICE(message) \
emit x_message(message, 2);

#define G_TASK_END \
emit x_results_ready();\
return;

#define CHANGE_DATA(type, parent_type) \
if(this->attached_to(soap::VariableType::parent_type))\
{\
	auto* parent = this->trace_back<parent_type>(1);\
	parent->type##s_[new_name] = std::move(*this->data());\
	parent->type##s_.erase(this->title_);\
	this->data_ = &parent->type##s_[new_name];\
	return;\
}

#define CHANGE_DATA_IMPL0(type, ...)

#define CHANGE_DATA_IMPL1(type, _1) \
CHANGE_DATA(type, _1);

#define CHANGE_DATA_IMPL2(type, _1, _2) \
CHANGE_DATA(type, _1);\
CHANGE_DATA(type, _2);

#define CHANGE_DATA_IMPL3(type, _1, _2, _3) \
CHANGE_DATA(type, _1);\
CHANGE_DATA_IMPL2(type, _2, _3);

#define CHANGE_DATA_IMPL4(type, _1, _2, _3, _4) \
CHANGE_DATA(type, _1);\
CHANGE_DATA_IMPL3(type, _2, _3, _4);

#define CHANGE_DATA_IMPL5(type, _1, _2, _3, _4, _5) \
CHANGE_DATA(type, _1);\
CHANGE_DATA_IMPL4(type, _2, _3, _4, _5);

#define CHANGE_DATA_IMPL6(type, _1, _2, _3, _4, _5, _6) \
CHANGE_DATA(type, _1);\
CHANGE_DATA_IMPL5(type, _2, _3, _4, _5, _6);

#define CHANGE_DATA_IMPL7(type, _1, _2, _3, _4, _5, _6, _7) \
CHANGE_DATA(type, _1);\
CHANGE_DATA_IMPL6(type, _2, _3, _4, _5, _6, _7);

#define CHANGE_DATA_FUNCTION(type, ...)\
void __change_data_name(const QString& new_name) override {\
	CONCAT_INDIRECT(CHANGE_DATA_IMPL, NARGS(__VA_ARGS__))(type, ##__VA_ARGS__);\
}

#define REMOVE_DATA(type, parent_type) \
if (this->attached_to(soap::VariableType::parent_type)) {\
	auto* parent = this->trace_back<parent_type>(1);\
	parent->type##s_.erase(this->title_);\
	return;\
}

#define REMOVE_DATA_IMPL0(type, ...)

#define REMOVE_DATA_IMPL1(type, _1) \
REMOVE_DATA(type, _1);

#define REMOVE_DATA_IMPL2(type, _1, _2) \
REMOVE_DATA(type, _1);\
REMOVE_DATA(type, _2);

#define REMOVE_DATA_IMPL3(type, _1, _2, _3) \
REMOVE_DATA(type, _1);\
REMOVE_DATA_IMPL2(type, _2, _3);

#define REMOVE_DATA_IMPL4(type, _1, _2, _3, _4) \
REMOVE_DATA(type, _1);\
REMOVE_DATA_IMPL3(type, _2, _3, _4);

#define REMOVE_DATA_IMPL5(type, _1, _2, _3, _4, _5) \
REMOVE_DATA(type, _1);\
REMOVE_DATA_IMPL4(type, _2, _3, _4, _5);

#define REMOVE_DATA_IMPL6(type, _1, _2, _3, _4, _5, _6) \
REMOVE_DATA(type, _1);\
REMOVE_DATA_IMPL5(type, _2, _3, _4, _5, _6);

#define REMOVE_DATA_IMPL7(type, _1, _2, _3, _4, _5, _6, _7) \
REMOVE_DATA(type, _1);\
REMOVE_DATA_IMPL6(type, _2, _3, _4, _5, _6, _7);

#define REMOVE_DATA_FUNCTION(type, ...)\
void __remove_data() override {\
	CONCAT_INDIRECT(REMOVE_DATA_IMPL, NARGS(__VA_ARGS__))(type, ##__VA_ARGS__);\
	delete this->data();\
}

#define DATA_FUNCTION(type, ...)\
CHANGE_DATA_FUNCTION(type, ##__VA_ARGS__);\
REMOVE_DATA_FUNCTION(type, ##__VA_ARGS__);

#define G_ITEM_CONSTRUCTION(type, text1, ...) \
type##Item(\
const QString& title,\
IndexTree* index_tree,\
type* data,\
PlotsSuite* draw_suite,\
InformationTextBrowser* information_area,\
SignalEmitter* signal_emitter\
) :\
	VariableItem(\
		data,\
		soap::VariableType::type,\
		title,\
		index_tree,\
		QStringList() << title << text1,\
		draw_suite,\
		information_area,\
		signal_emitter\
	) {\
	this->index_tree_ = this->index_tree_->append(title, this->data_type_, this->data_, this); \
	this->__check_data(); \
	this->__set_menu(); \
	connect(this->signal_emitter_, &SignalEmitter::x_update_interface, this, &type##Item::__s_update_interface);\
	this->__s_update_interface();\
}\
type* data() { \
return static_cast<type*>(this->data_); \
} \
void __duplicate_this() override{ \
	type* ptr = new type(*this->data()); \
	this->signal_emitter_->x_data_create_soon(ptr, this->data_type_, "Duplicated"); \
};\
DATA_FUNCTION(type, ##__VA_ARGS__)