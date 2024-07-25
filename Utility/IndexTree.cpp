#include "IndexTree.h"

IndexTree::IndexTree(
	const QString& name, 
	soap::VariableType type,
	void* data,
	void* item,
	IndexTree* parent
)
	:
	name_(name),
	type_(type),
	data_(data),
	item_(item),
	parent_(parent)
{};

IndexTree* IndexTree::append(
	const QString& name, 
	soap::VariableType type, 
	void* data,
	void* item
) {
	if (this->data_ == nullptr) {
		this->name_ = name;
		this->type_ = type;
		this->data_ = data;
		this->item_ = item;
		this->parent_ = nullptr;
		return this;
	}

	this->children_[name] = IndexTree(name, type, data, item, this);
	return &this->children_[name];
};

bool IndexTree::is_root() const {
	return this->parent_ == nullptr;
};

QString IndexTree::search_name(void* data, bool from_root) {

	IndexTree* it = this->search(data, from_root);

	if (it != nullptr) {
		return it->name_;
	}
	else {
		return {};
	}
};

IndexTree* IndexTree::search(void* data, bool from_root) {

	IndexTree* start = this;

	if (from_root) {
		start = this->root();
	}

	if (start->data_ == data) {
		return start;
	}

	for (auto& child : start->children_) {
		auto res = child.second.search(data, false);
		if (res != nullptr) {
			return res;
		}
	}

	return nullptr;
};

void* IndexTree::search_one(soap::VariableType type, bool from_root) {

	IndexTree* start = this;

	if (from_root) {
		start = this->root();
	}

	if (start->type_ == type) {
		return start->data_;
	}

	for (auto& child : start->children_) {
		auto sub_res = child.second.search_one(type, false);
		if (sub_res != nullptr) {
			return sub_res;
		}
	}

	return nullptr;
};

QStringList IndexTree::search(soap::VariableType type, bool from_root) {

	IndexTree* start = this;

	if (from_root) {
		start = this->root();
	}

	QStringList res;

	if (start->type_ == type) {
		res << start->name_;
	}

	for (auto& child : start->children_) {
		auto sub_res = child.second.search(type, false);
		if (!sub_res.isEmpty()) {
			res << sub_res;
		}
	}

	return res;
};

IndexTree* IndexTree::search(const QString& name, bool from_root) {

	IndexTree* start = this;

	if (from_root) {
		start = this->root();
	}

	if (start->name_ == name) {
		return start;
	}

	for (auto&& child : start->children_) {

		auto res = child.second.search(name, false);
		if (res != nullptr) {
			return res;
		}
	}

	return nullptr;
};