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

QString IndexTree::search_name(void* data) {
	IndexTree* it = this->search(data);

	if (it != nullptr) {
		return it->name_;
	}
	else {
		return QString();
	}
};

IndexTree* IndexTree::sub_search(void* data) {
	if (this->data_ == data) {
		return this;
	}

	for (auto& child : this->children_) {
		auto res = child.second.sub_search(data);
		if (res != nullptr) {
			return res;
		}
	}

	return nullptr;
};

IndexTree* IndexTree::search(void* data) {
	IndexTree* root = this->root();

	return root->sub_search(data);
};

QStringList IndexTree::sub_search_type(soap::VariableType type) {
	QStringList res;

	if (this->type_ == type) {
		res << this->name_;
	}

	for (auto& child : this->children_) {
		auto sub_res = child.second.sub_search_type(type);
		if (!sub_res.isEmpty()) {
			res << sub_res;
		}
	}

	return res;
};

QStringList IndexTree::search_type(soap::VariableType type) {
	IndexTree* root = this->root();

	return root->sub_search_type(type);
};

IndexTree* IndexTree::sub_search(const QString& name) {

	if (this->name_ == name) {
		return this;
	}

	for (auto&& child : this->children_) {

		auto res = child.second.sub_search(name);
		if (res != nullptr) {
			return res;
		}
	}

	return nullptr;
};

IndexTree* IndexTree::search(const QString& name) {

	IndexTree* root = this->root();

	return root->sub_search(name);
};