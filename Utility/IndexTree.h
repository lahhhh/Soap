#pragma once

#include "Identifier.h"

#include <QVector>
#include <map>


class IndexTree
{
public:
	IndexTree() = default;
	IndexTree(const IndexTree&) = default;
	IndexTree(IndexTree&&) = default;
	IndexTree& operator=(const IndexTree&) = default;
	IndexTree& operator=(IndexTree&&) = default;
	~IndexTree() = default;

	IndexTree(
		const QString& name,
		soap::VariableType type,
		void* data,
		void* item,
		IndexTree* parent
	);

	inline bool operator==(const IndexTree& rhs) const{
		return this->data_ == rhs.data_;
	}

	QString name_;

	soap::VariableType type_;

	std::map<QString, IndexTree> children_;

	void* data_ = nullptr;
	void* item_ = nullptr;

	IndexTree* parent_ = nullptr;

	IndexTree* append(const QString& name,
		soap::VariableType type,
		void* data,
		void* item
	);

	template<typename TraceType>
	TraceType* trace_back(int trace_layer) const {

		const IndexTree* data = this;

		for (int i = 0; i < trace_layer; ++i) {
			data = data->parent_;
		}

		return static_cast<TraceType*>(data->data_);
	}

	template<typename RootType>
	RootType* get_root() const {

		const IndexTree* root = this;

		while (root->parent_ != nullptr) {
			root = root->parent_;
		}

		return static_cast<RootType*>(root->data_);
	}

	inline soap::VariableType trace_back_type(int trace_layer) const {
		const IndexTree* data = this;

		for (int i = 0; i < trace_layer; ++i) {
			data = data->parent_;
		}

		return data->type_;
	}

	inline IndexTree* root() {
		IndexTree* root = this;

		while (root->parent_ != nullptr) {
			root = root->parent_;
		}

		return root;
	}

	IndexTree* sub_search(void* data);

	IndexTree* search(void* data);

	IndexTree* sub_search(const QString& name);

	IndexTree* search(const QString& name);

	QStringList sub_search_type(soap::VariableType type);

	QStringList search_type(soap::VariableType type);

	QString search_name(void* data);

	bool is_root() const;
};

