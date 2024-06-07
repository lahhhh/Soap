#include "SignalEmitter.h"

#include <QWidget>

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"
#include "DataFrame.h"
#include "StringVector.h"

SignalEmitter::SignalEmitter(QObject* parent) : QObject(parent)
{
	this->widget_ = new QWidget();
}

IndexTree* SignalEmitter::new_index_tree() {

	IndexTree* it = new IndexTree();
	
	this->index_ << it;
	
	return it;
};

void SignalEmitter::remove_index(IndexTree* index_tree) {

	QString name = index_tree->name_;
	
	auto it = this->search(name);

	if (it == nullptr) {
		return;
	}

	if (it->is_root()) {
		delete it;
		this->index_.removeOne(it);
	}
	else {
		it->parent_->children_.erase(it->name_);
	}
};

void SignalEmitter::lock(IndexTree* index_tree) {

	while (index_tree != nullptr) {
	
		this->variable_lock_[index_tree] = 1;
		
		index_tree = index_tree->parent_;
	}
};

bool SignalEmitter::is_lock_free(IndexTree* index_tree) const{

	while (index_tree != nullptr) {
	
		if (this->variable_lock_[index_tree] > 0) {
			return false;
		}
		
		index_tree = index_tree->parent_;
	}

	return true;
};

QList<IndexTree*> SignalEmitter::search(const QStringList& data) {

	return _Cs sapply(data, [this](const QString& data) {return this->search(data); });
};

IndexTree* SignalEmitter::search(const QString& data) {

	for (auto it : this->index_) {
	
		auto res = it->search(data);
		
		if (res != nullptr) {
			return res;
		}
	}
	return nullptr;
};

QList<IndexTree*> SignalEmitter::search(const QList<void*>& data) {

	return _Cs sapply(data, [this](void* data) {return this->search(data); });
};

void* SignalEmitter::get_item(const QString& name) {

	for (auto it : this->index_) {

		auto res = it->search(name);

		if (res != nullptr) {
			return res->item_;
		}
	}

	return nullptr;
};

void* SignalEmitter::get_item(void* data) {

	for (auto it : this->index_) {

		auto res = it->search(data);

		if (res != nullptr) {
			return res->item_;
		}
	}

	return nullptr;
};

IndexTree* SignalEmitter::search(void* data) {

	for (auto it : this->index_) {
	
		auto res = it->search(data);
		
		if (res != nullptr) {
			return res;
		}
	}

	return nullptr;
};

bool SignalEmitter::try_lock(IndexTree* index_tree) {

	if (!this->is_lock_free(index_tree)) {
		return false;
	}

	this->lock(index_tree);

	return true;
};

bool SignalEmitter::try_lock(const QList<IndexTree*>& data) {

	for (auto ptr : data) {
		if (!this->is_lock_free(ptr)) {
			return false;
		}
	}

	for (auto ptr : data) {
		this->lock(ptr);
	}

	return true;
};

void SignalEmitter::unlock(IndexTree* index_tree) {

	while (index_tree != nullptr) {
	
		this->variable_lock_[index_tree] = 0;
		
		index_tree = index_tree->parent_;
	}
};

void SignalEmitter::unlock(const QList<IndexTree*>& data) {

	for (const auto ptr : data) {
	
		this->unlock(ptr);
	}
};

std::map<QString, void*> SignalEmitter::get_type_variable_std(soap::VariableType type) const {

	std::map<QString, void*> ret;

	for (auto&& [name, info] : this->variable_information_) {
	
		if (info.first == type) {
		
			ret[name] = info.second;
		}
	}

	return ret;
};

QMap<QString, void*> SignalEmitter::get_type_variable(soap::VariableType type) const{

	QMap<QString, void*> ret;

	for (auto&& [name, info] : this->variable_information_) {

		if (info.first == type) {
		
			ret[name] = info.second;
		}
	}

	return ret;
};

SignalEmitter::~SignalEmitter() {

	delete this->widget_;

	for (auto&& [name, info] : this->top_level_variables_) {

		switch (info.first)
		{
		case soap::VariableType::SingleCellRna:
			delete static_cast<SingleCellRna*>(info.second);
			break;

		case soap::VariableType::SingleCellMultiome:
			delete static_cast<SingleCellMultiome*>(info.second);
			break;

		case soap::VariableType::DataField:
			delete static_cast<DataField*>(info.second);
			break;

		case soap::VariableType::DifferentialAnalysis:
			delete static_cast<DifferentialAnalysis*>(info.second);
			break;

		case soap::VariableType::DenseDouble:
			delete static_cast<DenseDouble*>(info.second);
			break;

		case soap::VariableType::DenseInt:
			delete static_cast<DenseInt*>(info.second);
			break;

		case soap::VariableType::SparseDouble:
			delete static_cast<SparseDouble*>(info.second);
			break;

		case soap::VariableType::SparseInt:
			delete static_cast<SparseInt*>(info.second);
			break;

		case soap::VariableType::Metadata:
			delete static_cast<Metadata*>(info.second);
			break;

		case soap::VariableType::Embedding:
			delete static_cast<Embedding*>(info.second);
			break;

		case soap::VariableType::Enrichment:
			delete static_cast<Enrichment*>(info.second);
			break;

		case soap::VariableType::DataFrame:
			delete static_cast<DataFrame*>(info.second);
			break;

		case soap::VariableType::GSEA:
			delete static_cast<GSEA*>(info.second);
			break;

		case soap::VariableType::CellChat:
			delete static_cast<CellChat*>(info.second);
			break;

		case soap::VariableType::CNV:
			delete static_cast<CNV*>(info.second);
			break;

		case soap::VariableType::Pando:
			delete static_cast<Pando*>(info.second);
			break;

		case soap::VariableType::GenomicRange:
			delete static_cast<GenomicRange*>(info.second);
			break;

		case soap::VariableType::MotifPosition:
			delete static_cast<MotifPosition*>(info.second);
			break;

		case soap::VariableType::Footprint:
			delete static_cast<Footprint*>(info.second);
			break;

		case soap::VariableType::CoverageTrack:
			delete static_cast<CoverageTrack*>(info.second);
			break;

		case soap::VariableType::Fragments:
			delete static_cast<Fragments*>(info.second);
			break;

		case soap::VariableType::VelocytoBase:
			delete static_cast<VelocytoBase*>(info.second);
			break;

		case soap::VariableType::VelocityEstimate:
			delete static_cast<VelocityEstimate*>(info.second);
			break;

		case soap::VariableType::ScveloEstimate:
			delete static_cast<ScveloEstimate*>(info.second);
			break;

		case soap::VariableType::StringVector:
			delete static_cast<StringVector*>(info.second);
			break;

		case soap::VariableType::IntegerVector:
			//delete static_cast<SingleCellRna*>(info.second);
			break;

		case soap::VariableType::NumericVector:
			//delete static_cast<NumericVector*>(info.second);
			break;

		case soap::VariableType::GeneName:
			delete static_cast<GeneName*>(info.second);
			break;

		default:
			break;
		}

	}
};

QString SignalEmitter::get_unique_name(const QString& name) {

	if (!this->variable_information_.contains(name) && 
		!this->occupied_names_.contains(name)
		){

		this->occupied_names_.insert(name);

		return name;
	}

	QString new_name = name.split(" | ")[0], ret = new_name;

	do {
		ret = new_name + " | " + QString::number(this->suffix_++);
	} 
	while (this->variable_information_.contains(ret) || this->occupied_names_.contains(ret));

	this->occupied_names_.insert(ret);

	return ret;
};

void SignalEmitter::update_information(const QString& name, soap::VariableType type, void* data, bool top_level) {

	this->variable_information_[name] = qMakePair(type, data);
	
	if (top_level) {
		this->top_level_variables_[name] = qMakePair(type, data);
	}
	
	this->occupied_names_.remove(name);
};

void SignalEmitter::remove_information(const QString& name) {

	this->variable_information_.erase(name);
	
	this->top_level_variables_.erase(name);
};

void* SignalEmitter::get_variable(const QString& name) {

	auto iter = this->variable_information_.find(name);

	if (iter != this->variable_information_.end()) {
		return iter->second.second;
	}
	else {
		return nullptr;
	}
};
