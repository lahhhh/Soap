#pragma once

#include "Identifier.h"

#include <QString>
#include <QSet>

#include "Custom.h"
#include "IndexTree.h"

class SignalEmitter : public QObject
{
	Q_OBJECT
public:
	explicit SignalEmitter(QObject* parent = nullptr);

	~SignalEmitter();

	QString get_unique_name(const QString& name);

	void remove_information(const QString& name);

	void update_information(const QString& name, soap::VariableType type, void* data, bool top_level = false);

	std::map< QString, QPair<soap::VariableType, void*> > variable_information_;

	std::map< QString, QPair<soap::VariableType, void*> > top_level_variables_;

	QSet<QString> occupied_names_;

	QMap<const IndexTree*, int> variable_lock_;

	QList<IndexTree*> index_;

	int suffix_ = 0;

	QWidget* widget_;

	void* get_variable(const QString& name);

	void remove_index(IndexTree* index_tree);

	void lock(IndexTree* index_tree);

	bool is_lock_free(IndexTree* index_tree) const;

	bool try_lock(IndexTree* index_tree);

	bool try_lock(const QList<IndexTree*>& data);

	void unlock(IndexTree* index_tree);

	void unlock(const QList<IndexTree*>& data);

	void* get_item(void* data);

	void* get_item(const QString& name);

	IndexTree* new_index_tree();

	IndexTree* search(void* data);

	QList<IndexTree*> search(const QList<void*>& data);

	IndexTree* search(const QString& data);

	QList<IndexTree*> search(const QStringList& data);

	QMap<QString, void*> get_type_variable(soap::VariableType type) const;
	std::map<QString, void*> get_type_variable_std(soap::VariableType type) const;

signals:

	void x_log(QString);

	void x_notice(QString);

	void x_warn(QString);

	void x_update_interface();

	void x_data_create_soon(void* data, soap::VariableType, QString name);

	void x_data_delete_soon(void* data, soap::VariableType = soap::VariableType::AnyVariable, void* item = nullptr);

	void x_data_edit_soon(void* data, soap::VariableType = soap::VariableType::AnyVariable, void* item = nullptr);
};

