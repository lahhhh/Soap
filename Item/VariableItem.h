#pragma once

#include <QTreeWidgetItem>

#include "Identifier.h"

#include "Custom.h"

#include "IndexTree.h"
#include "SignalEmitter.h"
#include "PlotsSuite.h"
#include "InformationTextBrowser.h"


class VariableItem :
    public QObject, public QTreeWidgetItem
{

public:
    VariableItem() = default;
    VariableItem(const VariableItem&) = delete;

    VariableItem(
        void* data,
        soap::VariableType data_type,
        const QString& title, 
        IndexTree* index_tree,
        const QStringList& item,         
        PlotsSuite* draw_suite, 
        InformationTextBrowser* information_area, 
        SignalEmitter* signal_emitter
    );

    virtual void __remove_information();

    virtual void __set_menu();

    virtual void __show_this();

    virtual void __data_delete_soon();

    virtual void __remove_data();

    virtual void __remove_this();

    virtual void __check_data() ;

    virtual void __clear_reserve_data();

    virtual void __change_data_name(const QString& new_name);

    virtual void __duplicate_this();

    bool __check_lock();

    virtual ~VariableItem() = default;

    template<typename TraceType>
    TraceType* trace_back(int trace_layer) const{

        return this->index_tree_->trace_back<TraceType>(trace_layer);
    }

    template<typename RootItemType>
    RootItemType* get_item(void* data) const {

        return static_cast<RootItemType*>(this->signal_emitter_->get_item(data));
    }

    template<typename RootType>
    RootType* get_root() const {

        return this->index_tree_->get_root<RootType>();
    }

    template<typename RootItemType>
    RootItemType* get_root_item() const {

        return static_cast<RootItemType*>(this->signal_emitter_->get_item(this->get_root<int>()));
    }

    bool attached_to(soap::VariableType type) const;

    bool stem_from(soap::VariableType type) const;

    bool is_atomic() const;

    void set_item(VariableItem* item);

	template <typename VariableType>
	void check_variable(std::map<QString, VariableType>& variable_map) {

        auto constexpr data_type = soap::type<VariableType>();

		QStringList original_titles = _Cs keys(variable_map);

		for (const auto& title : original_titles) {

			QString new_title = this->signal_emitter_->get_unique_name(title);

			if (new_title != title) {

				variable_map[new_title] = std::move(variable_map[title]);
				variable_map.erase(title);
			}

			auto data = &(variable_map[new_title]);

			auto item = new soap::item_type<VariableType>(
                new_title, 
                this->index_tree_, 
                data, this->draw_suite_, 
                this->information_area_, 
                this->signal_emitter_
            );

			this->set_item(item);
		}
	}

public slots:

    virtual void __s_update_interface();

    virtual void __s_rename();

    virtual void __s_duplicate();

    virtual void __s_delete_this();

    virtual void __s_export_as_item();

    virtual void __s_export_as_csv();

    virtual void __s_export_as_tsv();

    virtual void __s_receive_message(QString message, int mode);

    virtual void __s_work_finished();


public:  
    
    void* data_{ nullptr };
    soap::VariableType data_type_{ soap::VariableType::AnyVariable };

    QString title_;

    IndexTree* index_tree_{ nullptr };
    PlotsSuite* draw_suite_{ nullptr };
    InformationTextBrowser* information_area_{ nullptr };
    SignalEmitter* signal_emitter_{ nullptr };

    QMap<QString, QMenu*> menus_;
};

