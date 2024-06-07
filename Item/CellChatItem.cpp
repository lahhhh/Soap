#include "CellChatItem.h"

#include "MatrixWindow.h"
#include "CommonDialog.h"
#include "ItemIOWorker.h"

void CellChatItem::__set_menu() {

	CREATE_ROOT_MENU;

	ADD_MAIN_MENU("Statistics");

	ADD_ACTION("Number of interactions", "Statistics", s_show_interaction_numbers);
	ADD_ACTION("Weights of interactions", "Statistics", s_show_interaction_weights);
	ADD_ACTION("Interaction probability", "Statistics", s_show_interaction_probability);
	ADD_ACTION("Interaction P value", "Statistics", s_show_interaction_p_value);
	ADD_ACTION("Pathway probability", "Statistics", s_show_pathway_probability);

	ADD_MAIN_MENU("Export");

	ADD_ACTION("as Item", "Export", __s_export_as_item);

	ADD_MAIN_ACTION("Delete", __s_delete_this);
}

void CellChatItem::s_show_interaction_numbers() {
	MatrixWindow::show_matrix(
		&this->data()->counts_,
		this->data()->levels_,
		this->data()->levels_,
		"Number of interactions",
		this->signal_emitter_);
};

void CellChatItem::s_show_interaction_weights() {
	MatrixWindow::show_matrix(
		&this->data()->weights_,
		this->data()->levels_,
		this->data()->levels_,
		"Weights of interactions",
		this->signal_emitter_);
};

void CellChatItem::s_show_interaction_probability() {
	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::ComboBox},
		{ this->data()->ligand_receptor_index_}
	);

	if (res.isEmpty())return;

	int index = this->data()->ligand_receptor_index_.indexOf(res[0]);
	MatrixWindow::show_matrix(
		&this->data()->ligand_receptor_probability_[index],
		this->data()->levels_,
		this->data()->levels_,
		"Probability of " + res[0],
		this->signal_emitter_);
};

void CellChatItem::s_show_interaction_p_value() {
	QStringList res = CommonDialog::get_response(
		this->signal_emitter_, 
		"Select Interaction",
		{ "Interaction" },
		{ soap::InputStyle::ComboBox},
		{ this->data()->ligand_receptor_index_}
	);

	if (res.isEmpty())return;

	int index = this->data()->ligand_receptor_index_.indexOf(res[0]);

	MatrixWindow::show_matrix(
		&this->data()->ligand_receptor_p_value_[index],
		this->data()->levels_,
		this->data()->levels_,
		"P value of " + res[0],
		this->signal_emitter_);
};

void CellChatItem::s_show_pathway_probability() {
	QStringList res = CommonDialog::get_response(
		this->signal_emitter_,
		"Select Pathway",
		{ "Pathway" },
		{ soap::InputStyle::ComboBox},
		{ this->data()->pathway_index_}
	);

	if (res.isEmpty())return;

	MatrixWindow::show_matrix(
		&this->data()->pathway_probability_[res[0]],
		this->data()->levels_,
		this->data()->levels_,
		"Probability of " + res[0],
		this->signal_emitter_);
};
