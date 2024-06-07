#include "Fragments.h"

#include "Custom.h"

const std::pair< std::vector<int>, std::vector<int> >& 
Fragments::get_fragments(
    const QString& sequence_name,
    const QString& cell_name
) const {
    return this->data_.at(sequence_name)[this->cell_names_.indexOf(cell_name)];
};

void Fragments::finalize() {

	for (auto& [chr, data] : this->data_) {
		for (auto& [start, end] : data) {
			if (!start.empty()) {
				_Cs sort_by_first(start, end);
			}
		}
	}
};

void Fragments::adjust_length_by_cell_name_length() {
	const qsizetype length = this->cell_names_.size();

	for (auto& [chr, sequence] : this->data_) {
		sequence.resize(length);
	}
};