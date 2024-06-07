#include "MotifLocateWorker.h"

#include "Custom.h"
#include "FileIO.h"
#include "GenomeUtility.h"

#include "MOODS/moods_interface.h"

void MotifLocateWorker::run() {

	this->get_peak();

	if (!this->get_peak_sequence()) {
		G_TASK_WARN("Motif location Failed.");
		G_TASK_END;
	}

	this->get_base_background();

	this->get_result();

	emit x_motif_location_ready(this->mp_);

	G_TASK_END;
}

void MotifLocateWorker::get_result() {

	int n_peak = this->peak_names_.size();

	this->pattern_database_ = read_motif_database(this->database_name_);

	int n_motif = this->pattern_database_.size();

	MotifPosition mp(n_peak, n_motif);

	mp.peak_locations_ = this->peak_;
	mp.motifs_ = read_motif_database(this->database_name_);
	mp.peak_names_ = this->peak_names_;
	mp.motif_names_ = _Cs keys(mp.motifs_);

	match_motif(mp, this->peak_sequences_, this->sequence_background_, 5e-5, 7);

	this->mp_ = mp;
};

void MotifLocateWorker::get_base_background() {

	this->sequence_background_ = _Cs get_nucleic_acid_frequency(this->peak_sequences_);
}

bool MotifLocateWorker::get_peak_sequence() {
	
	if (this->species_ == soap::Species::Human) {
		this->genome_file_.set_sequence_file(FILE_HUMAN_GRCH38_2BIT);
	}
	else {
		G_TASK_NOTICE("Now only human genome is supported.");
		return false;
	}
	const qsizetype size = this->peak_.size();
	qsizetype miss_count = 0;

	for (qsizetype i = 0; i < size; ++i) {
		
		auto [sequence_name, start, end, strand] = this->peak_.at(i);
		QString sequence = this->genome_file_.get_sequence(sequence_name, start, end);
		if (sequence.isEmpty()) {
			++miss_count;
		}
		this->peak_sequences_ << sequence;
	}

	if (miss_count >= size / 2) {
		return false;
	}
	else {
		return true;
	}
}

void MotifLocateWorker::get_peak() {

	this->peak_ = _Cs stringlist_to_genomic_range(this->peak_names_);
};