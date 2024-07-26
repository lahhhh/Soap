#include "MacsCallPeakWorker.h"

#include <zlib.h>

#include "custom.h"
#include "Poisson.h"

void fix_coordinates(QVector<int>& vec) {
	int size = vec.size();
	for (int i = 0; i < size; ++i) {
		if (vec[i] < 0) {
			vec[i] = 0;
		}
		else {
			break;
		}
	}
}

struct PeakContent {
	QVector<int> above_cutoff_start_, above_cutoff_end_, above_cutoff_index_;
	QVector<double> treat_, control_;

	PeakContent() = default;

	PeakContent(int ts, int te, double tp, double cp, int ti) {
		this->above_cutoff_start_ << ts;
		this->above_cutoff_end_ << te;
		this->treat_ << tp;
		this->control_ << cp;
		this->above_cutoff_index_ << ti;
	}

	void append(int ts, int te, double tp, double cp, int ti) {
		this->above_cutoff_start_ << ts;
		this->above_cutoff_end_ << te;
		this->treat_ << tp;
		this->control_ << cp;
		this->above_cutoff_index_ << ti;
	}

	void clear() {
		this->above_cutoff_start_.clear();
		this->above_cutoff_end_.clear();
		this->treat_.clear();
		this->control_.clear();
		this->above_cutoff_index_.clear();
	}

	void refresh(int ts, int te, double tp, double cp, int ti) {
		this->clear();

		this->above_cutoff_start_ << ts;
		this->above_cutoff_end_ << te;
		this->treat_ << tp;
		this->control_ << cp;
		this->above_cutoff_index_ << ti;
	}

	qsizetype size() const {
		return this->above_cutoff_start_.size();
	}

	int peak_length() const {
		qsizetype size = this->size();
		if (size == 0) {
			return 0;
		}
		return this->above_cutoff_end_[size - 1] - this->above_cutoff_start_[0];
	}

	int start() const {
		return this->above_cutoff_start_[0];
	}

	int end() const {
		return this->above_cutoff_end_.last();
	}
};
static PeakContent peak_content;

GenomicRange MacsCallPeakWorker::call_peak(const QList<const Fragments*>& fragments_objects) {
	MacsCallPeakWorker worker(fragments_objects);

	if (!worker.detect_tag_size()) {
		return GenomicRange();
	}
	if (!worker.bed_file_names_.isEmpty()) {
		if (!worker.parse_bed_file()) {
			return GenomicRange();
		}
	}

	worker.parse_fragments_object();
	worker.filter_duplicates();

	try {
		worker.detect_peaks();
	}
	catch (...) {
		return GenomicRange();
	}

	GenomicRange genomic_range(worker.seqnames_, worker.ranges_, worker.strand_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_RANGE_NAME, worker.name_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_SCORE, worker.score_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_FOLD_CHANGE, worker.fold_change_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_P_VALUE, worker.p_value_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_Q_VALUE, worker.q_value_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_SUMMIT_POSITION, worker.summit_position_);

	genomic_range.finalize();

	return genomic_range;
};

void MacsCallPeakWorker::run() {
	if (!detect_tag_size()) {
		G_TASK_WARN("Fragments File is broken.");
		G_TASK_END;
	}

	if (!this->bed_file_names_.isEmpty()) {
		if (!parse_bed_file()) {
			G_TASK_WARN("Fragments File is broken.");
			G_TASK_END;
		}
	}

	this->parse_fragments_object();

	this->filter_duplicates();

	try {
		this->detect_peaks();
	}
	catch (std::exception& e) {
		G_TASK_WARN(QString::fromUtf8(e.what()));
		G_TASK_WARN("Problem Meeted. Please Save Your Data.");
		G_TASK_END;
	}
	catch (...) {
		G_TASK_WARN("Problem Meeted. Please Save Your Data.");
		G_TASK_END;
	}

	GenomicRange genomic_range(this->seqnames_, this->ranges_, this->strand_);

	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_RANGE_NAME, this->name_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_SCORE, this->score_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_FOLD_CHANGE, this->fold_change_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_P_VALUE, this->p_value_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_Q_VALUE, this->q_value_);
	genomic_range.metadata_.update(METADATA_GENOMIC_RANGE_MACS_SUMMIT_POSITION, this->summit_position_);

	genomic_range.finalize();

	emit x_genomic_range_ready(genomic_range);
	G_TASK_END;

}

void MacsCallPeakWorker::detect_peaks() {
	this->maximum_gap_ = this->tag_size_;
	this->call_peak_without_control();
};

void MacsCallPeakWorker::chromosome_call_peak(const QString& chromosome) {

	G_TASK_LOG("[MACS] Calling peaks in " + chromosome);

	auto& [p, treat, control] = chromosome_position_treat_control_[chromosome];
	auto score = calculate_q_score(treat, control);

	auto above_cutoff = _Cs greater_than(score, this->log_q_value_);

	if (above_cutoff.count() == 0)return;

	QVector<int> above_cutoff_index = _Cs which(above_cutoff);
	QVector<int> above_cutoff_end_position = _Cs sliced(p, above_cutoff);

	if (above_cutoff_index[0] == 0) {
		above_cutoff_index[0] = 1;
	}

	QVector<int> above_cutoff_start_position = _Cs reordered(p, _Cs minus(above_cutoff_index, 1));

	const int
		* pointer_to_above_cutoff_start = above_cutoff_start_position.data(),
		* pointer_to_above_cutoff_end = above_cutoff_end_position.data(),
		* pointer_to_above_cutoff_index = above_cutoff_index.data();

	const double
		* ptr_treat = treat.data(),
		* ptr_control = control.data();

	peak_content.append(
		*pointer_to_above_cutoff_start,
		*pointer_to_above_cutoff_end,
		ptr_treat[*pointer_to_above_cutoff_index],
		ptr_control[*pointer_to_above_cutoff_index],
		*pointer_to_above_cutoff_index
	);

	int last_p = *pointer_to_above_cutoff_end;

	++pointer_to_above_cutoff_start;
	++pointer_to_above_cutoff_end;
	++pointer_to_above_cutoff_index;

	int length = above_cutoff_start_position.size();


	for (int i = 1; i < length; ++i) {

		int tl = *pointer_to_above_cutoff_start - last_p;

		if (tl <= this->maximum_gap_) {

			peak_content.append(
				*pointer_to_above_cutoff_start,
				*pointer_to_above_cutoff_end,
				ptr_treat[*pointer_to_above_cutoff_index],
				ptr_control[*pointer_to_above_cutoff_index],
				*pointer_to_above_cutoff_index
			);
		}
		else {
			close_peak_without_subpeaks(chromosome, score);

			peak_content.refresh(
				*pointer_to_above_cutoff_start,
				*pointer_to_above_cutoff_end,
				ptr_treat[*pointer_to_above_cutoff_index],
				ptr_control[*pointer_to_above_cutoff_index],
				*pointer_to_above_cutoff_index
			);
		}

		last_p = *pointer_to_above_cutoff_end;

		++pointer_to_above_cutoff_start;
		++pointer_to_above_cutoff_end;
		++pointer_to_above_cutoff_index;
	}

	if (peak_content.size() > 0) {
		close_peak_without_subpeaks(chromosome, score);
	}
};

void MacsCallPeakWorker::close_peak_without_subpeaks(const QString& chromosome, const QVector<double>& score) {

	static int peak_count = 1;

	int peak_length = peak_content.peak_length();

	if (peak_length >= this->minimum_length_) {

		auto& [above_cutoff_start, above_cutoff_end, above_cutoff_index, treat, control] = peak_content;

		const qsizetype size = peak_content.size();

		double summit_value = 0;

		QVector<int> tsummit, tsummit_index;

		for (qsizetype i = 0; i < size; ++i) {

			int tstart = above_cutoff_start[i], tend = above_cutoff_end[i];
			double treat_score = treat[i];

			if (!summit_value || summit_value < treat_score) {
				tsummit << (tstart + tend) / 2;
				tsummit_index << i;
				summit_value = treat_score;
			}
			else if (summit_value == treat_score) {
				tsummit << (tstart + tend) / 2;
				tsummit_index << i;
			}
		}

		int middle_index = (tsummit.size() + 1) / 2 - 1;
		int summit_position = tsummit[middle_index];
		int summit_index = tsummit_index[middle_index];
		double summit_treat = treat[summit_index];
		double summit_control = control[summit_index];

		if ((summit_treat + 1) / (summit_control + 1) < 1 || this->log_q_value_ > score[peak_content.above_cutoff_index_[summit_index]]) {
			return;
		}

		double summit_p_score = this->p_score_map_[qMakePair((int)summit_treat, summit_control)];
		double summit_q_score = this->p_q_map_[summit_p_score];

		this->seqnames_ << chromosome;
		this->strand_ << '*';
		this->ranges_.append(peak_content.start(), peak_content.peak_length());
		this->name_ << "soap" + QString::number(peak_count);
		this->fold_change_ << (summit_treat + 1) / (summit_control + 1);
		this->p_value_ << summit_p_score;
		this->q_value_ << summit_q_score;
		this->score_ << (int)(summit_q_score * 10);
		this->summit_position_ << summit_position;

		++peak_count;
	}
};

QVector<double> MacsCallPeakWorker::calculate_q_score(const QVector<double>& treat, const QVector<double>& control) {
	const qsizetype size = treat.size();

	QVector<double> ret(size);

	for (int i = 0; i < size; ++i) {
		ret[i] = this->p_q_map_[get_p_score({ int(treat[i]), control[i] })];
	}

	return ret;
};

void MacsCallPeakWorker::call_peak_without_control() {
	this->minimum_length_ = this->d_ = this->ext_size_;
	this->lambda_background_ = (double)this->d_ * this->total_ / this->genome_size_;

	this->calculate_p_q_map();

	G_TASK_LOG("[MACS] Calling peaks...");

	for (auto iter = this->locations_.cbegin(); iter != this->locations_.cend(); ++iter) {
		
		QString chromosome = iter->first;

		this->chromosome_call_peak(chromosome);
	}
};

std::pair<QVector<int>, QVector<double>>
MacsCallPeakWorker::calculate_distribution(
	const QString& chromosome,
	int d,
	double scaling_factor,
	double baseline_value,
	bool directional,
	int end_shift)
{
	int five_shift, three_shift;
	if (directional) {
		five_shift = -end_shift;
		three_shift = d + end_shift;
	}
	else {
		five_shift = d / 2 - end_shift;
		three_shift = end_shift + d - d / 2;
	}

	QVector<int> start_position, end_position;
	auto& chr = this->locations_[chromosome];
	start_position << _Cs minus(chr[0], five_shift) << _Cs minus(chr[1], three_shift);
	end_position << _Cs add(chr[0], three_shift) << _Cs add(chr[1], five_shift);
	std::ranges::sort(start_position);
	std::ranges::sort(end_position);

	fix_coordinates(start_position);
	fix_coordinates(end_position);

	qsizetype length = start_position.size();
	QVector<int> return_position(2 * length + 1, 0); // index + 1
	QVector<double> return_value(2 * length + 1, 0);
	qsizetype index_start = 0, index_end = 0;

	int* pointer_to_return_position = return_position.data(), * pointer_to_start_position = start_position.data(), * pointer_to_end_position = end_position.data();
	double* pointer_to_return_value = return_value.data(), * val_start = return_value.data();
	qsizetype index = 0;
	int pile_up = 0;
	int previous_position = std::min(start_position[0], end_position[0]);
	if (previous_position != 0) {
		return_position[0] = previous_position;
		return_value[0] = std::max(0., baseline_value);
		++pointer_to_return_position;
		++pointer_to_return_value;
		++index;
	}

	previous_position = 0;
	while (index_start < length && index_end < length) {
		if (*pointer_to_start_position < *pointer_to_end_position) {
			int p = *pointer_to_start_position;
			if (p != previous_position) {
				++index;
				*pointer_to_return_position = p;
				*pointer_to_return_value = std::max(pile_up * scaling_factor, baseline_value);

				++pointer_to_return_position;
				++pointer_to_return_value;

				previous_position = p;
			}
			++pile_up;
			++index_start;
			++pointer_to_start_position;
		}
		else if (*pointer_to_start_position > *pointer_to_end_position) {
			int p = *pointer_to_end_position;
			if (p != previous_position) {
				++index;
				*pointer_to_return_position = p;
				*pointer_to_return_value = std::max(pile_up * scaling_factor, baseline_value);

				++pointer_to_return_position;
				++pointer_to_return_value;
				previous_position = p;
			}
			--pile_up;
			++index_end;
			++pointer_to_end_position;
		}
		else {
			++index_start;
			++index_end;
			++pointer_to_start_position;
			++pointer_to_end_position;
		}
	}
	if (index_end < length) {
		for (int i = index_end; i < length; ++i) {
			int p = *pointer_to_end_position;
			if (p != previous_position) {
				++index;
				*pointer_to_return_position = p;
				*pointer_to_return_value = std::max(pile_up * scaling_factor, baseline_value);

				++pointer_to_return_position;
				++pointer_to_return_value;
				previous_position = p;
			}
			--pile_up;
			++pointer_to_end_position;
		}
	}
	return_position.resize(index);
	return_value.resize(index);
	return qMakePair(return_position, return_value);
};

double MacsCallPeakWorker::get_p_score(const std::pair<int, double>& lo) {
	auto iter = this->p_score_map_.find(lo);
	if (iter != this->p_score_map_.end()) {
		return iter->second;
	}
	else {
		double val = -1.0 * poisson_cdf<false, true>(lo.first, lo.second);
		this->p_score_map_[lo] = val;
		return val;
	}
}

void MacsCallPeakWorker::calculate_p_q_map() {

	G_TASK_LOG("[MACS] Calculating p-q table...");

	double treat_scale = 1.0;
	std::unordered_map<double, double> p_score_stat;
	for (auto iter = locations_.begin(); iter != locations_.end(); ++iter) {
		QString chromosome = iter->first;
		int previous_position = 0;
		auto [treat_position, treat_value] = calculate_distribution(chromosome, this->d_, treat_scale, 0, true, this->end_shift_);
		auto [control_position, control_value] = calculate_distribution(chromosome, this->l_region_, double(this->d_) / this->l_region_, this->lambda_background_, false);
		auto& [position, treat, control] = make_pair(chromosome, treat_position, treat_value, control_position, control_value);

		treat_position.clear();
		treat_value.clear();
		control_position.clear();
		control_value.clear();

		int* pointer_to_position = position.data();
		double* ptr_treat = treat.data(), * ptr_control = control.data();

		previous_position = 0;
		qsizetype size = position.size();
		for (qsizetype i = 0; i < size; ++i) {
			double this_value = get_p_score(std::make_pair((int)(*ptr_treat), *ptr_control));
			double this_length = *pointer_to_position - previous_position;
			auto iter = p_score_stat.find(this_value);

			if (iter != p_score_stat.end()) {
				iter->second += this_length;
			}
			else {
				p_score_stat[this_value] = this_length;
			}
			previous_position = *pointer_to_position;
			++pointer_to_position;
			++ptr_treat;
			++ptr_control;
		}

	}

	std::size_t N = 0;
	for (auto iter = p_score_stat.cbegin(); iter != p_score_stat.cend(); ++iter) {
		N += iter->second;
	}
	std::size_t k = 1;
	double f = -log10(N);
	double previous_q_value = std::numeric_limits<int32_t>::max();

	QVector<double> unique_values;
	for (auto iter = p_score_stat.cbegin(); iter != p_score_stat.cend(); ++iter) {
		unique_values << iter->first;
	}

	unique_values = _Cs unique(unique_values);
	std::ranges::reverse(unique_values);

	qsizetype size = unique_values.size(), i;
	for (i = 0; i < size; ++i) {
		double v = unique_values[i];
		int l = p_score_stat[v];
		double q = v + std::log10(k) + f;
		if (q > previous_q_value) {
			q = previous_q_value;
		}
		if (q <= 0) {
			q = 0;
			break;
		}
		this->p_q_map_[v] = q;
		previous_q_value = q;
		k += l;
	}
	for (qsizetype j = i; j < size; ++j) {
		double v = unique_values[j];
		this->p_q_map_[v] = 0;
	}
};

std::tuple<QVector<int>, QVector<double>, QVector<double>>& MacsCallPeakWorker::make_pair(
	const QString& chromosome, 
	const QVector<int>& treat_position,
	const QVector<double>& treat_value, 
	const QVector<int>& control_position, 
	const QVector<double>& control_value
) {
	qsizetype length_treat = treat_position.size(), length_control = control_position.size();
	qsizetype maximum_length = length_treat + length_control;

	auto& [position, treat, control] = this->chromosome_position_treat_control_[chromosome];

	position.resize(maximum_length, 0);
	treat.resize(maximum_length, 0);
	control.resize(maximum_length, 0);

	const int* pointer_to_treat_position = treat_position.data(), * pointer_to_control_position = control_position.data();
	int* pointer_to_position = position.data();
	const double* pointer_to_treat_value = treat_value.data(), * pointer_to_control_value = control_value.data();
	double* ptr_treat = treat.data(), * ptr_control = control.data();
	qsizetype previous_position = 0, index_return = 0, index_treat = 0, index_control = 0;

	while (index_treat < length_treat && index_control < length_control) {
		*ptr_treat = *pointer_to_treat_value;
		*ptr_control = *pointer_to_control_value;

		if (*pointer_to_treat_position < *pointer_to_control_position) {
			*pointer_to_position = *pointer_to_treat_position;
			previous_position = *pointer_to_treat_position;

			++index_treat;
			++pointer_to_treat_position;
			++pointer_to_treat_value;
		}
		else if (*pointer_to_treat_position > *ptr_control) {
			*pointer_to_position = *pointer_to_control_position;
			previous_position = *pointer_to_control_position;

			++index_control;
			++pointer_to_control_position;
			++pointer_to_control_value;
		}
		else {
			*pointer_to_position = *pointer_to_treat_position;
			previous_position = *pointer_to_treat_position;

			++index_treat;
			++index_control;

			++pointer_to_control_position;
			++pointer_to_control_value;
			++pointer_to_treat_position;
			++pointer_to_treat_value;
		}
		++index_return;
		++pointer_to_position;
		++ptr_treat;
		++ptr_control;
	}
	position.resize(index_return);
	treat.resize(index_return);
	control.resize(index_return);
	return chromosome_position_treat_control_[chromosome];
};

void MacsCallPeakWorker::filter_duplicates() {
	qsizetype fragment_total = 0, fragment_filtered = 0;
	for (auto iter = locations_.begin(); iter != locations_.end(); ++iter) {
		for (auto& strand : iter->second) {
			fragment_total += strand.size();
			strand = _Cs unique(strand);
			fragment_filtered += strand.size();
		}
	}
	double proportion = ((double)fragment_filtered * 100 / fragment_total);
	G_TASK_LOG(QString::number(fragment_filtered) + " reads (" + QString::number(proportion, 'f', 2) + "%) remained after total " +
		QString::number(fragment_total) + " reads were filtered.");

	this->total_ = fragment_filtered;
};

bool MacsCallPeakWorker::detect_tag_size() {
	if (!this->bed_file_names_.isEmpty()) {
		return this->detect_tag_size_in_file();
	}
	else if (!this->fragments_objects_.isEmpty()) {
		return detect_tag_size_in_object();
	}
	return false;
}

bool MacsCallPeakWorker::detect_tag_size_in_object() {
	const auto& fragments = *this->fragments_objects_[0];
	for (const auto& [name, sequence] : fragments.data_) {
		if (_Cs sum(_Cs sapply(sequence, [](const auto& pair) { return pair.first.size(); })) > 100) {
			int tag_size = 0, count = 0;
			for (const auto& cell_fragments : sequence) {
				for (int j = 0; j < cell_fragments.first.size() && count < 100; ++j, ++count) {
					tag_size += cell_fragments.second[j] - cell_fragments.first[j];
				}
				if (count == 100) {
					this->tag_size_ = tag_size / 100;
					G_TASK_LOG("Estimated tag size : " + QString::number(this->tag_size_));
					return true;
				}
			}
		}
	}
	return false;
}

bool MacsCallPeakWorker::detect_tag_size_in_file() {
	auto bed_file = this->bed_file_names_[0].toUtf8();

	int buffer_length = 256;
	char buffer[256];

	gzFile fragments_file = gzopen(bed_file.data(), "rb");

	if (fragments_file == NULL) {
		return false;
	}

	bool end_of_file = false;
	while (!(end_of_file = gzgets(fragments_file, buffer, buffer_length) == 0))
	{
		if (buffer[0] != '#')break;
	}
	if (end_of_file) {
		gzclose(fragments_file);
		return false;
	}

	int length = 0;
	for (int i = 0; i < 100; ++i) {
		end_of_file = gzgets(fragments_file, buffer, buffer_length) == 0;
		if (end_of_file)break;
		const char* c = buffer;
		while (*c != '\t') {
			++c;
		}
		int start = _Cs atoi_specialized(&c);
		int end = _Cs atoi_specialized(&c);
		length += end - start;
	}
	gzclose(fragments_file);
	if (end_of_file) {
		return false;
	}
	this->tag_size_ = length / 100;
	G_TASK_LOG("Estimated tag size : " + QString::number(this->tag_size_));
	return true;
};

void MacsCallPeakWorker::parse_fragments_object() {

	for (const auto& fragments : this->fragments_objects_) {

		for (const auto& [name, data] : fragments->data_) {

			if (!this->locations_.contains(name)) {
				this->locations_[name] = QVector<QVector<int>>(2);
			}

			auto& location = this->locations_[name];

			const int n_cell = data.size();

			for (int i = 0; i < n_cell; ++i) {
				location[0] << _Cs cast<QVector>(data[i].first);
			}
		}
	}
};

bool MacsCallPeakWorker::parse_bed_file() {
	for (const auto& bed_file : this->bed_file_names_) {
		gzFile fragments_file = gzopen(bed_file.toUtf8().data(), "rb");

		if (fragments_file == NULL) {
			return false;
		}

		int buffer_length = 256;
		char buffer[256];

		bool not_end_of_file = false;
		while (not_end_of_file = gzgets(fragments_file, buffer, buffer_length) != 0)
		{
			if (buffer[0] != '#')break;
		}
		if (!not_end_of_file) {
			gzclose(fragments_file);
			return false;
		}

		const char* lead, * lag;
		do {
			lead = lag = buffer;
			while (*lead != '\t') {
				++lead;
			}
			QString chromosome_name = _Cs standardize_chromosome_name(QString::fromUtf8(lag, lead - lag));
			int start = _Cs atoi_specialized(&lead);
			++lead;
			int count = 2;
			while (count < 6) {

				while (*lead != '\t' && *lead != '\n') {
					++lead;
				}
				++count;

				if (*lead == '\n') {
					break;
				}
				++lead;
			}
			if (!this->locations_.contains(chromosome_name)) {
				this->locations_[chromosome_name] = QVector<QVector<int>>(2);
			}

			if (count < 6) {
				this->locations_[chromosome_name][0] << start;
			}
			else {
				if (*lead != '\n') {
					--lead;
				}
				--lead;
				if (*lead == '-') {
					this->locations_[chromosome_name][1] << start;
				}
				else {
					this->locations_[chromosome_name][0] << start;
				}
			}
		} while (gzgets(fragments_file, buffer, buffer_length) != 0);
		gzclose(fragments_file);
	}
	return true;
};