#include "RnaNormalizeWorker.h"

#include <QFile>

void RnaNormalizeWorker::run() {

	if (this->method_ == "FPKM") {
		this->fpkm();
	}
	if (this->method_ == "TPM") {
		this->tpm();
	}
	if (this->method_ == "Normalize1") {
		this->normalize_1();
	}

	G_TASK_END;
};

std::pair<QStringList, QVector<int>> RnaNormalizeWorker::get_transcript_length() {

	QFile file(FILE_HUMAN_TRANSCRIPT_LENGTH);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		return {};
	}

	QStringList gene_names;
	QVector<int> transcript_lengths;

	QTextStream in(&file);

	QString line = in.readLine();

	while (!line.isNull()) {

		auto s = line.split('\t');

		if (s.size() != 2) {
			return {};
		}

		QString gene_name = s[0];
		int length = s[1].toInt();

		gene_names << gene_name;
		transcript_lengths << length;

		line = in.readLine();
	}

	return std::make_pair(gene_names, transcript_lengths);
};

bool RnaNormalizeWorker::normalize_1() {

	DenseDouble normalized = this->counts_->to_double();

	Eigen::ArrayXd depth = normalized.mat_.colwise().sum();
	double f = _Cs median(depth);

	depth = f / depth;

	normalized.mat_.array().rowwise() *= depth.transpose();
	normalized.mat_ = (log1p(normalized.mat_.array())).eval();

	emit x_normalize_ready(normalized);

	return true;
};

bool RnaNormalizeWorker::fpkm() {

	auto [gene_names, transcript_lengths] = this->get_transcript_length();

	if (gene_names.isEmpty()) {
		G_TASK_WARN("Database Loading Faild");
		return false;
	}

	auto index = _Cs which(_Cs in(this->counts_->rownames_, gene_names));

	if (index.size() < (this->counts_->rows() / 2) || index.size() < 10000) {
		G_TASK_WARN("Too many genes not found in database.");
		return false;
	}

	DenseDouble normalized = this->counts_->to_double();
	normalized.row_reorder(index);
	transcript_lengths = _Cs reordered(transcript_lengths, _Cs index_of(normalized.rownames_, gene_names));

	Eigen::ArrayXd total_reads = normalized.mat_.colwise().sum();
	normalized.mat_.array().colwise() /= (_Cs cast<Eigen::ArrayX>(transcript_lengths)).cast<double>();
	normalized.mat_.array().rowwise() /= total_reads.transpose();
	normalized.mat_.array() *= 1e9;

	emit x_normalize_ready(normalized);

	return true;
};

bool RnaNormalizeWorker::tpm() {

	auto [gene_names, transcript_lengths] = this->get_transcript_length();

	if (gene_names.isEmpty()) {
		G_TASK_WARN("Database Loading Faild");
		return false;
	}

	auto index = _Cs which(_Cs in(this->counts_->rownames_, gene_names));

	if (index.size() < (this->counts_->rows() / 2) || index.size() < 10000) {
		G_TASK_WARN("Too many genes not found in database.");
		return false;
	}

	DenseDouble normalized = this->counts_->to_double();
	normalized.row_reorder(index);
	transcript_lengths = _Cs reordered(transcript_lengths, _Cs index_of(normalized.rownames_, gene_names));

	normalized.mat_.array().colwise() /= (_Cs cast<Eigen::ArrayX>(transcript_lengths)).cast<double>();
	Eigen::ArrayXd total_reads = normalized.mat_.colwise().sum();
	normalized.mat_.array().rowwise() /= total_reads.transpose();
	normalized.mat_.array() *= 1e6;

	emit x_normalize_ready(normalized);

	return true;
};