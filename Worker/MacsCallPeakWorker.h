#pragma once

#include "Identifier.h"

#include <unordered_map>
#include <QVector>
#include <QHash>

#include "SingleCellMultiome.h"
#include "GenomicRange.h"

namespace std
{

	template<> struct hash< std::pair<int, double> >
	{
		std::size_t operator()(std::pair<int, double> const& p) const
		{
			constexpr size_t _FNV_offset_basis = 14695981039346656037ULL;
			constexpr size_t _FNV_prime = 1099511628211ULL;

			size_t _Val = _FNV_offset_basis;

			const unsigned char* const int_begin = &(reinterpret_cast<const unsigned char&>(p.first));
			const unsigned char* const double_begin = &(reinterpret_cast<const unsigned char&>(p.second));

			constexpr std::size_t int_length = sizeof(int);
			constexpr std::size_t double_length = sizeof(double);

			for (size_t _Idx = 0; _Idx < int_length; ++_Idx) {
				_Val ^= static_cast<size_t>(int_begin[_Idx]);
				_Val *= _FNV_prime;
			}

			for (size_t _Idx = 0; _Idx < double_length; ++_Idx) {
				_Val ^= static_cast<size_t>(double_begin[_Idx]);
				_Val *= _FNV_prime;
			}

			return _Val;
		}
	};

} // namespace std

/*
* modified from MACS software package
* Zhang, Yong et al. ¡°Model-based analysis of ChIP-Seq (MACS).¡± 
* Genome biology vol. 9,9 (2008): R137. doi:10.1186/gb-2008-9-9-r137
*/

class MacsCallPeakWorker 
	: public QObject
{
	Q_OBJECT

public:
	MacsCallPeakWorker(const QList<const Fragments*>& fragments_objects) :fragments_objects_(fragments_objects)	{};

	void run();

	static GenomicRange call_peak(const QList<const Fragments*>& fragments_objects);

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_genomic_range_ready(GenomicRange);

private:

	QStringList bed_file_names_;

	QList<const Fragments*> fragments_objects_;

	long long genome_size_ = 2.7e9;

	int band_width_ = 300;
	int ext_size_ = 200;
	int end_shift_ = -100;
	int s_region_ = 1000;
	int l_region_ = 10000;
	int minimum_length_ = 200;
	int maximum_gap_ = 200;
	int tag_size_ = 200;
	int total_ = 0;
	int d_ = 200;

	double lambda_background_ = 0.;
	double log_q_value_ = -log10(0.05);

	RunLengthEncoding<QString> seqnames_;
	IRange ranges_;
	RunLengthEncoding<char> strand_;
	QStringList name_;

	QVector<double> fold_change_;
	QVector<double> p_value_;
	QVector<double> q_value_;

	QVector<int> summit_position_;
	QVector<int> score_;

	// strand_ direction ==> * : 0 , + : 0 , - : 1
	std::unordered_map<QString, QVector<QVector<int>>> locations_;

	std::unordered_map<QString, std::tuple<QVector<int>, QVector<double>, QVector<double>>> chromosome_position_treat_control_;

	std::unordered_map<std::pair<int, double>, double> p_score_map_;

	QVector<double> calculate_q_score(const QVector<double>& treat, const QVector<double>& control);

	std::unordered_map<double, double> p_q_map_;

	QPair<QVector<int>, QVector<double>> calculate_distribution(const QString& chromosome, int d, double scaling_factor,
		double baseline_value, bool directional, int end_shift = 0);

	bool parse_bed_file();

	void parse_fragments_object();

	void filter_duplicates();

	bool detect_tag_size();

	bool detect_tag_size_in_file();

	bool detect_tag_size_in_object();

	void detect_peaks();

	void call_peak_without_control();

	void chromosome_call_peak(const QString& chromosome);

	std::tuple<QVector<int>, QVector<double>, QVector<double>>& make_pair(const QString& chromosome, const QVector<int>& treat_position,
		const QVector<double>& treat_value, const QVector<int>& control_position, const QVector<double>& control_value);

	double get_p_score(const std::pair<int, double>& lo);

	void calculate_p_q_map();

	void close_peak_without_subpeaks(const QString& chromosome, const QVector<double>& score);
};


