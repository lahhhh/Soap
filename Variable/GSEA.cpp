#include "Gsea.h"
#include "Custom.h"

GSEA::GSEA(
	const QString& database,
	const QStringList& comparison, 
	const QMap<QString, double>& es, 
	const QMap<QString, double>& nes, 
	const QMap<QString, double>& p, 
	const QMap<QString, double>& fdr,
	const QMap<QString, double>& fwer, 
	const QMap<QString, int>& size, 
	const QVector<double>& correlations,
	const QMap<QString, QVector<int> >& gene_location, 
	const QMap<QString, QVector<double>>& point_x, 
	const QMap<QString, QVector<double>>& point_y
) : 
	database_(database),
	comparison_(comparison), 
	pathways_(es.keys()),
	correlations_(correlations), 
	gene_location_(gene_location), 
	point_x_(point_x), 
	point_y_(point_y)
{
	int length = correlations.size();
	int pos = ceil(correlations[0] * 10), neg = floor(correlations[length - 1] * 10);
	int length2 = pos - neg + 1;

	QVector<QColor> positive_color{ 
		QColor(255, 192, 229, 221), 
		QColor(255, 137, 137, 221), 
		QColor(255, 112, 128, 221),
		QColor(255, 90, 90, 221), 
		QColor(239, 64, 64, 221), 
		QColor(195, 22, 11, 221) 
	};

	QVector<QColor> negative_color{
		QColor(213, 213, 255, 221), 
		QColor(199, 193, 255, 221), 
		QColor(138, 138, 255, 221),
		QColor(107, 88, 239, 221),
		QColor(39, 0, 209, 221), 
		QColor(69, 0, 173, 221) 
	};

	
	this->bar_points_ = QVector<double>(length2);
	this->bar_points_[0] = 0;
	this->bar_points_[length2 - 1] = length - 1;
	for (int i = 0; i < length2 - 1; ++i) {
		double loc = (pos - i) * 0.1;
		for (int j = 0; j < length - 1; ++j) {
			if (this->correlations_[j] >= loc && this->correlations_[j + 1] < loc) {
				this->bar_points_[i] = j;
				break;
			}
		}
	}
	for (int i = 0; i < length2 - 1; ++i) {
		double color_degree = nearbyint(this->correlations_[this->bar_points_[i]] * 10);
		if (color_degree > 0) {
			this->bar_colors_ << positive_color[std::max(std::min((int)color_degree - 1, 5), 0)];
		}
		else {
			this->bar_colors_ << negative_color[std::max(std::min((int)(-color_degree), 5), 0)];
		}
	}
	this->mat_.set_rownames(this->pathways_);
	QVector<double> enrichment_score, normalized_enrichment_score, p_value, family_wise_error_rate, false_discovery_rate;
	QVector<int> pathway_size;
	for (const auto& path : this->pathways_) {
		enrichment_score << es[path];
		normalized_enrichment_score << nes[path];
		p_value << p[path];
		pathway_size << size[path];
		family_wise_error_rate << fwer[path];
		false_discovery_rate << fdr[path];
	}
	this->mat_.update(METADATA_GSEA_ENRICHMENT_SCORE, enrichment_score);
	this->mat_.update(METADATA_GSEA_NORMALIZED_ENRICHMENT_SCORE, normalized_enrichment_score);
	this->mat_.update(METADATA_GSEA_P_VALUE, p_value);
	this->mat_.update(METADATA_GSEA_PATHWAY_SIZE, pathway_size);
	this->mat_.update(METADATA_GSEA_FALSE_DISCOVERY_RATE, false_discovery_rate);
	this->mat_.update(METADATA_GSEA_FAMILY_WISE_ERROR_RATE, family_wise_error_rate);
	this->mat_.row_reorder(_Cs order(this->mat_.get_const_double_reference(METADATA_GSEA_P_VALUE)));
}

bool GSEA::is_empty() const {
	return this->pathways_.isEmpty();
};

QVector<double> GSEA::get_nes(const QStringList& pathways, const QString& filter_type, double threshold) const {
	int size = pathways.size();
	auto index = _Cs index_of(pathways, this->mat_.rownames_);
	const QVector<double>& nes = this->mat_.get_const_double_reference(METADATA_GSEA_NORMALIZED_ENRICHMENT_SCORE);
	QVector<double> ret(size, 0);
	for (int i = 0; i < size; ++i) {
		if (index[i] != -1) {
			ret[i] = nes[index[i]];
		}
	}
	QVector<double>* filter;
	if (filter_type == "P") {
		filter = (QVector<double>*)this->mat_.at(METADATA_GSEA_P_VALUE);
		for (int i = 0; i < size; ++i) {
			if (index[i] != -1) {
				if (filter->at(index[i]) > threshold) {
					ret[i] = 0;
				}
			}
		}
	}
	else if (filter_type == "FDR") {
		filter = (QVector<double>*)this->mat_.at(METADATA_GSEA_FALSE_DISCOVERY_RATE);
		for (int i = 0; i < size; ++i) {
			if (index[i] != -1) {
				if (filter->at(index[i]) > threshold) {
					ret[i] = 0;
				}
			}
		}
	}
	else if (filter_type == "FWER") {
		filter = (QVector<double>*)this->mat_.at(METADATA_GSEA_FAMILY_WISE_ERROR_RATE);
		for (int i = 0; i < size; ++i) {
			if (index[i] != -1) {
				if (filter->at(index[i]) > threshold) {
					ret[i] = 0;
				}
			}
		}
	}
	return ret;
};