#pragma once

#include "Identifier.h"

#include <QColor>

#include "CustomMatrix.h"

class GSEA
{
public:

	enum class DataType : int { Plain = 0 };

	G_CLASS_FUNCTION_DEFAULT(GSEA);

	GSEA(
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
	);


	bool is_empty() const;

	DataType data_type_{ DataType::Plain };

	CustomMatrix mat_;

	QString database_;
	QStringList comparison_;
	QStringList pathways_;

	QVector<double> correlations_;

	QVector<double> bar_points_;
	QVector<QColor> bar_colors_;

	QMap<QString, QVector<int>> gene_location_;

	QMap<QString, QVector<double>> point_x_;
	QMap<QString, QVector<double>> point_y_;

	QVector<double> get_nes(const QStringList& pathways_, const QString& filter_type, double threshold) const;

	G_SET_IDENTIFIER("GSEA");

};
