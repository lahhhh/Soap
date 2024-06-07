#pragma once

#include "Identifier.h"

#include <QDialog>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>
#include <QHBoxLayout>

#include "Switch.h"

#include "SingleCellRna.h"
#include "SingleCellMultiome.h"

#include "FeatureHandler.h"

static std::tuple<double, double, double, double, double, double, double, double> 
get_statistics(const QVector<double>& data);

class StatisticsDialog 
	: public QDialog
{
	Q_OBJECT

public:

	StatisticsDialog(BulkRna* bulk_rna);
	StatisticsDialog(SingleCellRna* single_cell_rna);
	StatisticsDialog(SingleCellAtac* single_cell_atac);
	StatisticsDialog(SingleCellMultiome* single_cell_multiome);

	static void get_response(BulkRna* bulk_rna);
	static void get_response(SingleCellRna* single_cell_rna);
	static void get_response(SingleCellAtac* single_cell_atac);
	static void get_response(SingleCellMultiome* single_cell_multiome);

	FeatureHandler handler_;

	QGridLayout* main_layout_ = nullptr;

	QHBoxLayout* feature_layout_ = nullptr;

	QLabel* feature_label_ = nullptr;
	QLineEdit* feature_lineedit_ = nullptr;
	QPushButton* refresh_button_ = nullptr;

	QLabel* normalize_label_ = nullptr;
	Switch* normalize_switch_ = nullptr;

	QLabel* mean_label_ = nullptr;
	QLabel* mean_view_label_ = nullptr;

	QLabel* max_label_ = nullptr;
	QLabel* max_view_label_ = nullptr;

	QLabel* min_label_ = nullptr;
	QLabel* min_view_label_ = nullptr;

	QLabel* p25_label_ = nullptr;
	QLabel* p25_view_label_ = nullptr;

	QLabel* p50_label_ = nullptr;
	QLabel* p50_view_label_ = nullptr;

	QLabel* p75_label_ = nullptr;
	QLabel* p75_view_label_ = nullptr;

	QLabel* var_label_ = nullptr;
	QLabel* var_view_label_ = nullptr;

	QLabel* sd_label_ = nullptr;
	QLabel* sd_view_label_ = nullptr;

private:

	void set_layout();

	void refresh_interface(const std::tuple<double, double, double, double, double, double, double, double>& data);

private slots:

	void s_refresh();

};

