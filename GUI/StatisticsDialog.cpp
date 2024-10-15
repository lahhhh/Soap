#include "StatisticsDialog.h"

#include "MessageDialog.h"

#include "SoapGUI.h"

static std::tuple<double, double, double, double, double, double, double, double>  
get_statistics(const QVector<double>& data) {
	auto sorted = custom::sorted(data);
	double mean = custom::mean(sorted);
	
	const std::size_t size = sorted.size();
	double max = sorted.last(), min = sorted.first();
	double p25 = sorted[static_cast<std::size_t>(0.25 * size)];
	double p50 = sorted[static_cast<std::size_t>(0.50 * size)];
	double p75 = sorted[static_cast<std::size_t>(0.75 * size)];

	double var = 0;
	for (auto iter = sorted.cbegin(); iter != sorted.cend(); ++iter) {
		var += (*iter - mean) * (*iter - mean);
	}
	var /= (size - 1);
	double sd = std::sqrt(var);

	return std::make_tuple(mean, max, min, p25, p50, p75, var, sd);
};

StatisticsDialog::StatisticsDialog(BulkRna* bulk_rna) :
	handler_(bulk_rna)
{
	this->set_layout();
};

StatisticsDialog::StatisticsDialog(SingleCellRna* single_cell_rna) :
	handler_(single_cell_rna)
{
	this->set_layout();
};

StatisticsDialog::StatisticsDialog(SingleCellAtac* single_cell_atac) :
	handler_(single_cell_atac)
{
	this->set_layout();
};

StatisticsDialog::StatisticsDialog(SingleCellMultiome* single_cell_multiome) :
	handler_(single_cell_multiome)
{
	this->set_layout();
};

void StatisticsDialog::get_response(BulkRna* bulk_rna) {
	StatisticsDialog dlg(bulk_rna);
};

void StatisticsDialog::get_response(SingleCellRna* single_cell_rna) {
	StatisticsDialog dlg(single_cell_rna);
};

void StatisticsDialog::get_response(SingleCellAtac* single_cell_atac) {
	StatisticsDialog dlg(single_cell_atac);
};

void StatisticsDialog::get_response(SingleCellMultiome* single_cell_multiome) {
	StatisticsDialog dlg(single_cell_multiome);
};

void StatisticsDialog::refresh_interface(const std::tuple<double, double, double, double, double, double, double, double>& data) {
	auto [mean, max, min, p25, p50, p75, var, sd] = data;

	this->mean_view_label_->setText(QString::number(mean, 'f', 2));
	this->mean_view_label_->adjustSize();

	this->max_view_label_->setText(QString::number(max, 'f', 2));
	this->max_view_label_->adjustSize();

	this->min_view_label_->setText(QString::number(min, 'f', 2));
	this->min_view_label_->adjustSize();

	this->p25_view_label_->setText(QString::number(p25, 'f', 2));
	this->p25_view_label_->adjustSize();

	this->p50_view_label_->setText(QString::number(p50, 'f', 2));
	this->p50_view_label_->adjustSize();

	this->p75_view_label_->setText(QString::number(p75, 'f', 2));
	this->p75_view_label_->adjustSize();

	this->var_view_label_->setText(QString::number(var, 'f', 2));
	this->var_view_label_->adjustSize();

	this->sd_view_label_->setText(QString::number(sd, 'f', 2));
	this->sd_view_label_->adjustSize();
};

void StatisticsDialog::s_refresh() {

	QString feature = this->feature_lineedit_->text();
	bool normalized = this->normalize_switch_->value_;
	auto res = this->handler_.get_data(QUERY_INFO{ feature, normalized, false });

	if (!res.is_valid()) {
		MessageDialog::get_response("Error", "No Data Found.");
		return;
	}

	if (res.is_continuous()) {
		this->refresh_interface(get_statistics(res.get_continuous()));
	}
	else {
		MessageDialog::get_response("Error", "Data Type is not continuous.");
	}
};

void StatisticsDialog::set_layout() {
	this->main_layout_ = new QGridLayout;

	int row = 0;

	G_SET_LABEL(this->feature_label_, "Feature", soap::MiddleSize);
	this->main_layout_->addWidget(this->feature_label_, row, 0);

	G_SET_EMPTY_LINEEDIT(this->feature_lineedit_, soap::MiddleSize);
	G_SET_BUTTON(this->refresh_button_, "Refresh", soap::MiddleSize);
	QHBoxLayout* row_layout = new QHBoxLayout;
	row_layout->addWidget(this->feature_lineedit_);
	row_layout->addWidget(this->refresh_button_);
	this->main_layout_->addLayout(row_layout, row, 1);
	connect(this->refresh_button_, &QPushButton::clicked, this, &StatisticsDialog::s_refresh);

	++row;

	G_SET_LABEL(this->normalize_label_, "Normalized", soap::MiddleSize);
	G_SET_SWITCH(this->normalize_switch_, false, this->normalize_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->normalize_label_, row, 0);
	this->main_layout_->addWidget(this->normalize_switch_, row, 1);

	++row;

	G_SET_LABEL(this->mean_label_, "Mean Value", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->mean_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->mean_label_, row, 0);
	this->main_layout_->addWidget(this->mean_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->max_label_, "Max Value", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->max_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->max_label_, row, 0);
	this->main_layout_->addWidget(this->max_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->min_label_, "Min Value", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->min_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->min_label_, row, 0);
	this->main_layout_->addWidget(this->min_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->p25_label_, "25% Value", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->p25_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->p25_label_, row, 0);
	this->main_layout_->addWidget(this->p25_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->p50_label_, "50%(median) Value", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->p50_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->p50_label_, row, 0);
	this->main_layout_->addWidget(this->p50_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->p75_label_, "75% Value", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->p75_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->p75_label_, row, 0);
	this->main_layout_->addWidget(this->p75_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->var_label_, "Variance", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->var_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->var_label_, row, 0);
	this->main_layout_->addWidget(this->var_view_label_, row, 1);

	++row;

	G_SET_LABEL(this->sd_label_, "Standard Deviation", soap::MiddleSize);
	G_SET_EMPTY_LABEL(this->sd_view_label_, soap::MiddleSize);
	this->main_layout_->addWidget(this->sd_label_, row, 0);
	this->main_layout_->addWidget(this->sd_view_label_, row, 1);

	this->setLayout(this->main_layout_);

	auto feature_names = this->handler_.get_feature_names();
	feature_names.numeric_names;

	this->feature_lineedit_->setCompleter(new QCompleter(feature_names.numeric_names, this));

	G_SET_ICON;
	this->setWindowTitle("Statistics");
	this->exec();
};