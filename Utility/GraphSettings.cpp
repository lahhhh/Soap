#include "GraphSettings.h"

#include "CustomPlot.h"

GraphSettings::GraphSettings() {

	this->string_info_[TITLE] = "";
	this->bool_info_[IS_TITLE_ACTIVE] = false;

	this->font_info_[TITLE_FONT] = DEFAULT_TITLE_FONT;
	this->bool_info_[IS_TITLE_FONT_ACTIVE] = false;

	this->integer_info_[AXIS_STYLE] = 0;
	
	this->bool_info_[IS_AXIS_STYLE_ACTIVE] = false;

	this->string_info_[BOTTOM_TITLE] = "";
	this->bool_info_[IS_BOTTOM_TITLE_ACTIVE] = false;

	this->font_info_[BOTTOM_TITLE_FONT] = DEFAULT_BOTTOM_TITLE_FONT;
	this->bool_info_[IS_BOTTOM_TITLE_FONT_ACTIVE] = false;

	this->font_info_[BOTTOM_LABEL_FONT] = DEFAULT_BOTTOM_LABEL_FONT;
	this->bool_info_[IS_BOTTOM_LABEL_FONT_ACTIVE] = false;

	this->integer_info_[BOTTOM_LABEL_ANGLE] = 0;
	this->bool_info_[IS_BOTTOM_LABEL_ANGLE_ACTIVE] = false;

	this->string_info_[LEFT_TITLE] = "";
	this->bool_info_[IS_LEFT_TITLE_ACTIVE] = false;

	this->font_info_[LEFT_TITLE_FONT] = DEFAULT_LEFT_TITLE_FONT;
	this->bool_info_[IS_LEFT_TITLE_FONT_ACTIVE] = false;

	this->font_info_[LEFT_LABEL_FONT] = DEFAULT_LEFT_LABEL_FONT;
	this->bool_info_[IS_LEFT_LABEL_FONT_ACTIVE] = false;

	this->integer_info_[LEFT_LABEL_ANGLE] = 0;
	this->bool_info_[IS_LEFT_LABEL_ANGLE_ACTIVE] = false;

	this->string_list_info_[LEGEND_TITLES] = {""};
	this->bool_info_[IS_LEGEND_TITLES_ACTIVE] = false;

	this->font_info_[LEGEND_TITLE_FONT] = DEFAULT_LEGEND_TITLE_FONT;
	this->bool_info_[IS_LEGEND_TITLE_FONT_ACTIVE] = false;

	this->font_info_[LEGEND_LABEL_FONT] = DEFAULT_LEGEND_LABEL_FONT;
	this->bool_info_[IS_LEGEND_LABEL_FONT_ACTIVE] = false;

	this->font_info_[LEGEND_TICK_LABEL_FONT] = DEFAULT_LEGEND_TICK_LABEL_FONT;
	this->bool_info_[IS_LEGEND_TICK_LABEL_FONT_ACTIVE] = false;

	this->integer_info_[LEGEND_COLUMN_WIDTH] = 0;
	this->bool_info_[IS_LEGEND_COLUMN_WIDTH_ACTIVE] = false;

	this->integer_info_[LEGEND_ROW_WIDTH] = 0;
	this->bool_info_[IS_LEGEND_ROW_WIDTH_ACTIVE] = false;

	this->bool_info_[IS_LEGEND_TICK_SHOWN] = false;
	this->bool_info_[IS_LEGEND_TICK_SHOWN_ACTIVE] = false;

	this->color_info_[GRADIENT_LOW_COLOR] = DEFAULT_GRADIENT_LOW_COLOR;
	this->bool_info_[IS_GRADIENT_LOW_COLOR_ACTIVE] = false;

	this->color_info_[GRADIENT_MIDDLE_COLOR] = DEFAULT_GRADIENT_MIDDLE_COLOR;
	this->bool_info_[IS_GRADIENT_MIDDLE_COLOR_ACTIVE] = false;

	this->color_info_[GRADIENT_HIGH_COLOR] = DEFAULT_GRADIENT_HIGH_COLOR;
	this->bool_info_[IS_GRADIENT_HIGH_COLOR_ACTIVE] = false;

	this->bool_info_[IS_SCATTER_LABELED] = false;
	this->bool_info_[IS_SCATTER_LABELED_ACTIVE] = false;

	this->integer_info_[SCATTER_POINT_SIZE] = DEFAULT_SCATTER_POINT_SIZE;
	this->bool_info_[IS_SCATTER_POINT_SIZE_ACTIVE] = false;

	this->font_info_[SCATTER_LABEL_FONT] = DEFAULT_SCATTER_LABEL_FONT;
	this->bool_info_[IS_SCATTER_LABEL_FONT_ACTIVE] = false;

	this->bool_info_[IS_PALATTE_ACTIVE] = false;

	this->bool_info_[IS_TRANSPARENT_BACKGROUND] = false;
	this->bool_info_[IS_TRANSPARENT_BACKGROUND_ACTIVE] = false;

	this->bool_info_[USE_FACET_VIOLIN_PLOT] = true;
	this->bool_info_[USE_FACET_VIOLIN_PLOT_ACTIVE] = false;

	this->bool_info_[USE_BOXPLOT] = false;
	this->bool_info_[USE_BOXPLOT_ACTIVE] = false;

	this->bool_info_[BOXPLOT_DRAW_OUTLIER] = false;
	this->bool_info_[BOXPLOT_DRAW_OUTLIER_ACTIVE] = false;
};

void GraphSettings::apply_palette(const QStringList& names, QVector<QColor>& colors) const {
	if (!this->bool_info_[IS_PALATTE_ACTIVE] || !this->active())return;

	const qsizetype size = names.size();

	auto end = this->palette_.cend();
	for (qsizetype i = 0; i < size; ++i) {
		auto iter = this->palette_.find(names[i]);
		if (iter != end) {
			colors[i] = *iter;
		}
	}
}

void GraphSettings::apply_palette(const QString& name, QColor& color) const {

	if (!this->bool_info_[IS_PALATTE_ACTIVE] || !this->active())return;

	auto iter = this->palette_.find(name);
	if (iter != this->palette_.cend()) {
		color = *iter;
	}
};

QList<QColor> GraphSettings::palette(const QStringList& levels) const {

	auto ulevels = custom::unique(levels);

	auto colors = custom_plot::utility::kmeans_palette(ulevels.size());

	this->apply_palette(ulevels, colors);

	return custom::sapply(levels, [&ulevels, &colors](auto&& level) {return colors[ulevels.indexOf(level)]; });
};

QMap<QString, QColor> GraphSettings::palette_map(const QStringList& levels) const {

	auto ulevels = custom::unique(levels);

	auto colors = custom_plot::utility::kmeans_palette(ulevels.size());

	this->apply_palette(ulevels, colors);

	QMap<QString, QColor> pal;

	int n_level = ulevels.size();

	for (int i = 0; i < n_level; ++i) {
		pal[ulevels[i]] = colors[i];
	}

	return pal;
};


bool GraphSettings::active() const {
	return this->bool_info_[SETTINGS_ACTIVATED];
}

QString GraphSettings::get_title(const QString& title) const {
	return this->active() && this->bool_info_[IS_TITLE_ACTIVE] ? this->string_info_[TITLE] : title;
};

QString GraphSettings::get_legend_title(const QString& title, unsigned int index) const {
	if (!this->active() || !this->bool_info_[IS_LEGEND_TITLES_ACTIVE]) {
		return title;
	}
	QStringList legend_titles = this->string_list_info_[LEGEND_TITLES];
	if (legend_titles.size() <= index) {
		return title;
	}
	else {
		return legend_titles[index];
	}
};

QFont GraphSettings::get_title_font(const QFont& title_font) const {
	return this->active() && this->bool_info_[IS_TITLE_FONT_ACTIVE] ?
		this->font_info_[TITLE_FONT] : title_font;
};

GraphSettings::AxisStyle GraphSettings::get_axis_style(AxisStyle style) const {
	return this->active() && this->bool_info_[IS_AXIS_STYLE_ACTIVE] ?
		(AxisStyle)this->integer_info_[AXIS_STYLE] : style;
};

QFont GraphSettings::get_bottom_label_font(const QFont& bottom_label_font) const {
	return this->active() && this->bool_info_[IS_BOTTOM_LABEL_FONT_ACTIVE] ?
		this->font_info_[BOTTOM_LABEL_FONT] : bottom_label_font;
};

QFont GraphSettings::get_left_label_font(const QFont& left_label_font) const {
	return this->active() && this->bool_info_[IS_LEFT_LABEL_FONT_ACTIVE] ?
		this->font_info_[LEFT_LABEL_FONT] : left_label_font;
};

QFont GraphSettings::get_label_font(QCPAxis::AxisType type) const {
	if (type == QCPAxis::atLeft || QCPAxis::atRight) {
		return this->get_left_label_font();
	}
	else {
		return this->get_bottom_label_font();
	}
};

int GraphSettings::get_bottom_label_angle(const int bottom_label_angle) const {
	return this->active() && this->bool_info_[IS_BOTTOM_LABEL_ANGLE_ACTIVE] ?
		this->integer_info_[BOTTOM_LABEL_ANGLE] : bottom_label_angle;
};

int GraphSettings::get_left_label_angle(const int left_label_angle) const {
	return this->active() && this->bool_info_[IS_LEFT_LABEL_ANGLE_ACTIVE] ?
		this->integer_info_[LEFT_LABEL_ANGLE] : left_label_angle;
};

int GraphSettings::get_label_angle(QCPAxis::AxisType type) const {
	if (type == QCPAxis::atLeft || type == QCPAxis::atRight) {
		return this->get_left_label_angle();
	}
	else {
		return this->get_bottom_label_angle();
	}
};

QFont GraphSettings::get_bottom_title_font(const QFont& bottom_title_font) const {
	return this->active() && this->bool_info_[IS_BOTTOM_TITLE_FONT_ACTIVE] ?
		this->font_info_[BOTTOM_TITLE_FONT] : bottom_title_font;
};

QFont GraphSettings::get_left_title_font(const QFont& left_title_font) const {
	return this->active() && this->bool_info_[IS_LEFT_TITLE_FONT_ACTIVE] ?
		this->font_info_[LEFT_TITLE_FONT] : left_title_font;
};

QColor GraphSettings::get_gradient_low_color(const QColor& color) const {
	return this->active() && this->bool_info_[IS_GRADIENT_LOW_COLOR_ACTIVE] ?
		this->color_info_[GRADIENT_LOW_COLOR] : color;
};

QColor GraphSettings::get_gradient_middle_color(const QColor& color) const {
	return this->active() && this->bool_info_[IS_GRADIENT_MIDDLE_COLOR_ACTIVE] ?
		this->color_info_[GRADIENT_MIDDLE_COLOR] : color;
};

QColor GraphSettings::get_gradient_high_color(const QColor& color) const {
	return this->active() && this->bool_info_[IS_GRADIENT_HIGH_COLOR_ACTIVE] ?
		this->color_info_[GRADIENT_HIGH_COLOR] : color;
};

QFont GraphSettings::get_legend_title_font(const QFont& legend_title_font) const {
	return this->active() && this->bool_info_[IS_LEGEND_TITLE_FONT_ACTIVE] ?
		this->font_info_[LEGEND_TITLE_FONT] : legend_title_font;
};

QFont GraphSettings::get_legend_label_font(const QFont& legend_label_font) const {
	return this->active() && this->bool_info_[IS_LEGEND_LABEL_FONT_ACTIVE] ?
		this->font_info_[LEGEND_LABEL_FONT] : legend_label_font;
};

int GraphSettings::get_legend_column_width(const int legend_column_width) const {
	return this->active() && this->bool_info_[IS_LEGEND_COLUMN_WIDTH_ACTIVE] ?
		this->integer_info_[LEGEND_COLUMN_WIDTH] : legend_column_width;
};

int GraphSettings::get_legend_row_width(const int legend_row_width) const {
	return this->active() && this->bool_info_[IS_LEGEND_ROW_WIDTH_ACTIVE] ?
		this->integer_info_[LEGEND_ROW_WIDTH] : legend_row_width;
};

int GraphSettings::get_scatter_point_size(int scatter_size) const {
	return this->active() && this->bool_info_[IS_SCATTER_POINT_SIZE_ACTIVE] ?
		this->integer_info_[SCATTER_POINT_SIZE] : scatter_size;
};

bool GraphSettings::is_scatter_labeled() const {
	return this->active() && this->bool_info_[IS_SCATTER_LABELED_ACTIVE] ?
		this->bool_info_[IS_SCATTER_LABELED] : DEFAULT_SCATTER_LABELED;
};

QFont GraphSettings::get_scatter_label_font() const {
	return this->active() && this->bool_info_[IS_SCATTER_LABEL_FONT_ACTIVE] ?
		this->font_info_[SCATTER_LABEL_FONT] : DEFAULT_SCATTER_LABEL_FONT;
};

bool GraphSettings::use_boxplot() const {
	return this->active() && this->bool_info_[USE_BOXPLOT_ACTIVE] ?
		this->bool_info_[USE_BOXPLOT] : DEFAULT_USE_BOXPLOT;
};

bool GraphSettings::boxplot_draw_outlier() const {
	return this->active() && this->bool_info_[BOXPLOT_DRAW_OUTLIER_ACTIVE] ?
		this->bool_info_[BOXPLOT_DRAW_OUTLIER] : DEFAULT_BOXPLOT_DRAW_OUTLIER;
};

bool GraphSettings::use_facet_violin_plot() const {
	return this->active() && this->bool_info_[USE_FACET_VIOLIN_PLOT_ACTIVE] ?
		this->bool_info_[USE_FACET_VIOLIN_PLOT] : DEFAULT_USE_FACET_VIOLIN_PLOT;
};

bool GraphSettings::is_legend_tick_shown() const {
	return this->active() && this->bool_info_[IS_LEGEND_TICK_SHOWN_ACTIVE] ?
		this->bool_info_[IS_LEGEND_TICK_SHOWN] : DEFAULT_LEGEND_TICK_SHOWN;
};

bool GraphSettings::is_transparent_background() const {
	return this->active() && this->bool_info_[IS_TRANSPARENT_BACKGROUND_ACTIVE] ?
		this->bool_info_[IS_TRANSPARENT_BACKGROUND] : DEFAULT_TRANSPARENT_BACKGROUND;
}

void GraphSettings::set_activated(bool activated) {
	this->bool_info_[SETTINGS_ACTIVATED] = activated;
};

QFont GraphSettings::get_legend_tick_label_font(const QFont& legend_tick_label_font) const {
	return this->active() && this->bool_info_[IS_LEGEND_TICK_LABEL_FONT_ACTIVE] ?
		this->font_info_[LEGEND_TICK_LABEL_FONT] : legend_tick_label_font;
};

QString GraphSettings::get_bottom_title(const QString& title) const {
	return this->active() && this->bool_info_[IS_BOTTOM_TITLE_ACTIVE] ? this->string_info_[BOTTOM_TITLE] : title;
};

QString GraphSettings::get_left_title(const QString& title) const {
	return this->active() && this->bool_info_[IS_LEFT_TITLE_ACTIVE] ? this->string_info_[LEFT_TITLE] : title;
};