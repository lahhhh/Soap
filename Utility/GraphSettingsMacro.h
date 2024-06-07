#pragma once


/* string */

#define TITLE "title"

#define LEGEND_TITLES "legend titles"

#define BOTTOM_TITLE "bottom title"
#define LEFT_TITLE "left title"

/* font */

#define TITLE_FONT "title font"

#define LEGEND_TITLE_FONT "legend title font"
#define LEGEND_LABEL_FONT "legend label font"

#define BOTTOM_TITLE_FONT "bottom title font"
#define LEFT_TITLE_FONT "left title font"

#define BOTTOM_LABEL_FONT "bottom label font"
#define LEFT_LABEL_FONT "left label font"

#define LEGEND_TICK_LABEL_FONT "legend tick label font"

#define SCATTER_LABEL_FONT "scatter label font"

/* integer */

#define BOTTOM_LABEL_ANGLE "bottom label angle"
#define LEFT_LABEL_ANGLE "left label angle"

#define LEGEND_COLUMN_WIDTH "legend column width"
#define LEGEND_ROW_WIDTH "legend row width"

#define SCATTER_POINT_SIZE "scatter point size"

#define AXIS_STYLE "axis style"

/* bool */

#define SETTINGS_ACTIVATED "settings activated"

#define IS_TITLE_ACTIVE "is title active"
#define IS_TITLE_FONT_ACTIVE "is title font active"

#define IS_AXIS_STYLE_ACTIVE "is axis style active"

#define IS_LEGEND_TITLES_ACTIVE "is legend titles active"

#define IS_LEGEND_TITLE_FONT_ACTIVE "is legend title font active"
#define IS_LEGEND_LABEL_FONT_ACTIVE "is legend label font active"

#define IS_BOTTOM_TITLE_ACTIVE "is bottom title active"
#define IS_LEFT_TITLE_ACTIVE "is left title active"

#define IS_BOTTOM_TITLE_FONT_ACTIVE "is bottom title font active"
#define IS_LEFT_TITLE_FONT_ACTIVE "is left title font active"

#define IS_BOTTOM_LABEL_FONT_ACTIVE "is bottom label font active"
#define IS_LEFT_LABEL_FONT_ACTIVE "is left label font active"

#define IS_BOTTOM_LABEL_ANGLE_ACTIVE "is bottom label angle active"
#define IS_LEFT_LABEL_ANGLE_ACTIVE "is left label angle active"

#define IS_LEGEND_TICK_LABEL_FONT_ACTIVE "is legend tick label font active"

#define IS_GRADIENT_HIGH_COLOR_ACTIVE "is gradient high color active"
#define IS_GRADIENT_MIDDLE_COLOR_ACTIVE "is gradient middle color active"
#define IS_GRADIENT_LOW_COLOR_ACTIVE "is gradient low color active"

#define IS_LEGEND_COLUMN_WIDTH_ACTIVE "is legend column width active"
#define IS_LEGEND_ROW_WIDTH_ACTIVE "is legend row width active"

#define IS_SCATTER_LABELED "is scatter labeled"
#define IS_SCATTER_LABELED_ACTIVE "is scatter labeled active"

#define IS_LEGEND_TICK_SHOWN "is legend tick shown"
#define IS_LEGEND_TICK_SHOWN_ACTIVE "is legend tick shown active"

#define IS_SCATTER_POINT_SIZE_ACTIVE "is scatter point size active"

#define IS_SCATTER_LABEL_FONT_ACTIVE "is scatter label font active"

#define IS_PALATTE_ACTIVE "is palatte active"

#define IS_TRANSPARENT_BACKGROUND "is transparent background"
#define IS_TRANSPARENT_BACKGROUND_ACTIVE "is transparent background active"

#define USE_FACET_VIOLIN_PLOT "use facet violin plot"
#define USE_FACET_VIOLIN_PLOT_ACTIVE "use facet violin plot active"

#define USE_BOXPLOT "use boxplot"
#define USE_BOXPLOT_ACTIVE "use boxplot active"

#define BOXPLOT_DRAW_OUTLIER "boxplot draw outlier"
#define BOXPLOT_DRAW_OUTLIER_ACTIVE "boxplot draw outlier active"


/* color */

#define GRADIENT_HIGH_COLOR "gradient high color"
#define GRADIENT_MIDDLE_COLOR "gradient middle color"
#define GRADIENT_LOW_COLOR "gradient low color"

/* defaults */

#define DEFAULT_TITLE_FONT QFont("Arial", 20, QFont::Bold)

#define DEFAULT_BOTTOM_TITLE_FONT QFont("Arial", 18, QFont::Bold)
#define DEFAULT_LEFT_TITLE_FONT QFont("Arial", 18, QFont::Bold)

#define DEFAULT_BOTTOM_LABEL_FONT QFont("Arial", 15)
#define DEFAULT_LEFT_LABEL_FONT QFont("Arial", 15)

#define DEFAULT_LEGEND_TITLE_FONT QFont("Arial", 18)
#define DEFAULT_LEGEND_LABEL_FONT QFont("Arial", 15)

#define DEFAULT_SCATTER_LABEL_FONT QFont("Arial", 15)
#define DEFAULT_LEGEND_TICK_LABEL_FONT QFont("Arial", 10)

#define DEFAULT_AXIS_STYLE AxisStyle::Simple

#define DEFAULT_BOTTOM_LABEL_ANGLE 0
#define DEFAULT_LEFT_LABEL_ANGLE 0

#define DEFAULT_GRADIENT_LOW_COLOR QColor(221, 221, 221)
#define DEFAULT_GRADIENT_MIDDLE_COLOR QColor(110, 110, 238)
#define DEFAULT_GRADIENT_HIGH_COLOR QColor(0, 0, 255)

#define DEFAULT_LEGEND_COLUMN_WIDTH 0
#define DEFAULT_LEGEND_ROW_WIDTH 30

#define DEFAULT_SCATTER_POINT_SIZE 2

#define DEFAULT_SCATTER_LABELED false
#define DEFAULT_LEGEND_TICK_SHOWN false
#define DEFAULT_TRANSPARENT_BACKGROUND false
#define DEFAULT_USE_FACET_VIOLIN_PLOT true
#define DEFAULT_USE_BOXPLOT false
#define DEFAULT_BOXPLOT_DRAW_OUTLIER false