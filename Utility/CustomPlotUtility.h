#pragma once

#include "Identifier.h"
#include <qCustomPlot.h>
#include "SoapTextElement.h"

#include "Custom.h"

namespace custom_plot {

	namespace utility {

		template <typename ContainerType>
			requires custom::is_specific_container<ContainerType, double>
		QCPRange get_range(const ContainerType& vec, double margin_proportion = 0.1) {

			auto [min_val, max_val] = ::std::ranges::minmax(vec);
			double val_span = max_val - min_val;
			double val_min_limit = min_val - margin_proportion * val_span, val_maval_limit = max_val + margin_proportion * val_span;

			return QCPRange(val_min_limit, val_maval_limit);
		}

		template <typename ContainerType1, typename ContainerType2>
			requires custom::is_specific_container<ContainerType1, double>&& custom::is_specific_container<ContainerType2, double>
		std::pair<QCPRange, QCPRange> get_range(const ContainerType1& x, const ContainerType2& y, double margin_proportion = 0.1) {

			return std::make_pair(get_range(x, margin_proportion), get_range(y, margin_proportion));
		}

		QCPRange get_range(double min_val, double max_val, double margin_proportion = 0.1);

		// return [1:-2]
		Eigen::ArrayXd catmull(const Eigen::ArrayXd& x, int n_interpolation);

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > violin_smooth(
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			double center_x,
			int n_interpolation = 32
		);

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > left_violin_smooth(
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			double center_x,
			int n_interpolation = 32
		);

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > right_violin_smooth(
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			double center_x,
			int n_interpolation = 32
		);

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > violin_curve(
			const Eigen::ArrayXd& orig, double center = 0.0, int bins = 16, double zero_space = 0.01);

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > left_violin_curve(
			const Eigen::ArrayXd& orig, double center = 0.0, int bins = 16, double zero_space = 0.01);

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > right_violin_curve(
			const Eigen::ArrayXd& orig, double center = 0.0, int bins = 16, double zero_space = 0.01);

		std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXd> histogram(
			const Eigen::ArrayXd& orig, int bins = 16, double zero_space = 0.01);

		QColor float_to_rgb(double r, double g, double b, int alpha);

		QColor hsv_to_rgb(double h, double s, double v, int alpha);

		// skip the first step. RGB -> [0,1]
		Eigen::MatrixXd rgb_to_lab(const Eigen::MatrixXd& rgb);

		Eigen::MatrixXi lab_to_rgb(const Eigen::MatrixXd& lab);

		Eigen::MatrixXi lab_to_srgb(const Eigen::MatrixXd& lab);

		QVector<QColor> rainbow(int n, int alpha = 255);

		int get_max_text_width(const QStringList& text, const QFont& font);
		int get_max_text_height(const QStringList& text, const QFont& font);

		QVector<QColor> kmeans_palette(int n, int alpha = 255, int m = 2000, int random_state = 1997);

		QColor gradient_color(double proportion, const QColor& low_color, const QColor& high_color);

		std::pair<double, double> calculate_text_size(
			const QString& text,
			const QFont& font,
			double xrange,
			double yrange,
			double rect_width,
			double rect_height);
	};
};

