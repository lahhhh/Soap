#include "CustomPlotUtility.h"

#include "ViolinPlot.h"

namespace custom_plot {

	namespace utility {

		QCPRange get_range(
			double min_val, 
			double max_val,
			double margin_proportion) 
		{

			double val_span = max_val - min_val;
			double val_min_limit = min_val - margin_proportion * val_span;
			double val_maval_limit = max_val + margin_proportion * val_span;

			return QCPRange(val_min_limit, val_maval_limit);
		};

		Eigen::ArrayXd catmull(const Eigen::ArrayXd& x, int n_interpolation) {

			int len_x = x.size();
			int len_smoothed = (len_x - 3) * (n_interpolation + 1) + 1;

			Eigen::ArrayXd smoothed(len_smoothed);
			
			Eigen::ArrayXd u1 = Eigen::ArrayXd::LinSpaced(n_interpolation + 2, 0, 1).segment(0, n_interpolation + 1);
			Eigen::ArrayXd u2 = u1 * u1, u3 = u1 * u1 * u1;

			double x1, x2, x3, x4;
			for (int i = 0; i < len_x - 3; ++i) {

				x1 = -x[i] + 3 * x[i + 1] - 3 * x[i + 2] + x[i + 3];
				x2 = 2 * x[i] - 5 * x[i + 1] + 4 * x[i + 2] - x[i + 3];
				x3 = -x[i] + x[i + 2];
				x4 = 2 * x[i + 1];

				for (int j = 0; j <= n_interpolation; ++j) {
					smoothed[i * (n_interpolation + 1) + j] = 0.5 * (x1 * u3[j] + x2 * u2[j] + x3 * u1[j] + x4);
				}
			}
			smoothed[(len_x - 3) * (n_interpolation + 1)] = x[len_x - 2];

			return smoothed;
		};

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > violin_smooth(
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			double x_center,
			int n_interpolation
		) {
			auto [left_x, left_y] = left_violin_smooth(x, y, x_center, n_interpolation);
			Eigen::ArrayXd reverse_x = (2 * x_center - left_x).reverse();
			Eigen::ArrayXd reverse_y = left_y.reverse();

			return std::make_pair(custom::concatenated(left_x, reverse_x), custom::concatenated(left_y, reverse_y));
		};

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > left_violin_smooth(
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			double x_center,
			int n_interpolation
		) {
			const std::size_t original_length = x.size();

			Eigen::ArrayXd new_x = Eigen::ArrayXd::Constant(original_length + 4, x_center), new_y(original_length + 4);

			new_y[1] = 2 * y[0] - y[1];
			new_y[0] = 3 * y[0] - 2 * y[1];
			new_y[original_length + 2] = 2 * y[original_length - 1] - y[original_length - 2];
			new_y[original_length + 3] = 3 * y[original_length - 1] - 2 * y[original_length - 2];
			new_x.segment(2, original_length) = x;
			new_y.segment(2, original_length) = y;

			Eigen::ArrayXd smooth_x = utility::catmull(new_x, n_interpolation);
			Eigen::ArrayXd smooth_y = utility::catmull(new_y, n_interpolation);

			smooth_x = custom::set_upper_bound(smooth_x, x_center);
			return std::make_pair(smooth_x, smooth_y);
		};

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > right_violin_smooth(
			const Eigen::ArrayXd& x,
			const Eigen::ArrayXd& y,
			double x_center,
			int n_interpolation
		) {
			auto [left_x, left_y] = left_violin_smooth(x, y, x_center, n_interpolation);

			Eigen::ArrayXd reverse_x = (2 * x_center - left_x);

			return std::make_pair(reverse_x, left_y);
		}

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > violin_curve(
			const Eigen::ArrayXd& data, double center, int n, double zero_space) {

			auto [min, max] = std::ranges::minmax(data);

			Eigen::ArrayXd p;

			if (min == max) {
				p = Eigen::ArrayXd::LinSpaced(n, min - zero_space, max + zero_space);
			}
			else {
				p = Eigen::ArrayXd::LinSpaced(n, min, max);
			}
			auto kde = evaluate_KDE(data, p);
			Eigen::ArrayXd left = (center - ((kde / kde.maxCoeff()) * 0.8 ));
			return custom_plot::utility::violin_smooth(left, p, center);
		}

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > left_violin_curve(
			const Eigen::ArrayXd& data, double center, int n, double zero_space) {

			auto [min, max] = std::ranges::minmax(data);

			Eigen::ArrayXd p;

			if (min == max) {
				p = Eigen::ArrayXd::LinSpaced(n, min - zero_space, max + zero_space);
			}
			else {
				p = Eigen::ArrayXd::LinSpaced(n, min, max);
			}

			auto kde = evaluate_KDE(data, p);

			Eigen::ArrayXd left = (center - ((kde / kde.maxCoeff()) * 0.8));

			return custom_plot::utility::left_violin_smooth(left, p, center);
		}

		std::pair<Eigen::ArrayXd, Eigen::ArrayXd > right_violin_curve(
			const Eigen::ArrayXd& data, double center, int n, double zero_space) {

			auto [min, max] = std::ranges::minmax(data);

			Eigen::ArrayXd p;

			if (min == max) {
				p = Eigen::ArrayXd::LinSpaced(n, min - zero_space, max + zero_space);
			}
			else {
				p = Eigen::ArrayXd::LinSpaced(n, min, max);
			}

			auto kde = evaluate_KDE(data, p);

			Eigen::ArrayXd left = (center - ((kde / kde.maxCoeff()) * 0.8));

			return custom_plot::utility::right_violin_smooth(left, p, center);
		}

		std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXd> histogram(
			const Eigen::ArrayXd& arr, int bins, double zero_space) {

			double amin = arr.minCoeff(), amax = arr.maxCoeff();

			double d = (amax - amin) / bins;

			if (d == 0) {
				Eigen::ArrayXi counts = Eigen::ArrayXi::Zero(bins);

				counts[bins / 2] = arr.size();

				Eigen::ArrayXd edges = Eigen::ArrayXd::LinSpaced(bins + 1, amin - bins * zero_space, amin + bins * zero_space);

				Eigen::ArrayXd centers = Eigen::ArrayXd::LinSpaced(bins, amin - (bins - 1) * zero_space, amin + (bins - 1) * zero_space);

				return std::make_tuple(counts, edges, centers);
			}
			else {
				Eigen::ArrayXd edges = Eigen::ArrayXd::LinSpaced(bins + 1, amin, amax);

				Eigen::ArrayXi counts = Eigen::ArrayXi::Zero(bins + 1);

				for (int i = 0; i < arr.size(); ++i) {
					++counts[floor((arr[i] - amin) / d)];
				}

				counts[bins - 1] += counts[bins];

				Eigen::ArrayXd centers = (edges.segment(0, bins) + edges.segment(1, bins)) / 2;

				return std::make_tuple(counts.segment(0, bins), edges, centers);
			}
		};


		QColor float_to_rgb(double r, double g, double b, int alpha) {

			return QColor((int)(r * 255) % 256, (int)(g * 255) % 256, (int)(b * 255) % 256, alpha);
		};

		QColor hsv_to_rgb(double h, double s, double v, int alpha) {

			double t;
			double f = modf(h * 6.0, &t);
			int i = ((int)t) % 6;
			double p = v * (1 - s);
			double q = v * (1 - s * f);
			t = v * (1 - (s * (1 - f)));

			switch (i) {
			case 0:	return float_to_rgb(v, t, p, alpha);
			case 1:	return float_to_rgb(q, v, p, alpha);
			case 2:	return float_to_rgb(p, v, t, alpha);
			case 3:	return float_to_rgb(p, q, v, alpha);
			case 4:	return float_to_rgb(t, p, v, alpha);
			case 5:	return float_to_rgb(v, p, q, alpha);
			default:
				return QColor(0, 0, 0, alpha);
			}
		};

		Eigen::MatrixXd rgb_to_lab(const Eigen::MatrixXd& rgb)
		{
			const int n_color = rgb.rows();
			Eigen::ArrayXd R = rgb.col(0), G = rgb.col(1), B = rgb.col(2), X, Y, Z, FX(n_color), FY(n_color), FZ(n_color), l(n_color), a, b;
			X = 0.412453 * R + 0.357580 * G + 0.180423 * B;
			Y = 0.212671 * R + 0.715160 * G + 0.072169 * B;
			Z = 0.019334 * R + 0.119193 * G + 0.950227 * B;

			X /= 0.95047;
			Y /= 1.0;
			Z /= 1.08883;

			for (int i = 0; i < n_color; ++i) {
				FX[i] = X[i] > 0.008856 ? std::pow(X[i], 1.0 / 3) : (7.787 * X[i] + 0.137931);
				FY[i] = Y[i] > 0.008856 ? std::pow(Y[i], 1.0 / 3) : (7.787 * Y[i] + 0.137931);
				FZ[i] = Z[i] > 0.008856 ? std::pow(Z[i], 1.0 / 3) : (7.787 * Z[i] + 0.137931);
				l[i] = Y[i] > 0.008856 ? (116.0 * FY[i] - 16.0) : (903.3 * Y[i]);
			}

			a = 500 * (FX - FY);
			b = 200 * (FY - FZ);
			Eigen::MatrixXd ret(n_color, 3);
			ret << l, a, b;
			return ret;
		}


		Eigen::MatrixXi lab_to_rgb(const Eigen::MatrixXd& lab)
		{
			const int n_color = lab.rows();

			Eigen::ArrayXd
				l = lab.col(0),
				a = lab.col(1),
				b = lab.col(2),
				X(n_color),
				Y(n_color),
				Z(n_color),
				fx(n_color),
				fy(n_color),
				fz(n_color),
				R, G, B;

			fy = pow((l + 16) / 116, 3);

			for (int i = 0; i < n_color; ++i) {

				if (fy[i] < 0.008856)fy[i] = l[i] / 903.3;

				Y[i] = fy[i];

				if (fy[i] > 0.08856)fy[i] = std::pow(fy[i], 1.0 / 3);
				else fy[i] = 7.787 * fy[i] + 16.0 / 116;

				fx[i] = a[i] / 500 + fy[i];

				if (fx[i] > 0.206893)X[i] = std::pow(fx[i], 3);
				else X[i] = (fx[i] - 16.0 / 116) / 7.787;

				fz[i] = fy[i] - b[i] / 200;

				if (fz[i] > 0.206893)Z[i] = std::pow(fz[i], 3);
				else Z[i] = (fz[i] - 16.0 / 116) / 7.787;
			}

			X = X * (0.950456 * 255);

			Y = Y * 255;

			Z = Z * (1.088754 * 255);

			R = 3.240479 * X - 1.53715 * Y - 0.498535 * Z;

			G = -0.969256 * X + 1.875992 * Y + 0.041556 * Z;

			B = 0.055648 * X - 0.204043 * Y + 1.057311 * Z;

			Eigen::ArrayXi ir, ig, ib;

			ir = R.cast<int>();
			ig = G.cast<int>();
			ib = B.cast<int>();

			for (int i = 0; i < n_color; ++i) {
				ir[i] = ir[i] > 255 ? 255 : (ir[i] < 0 ? 0 : ir[i]);
				ig[i] = ig[i] > 255 ? 255 : (ig[i] < 0 ? 0 : ig[i]);
				ib[i] = ib[i] > 255 ? 255 : (ib[i] < 0 ? 0 : ib[i]);
			}

			Eigen::MatrixXi ret(n_color, 3);
			ret << ir, ig, ib;
			return ret;
		}

		Eigen::MatrixXi lab_to_srgb(const Eigen::MatrixXd& lab)
		{
			const int n_color = lab.rows();
			Eigen::ArrayXd
				l = lab.col(0),
				a = lab.col(1),
				b = lab.col(2),
				X(n_color),
				Y(n_color),
				Z(n_color),
				fx(n_color),
				fy(n_color),
				fz(n_color),
				R, G, B;

			fy = pow((l + 16) / 116, 3);

			for (int i = 0; i < n_color; ++i) {
				if (fy[i] < 0.008856)fy[i] = l[i] / 903.3;

				Y[i] = fy[i];

				if (fy[i] > 0.08856)fy[i] = std::pow(fy[i], 1.0 / 3);
				else fy[i] = 7.787 * fy[i] + 16.0 / 116;

				fx[i] = a[i] / 500 + fy[i];

				if (fx[i] > 0.206893)X[i] = std::pow(fx[i], 3);
				else X[i] = (fx[i] - 16.0 / 116) / 7.787;

				fz[i] = fy[i] - b[i] / 200;

				if (fz[i] > 0.206893)Z[i] = std::pow(fz[i], 3);
				else Z[i] = (fz[i] - 16.0 / 116) / 7.787;
			}

			X *= 95.0456;

			Y *= 100;

			Z *= 108.8754;

			R = 0.032406 * X - 0.015371 * Y - 0.0049895 * Z;
			G = -0.0096891 * X + 0.018757 * Y + 0.00041914 * Z;
			B = 0.00055708 * X - 0.0020401 * Y + 0.01057 * Z;

			for (int i = 0; i < n_color; ++i) {
				if (R[i] <= 0.00313)R[i] = R[i] * 12.92;
				else R[i] = exp(log(R[i]) / 2.4) * 1.055 - 0.055;

				if (G[i] <= 0.00313)G[i] = G[i] * 12.92;
				else G[i] = exp(log(G[i]) / 2.4) * 1.055 - 0.055;

				if (B[i] <= 0.00313)B[i] = B[i] * 12.92;
				else B[i] = exp(log(B[i]) / 2.4) * 1.055 - 0.055;
			}

			Eigen::ArrayXi ir, ig, ib;

			ir = round(R * 255).cast<int>();
			ig = round(G * 255).cast<int>();
			ib = round(B * 255).cast<int>();

			for (int i = 0; i < n_color; ++i) {
				ir[i] = ir[i] > 255 ? 255 : (ir[i] < 0 ? 0 : ir[i]);
				ig[i] = ig[i] > 255 ? 255 : (ig[i] < 0 ? 0 : ig[i]);
				ib[i] = ib[i] > 255 ? 255 : (ib[i] < 0 ? 0 : ib[i]);
			}

			Eigen::MatrixXi ret(n_color, 3);
			ret << ir, ig, ib;
			return ret;
		}

		QVector<QColor> rainbow(int n, int alpha) {

			if (n <= 1)
			{
				return QVector<QColor>{ QColor(0, 0, 0, alpha)};
			}

			QVector<QColor> res(n);

			for (int i = 0; i < n; ++i) {
				
				double h = i / (double)n;

				res[i] = hsv_to_rgb(h, 1.0, 1.0, alpha);
			}

			return res;
		};

		int get_max_text_width(const QStringList& text, const QFont& font) {

			QFontMetrics fm(font);

			int width = 0;
			for (const auto& str : text) {
				int w = fm.boundingRect(str).width();
				width = width < w ? w : width;
			}
			return width;
		};

		int get_max_text_height(const QStringList& text, const QFont& font) {

			QFontMetrics fm(font);

			int height = 0;
			for (const auto& str : text) {
				int w = fm.boundingRect(str).height();
				height = height < w ? w : height;
			}

			return height;
		};

		QVector<QColor> kmeans_palette(int n, int alpha, int m, int random_state) {

			QVector<QColor> ret(n);

			if (n < 1) {
			
				return ret;
			}
			else if (n == 1) {

				ret << Qt::black;

				return ret;
			}

			srand(random_state);

			Eigen::MatrixXd rgb_matrix = (Eigen::MatrixXd::Random(m, 3).array() + 1.0) / 2;

			auto [lab_matrix, _] = custom::kmeans_hartigan_wong_mt(rgb_to_lab(rgb_matrix).transpose(), n);
			Eigen::MatrixXi srgb_matrix = lab_to_srgb(lab_matrix.transpose());

			for (int i = 0; i < n; ++i) {
				ret[i] = QColor(srgb_matrix(i, 0), srgb_matrix(i, 1), srgb_matrix(i, 2), alpha);
			}
			return ret;
		};

		QColor gradient_color(double proportion, const QColor& low_color, const QColor& high_color) {

			return QColor(low_color.red() + (int)(proportion * (high_color.red() - low_color.red())),
				low_color.green() + (int)(proportion * (high_color.green() - low_color.green())),
				low_color.blue() + (int)(proportion * (high_color.blue() - low_color.blue())));
		};


		std::pair<double, double> calculate_text_size(
			const QString& text, 
			const QFont& font, 
			double xrange, 
			double yrange,
			double rect_width,
			double rect_height)
		{
			// Step 1: Calculate text size in pixels using QFontMetrics
			QFontMetrics fm(font);

			auto [w, h] = fm.boundingRect(text).size();

			// Calculate text size in plot range units
			return { w / rect_width * xrange, h / rect_height * yrange };
		}

		std::pair<QVector<double>, QVector<double>> circos_curve(
			double x_center,
			double y_center,
			double x_start,
			double y_start,
			double x_end,
			double y_end,
			int n_segment
		) {

			double x_control = x_center;
			double y_control = y_center;

			QVector<double> xs, ys;

			int n_point = n_segment + 1;
			for (int i = 0; i < n_point; ++i) {

				double t = static_cast<double>(i) / n_point;

				double x = (1 - t) * (1 - t) * x_start + 2 * (1 - t) * t * x_control + t * t * x_end;
				double y = (1 - t) * (1 - t) * y_start + 2 * (1 - t) * t * y_control + t * t * y_end;

				xs << x;
				ys << y;
			}

			return { xs, ys };
		};

		std::pair<QVector<double>, QVector<double>> arc(
			double x_center,
			double y_center,
			double radius,
			double angle_a,
			double angle_b,
			int n_segment
		) {

			auto angles = custom::linspaced(n_segment + 1, angle_a, angle_b);

			auto xs = custom::sapply(angles, [x_center, radius](auto angle) {
				return x_center + radius * std::cos(angle);
			});

			auto ys = custom::sapply(angles, [y_center, radius](auto angle) {
				return y_center + radius * std::sin(angle);
			});

			return { xs, ys };
		};
	};
};