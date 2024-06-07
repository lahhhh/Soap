#include "CustomDebug.h"

#include <QDebug>

namespace custom {

	void print(const Eigen::MatrixXd& mat, int max_row) {
		
		qDebug() << "-----------------| print start |-----------------";

		const int nrow = mat.rows(), ncol = mat.cols();
		if (max_row == -1 || max_row > nrow) {
			max_row = nrow;
		}
		qDebug() << "Eigen Dense Double Matrix | nrow : " + QString::number(nrow) + " | ncol : " + QString::number(ncol);

		for (int i = 0; i < max_row; ++i) {
			QString row_str = "Row " + QString::number(i) + " |";
			for (int j = 0; j < ncol; ++j) {
				row_str += "    " + QString::number(mat(i, j));
			}
			qDebug() << row_str;
		}

		qDebug() << "-----------------| print end |-----------------";
	}

	void print(const Eigen::ArrayXd& arr) {
		qDebug() << "-----------------| print start |-----------------";

		std::size_t size = arr.size();
		qDebug() << "Eigen Dense Double Matrix | length : " + QString::number(size);

		QString row_str = "Array : ";
		for (int i = 0; i < size; ++i) {
			row_str += "    " + QString::number(arr[i]);
		}
		qDebug() << row_str;

		qDebug() << "-----------------| print end |-----------------";
	};

}