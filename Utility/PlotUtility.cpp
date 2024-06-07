#include "PlotUtility.h"

QRect get_align_rect(int width, int height, Qt::Alignment alignment) {

	QRect ret(0, 0, width, height);

	if (alignment & Qt::AlignLeft) {
		ret.moveLeft(0);
	}

	if (alignment & Qt::AlignRight) {
		ret.moveRight(0);
	}

	if (alignment & Qt::AlignHCenter) {
		ret.moveLeft(-width / 2);
	}

	if (alignment & Qt::AlignVCenter) {
		ret.moveLeft(-height / 2);
	}

	if (alignment & Qt::AlignTop) {
		ret.moveTop(0);
	}

	if (alignment & Qt::AlignBottom) {
		ret.moveBottom(0);
	}

	return ret;
};