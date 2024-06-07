#include "IRange.h"

void IRange::append(int start, int width) {
	this->start_ << start;
	this->width_ << width;
}

qsizetype IRange::size() const {
	return this->start_.size();
};

void IRange::append(const IRange& range) {
	this->start_ << range.start_;
	this->width_ << range.width_;
};

void IRange::clear() {
	this->start_.clear();
	this->width_.clear();
}

IRange::IRange(const QVector<int>& start, const QVector<int>& width):
	start_(start),
	width_(width)
{

};