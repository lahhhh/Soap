#include "Switch.h"

Animator::Animator(QObject* target, QObject* parent) : QVariantAnimation(parent) {
	setTargetObject(target);
}

Animator::~Animator() {
	stop();
}

QObject* Animator::targetObject() const {
	return target.data();
}

void Animator::setTargetObject(QObject* _target) {
	if (target.data() == _target)
		return;
	target = _target;
}

void Animator::updateCurrentValue(const QVariant& value) {
	Q_UNUSED(value);

	if (!target.isNull()) {
		auto update = QEvent(QEvent::StyleAnimationUpdate);
		update.setAccepted(false);
		QCoreApplication::sendEvent(target.data(), &update);
		if (!update.isAccepted())
			stop();
	}
}

void Animator::updateState(QAbstractAnimation::State newState, QAbstractAnimation::State oldState) {
	QVariantAnimation::updateState(newState, oldState);
}

void Animator::setup(int duration, QEasingCurve easing) {
	setDuration(duration);
	setEasingCurve(easing);
}

void Animator::interpolate(const QVariant& _start, const QVariant& end) {
	setStartValue(_start);
	setEndValue(end);
	start();
}

void Animator::setCurrentValue(const QVariant& value) {
	setStartValue(value);
	setEndValue(value);
	updateCurrentValue(currentValue());
}



SelectionControl::SelectionControl(QWidget* parent) : QAbstractButton(parent) {
	setObjectName("SelectionControl");
	setCheckable(true);
}

void SelectionControl::enterEvent(QEnterEvent* e) {
	setCursor(Qt::PointingHandCursor);
	QAbstractButton::enterEvent(e);
}

Qt::CheckState SelectionControl::checkState() const {
	return isChecked() ? Qt::Checked : Qt::Unchecked;
}

void SelectionControl::checkStateSet() {
	const auto state = checkState();
	emit stateChanged(state);
	toggle(state);
}

void SelectionControl::nextCheckState() {
	QAbstractButton::nextCheckState();
	SelectionControl::checkStateSet();
}

void Switch::set_status(bool checked) {
	if (checked) {
		this->value_ = true;
		track_brush_animation_->setStartValue(colorFromOpacity(this->style_.track_on_brush_, this->style_.track_on_opacity_));
		track_brush_animation_->setEndValue(colorFromOpacity(this->style_.track_on_brush_, this->style_.track_on_opacity_));
		thumb_brush_animation_->setStartValue(colorFromOpacity(this->style_.thumb_on_brush_, this->style_.thumb_on_opacity_));
		thumb_brush_animation_->setEndValue(colorFromOpacity(this->style_.thumb_on_brush_, this->style_.thumb_on_opacity_));
		this->setChecked(true);
	}
	else {
		this->value_ = false;
		track_brush_animation_->setStartValue(colorFromOpacity(this->style_.track_off_brush_, this->style_.track_off_opacity_));
		track_brush_animation_->setEndValue(colorFromOpacity(this->style_.track_off_brush_, this->style_.track_off_opacity_));
		thumb_brush_animation_->setStartValue(colorFromOpacity(this->style_.thumb_off_brush_, this->style_.thumb_off_opacity_));
		thumb_brush_animation_->setEndValue(colorFromOpacity(this->style_.thumb_off_brush_, this->style_.thumb_off_opacity_));
		this->setChecked(false);
	}
};

void Switch::init() {
	/* setup animations */
	thumb_brush_animation_ = new Animator{ this, this };
	track_brush_animation_ = new Animator{ this, this };
	thumb_position_animation_ = new Animator{ this, this };
	thumb_position_animation_->setup(this->style_.thumb_position_animation_.duration, this->style_.thumb_position_animation_.easing);
	track_brush_animation_->setup(this->style_.track_brush_animation_.duration, this->style_.track_brush_animation_.easing);
	thumb_brush_animation_->setup(this->style_.thumb_brush_animation_.duration, this->style_.thumb_brush_animation_.easing);
	/* set init values */
	if (this->value_) {
		track_brush_animation_->setStartValue(colorFromOpacity(this->style_.track_on_brush_, this->style_.track_on_opacity_));
		track_brush_animation_->setEndValue(colorFromOpacity(this->style_.track_on_brush_, this->style_.track_on_opacity_));
		thumb_brush_animation_->setStartValue(colorFromOpacity(this->style_.thumb_on_brush_, this->style_.thumb_on_opacity_));
		thumb_brush_animation_->setEndValue(colorFromOpacity(this->style_.thumb_on_brush_, this->style_.thumb_on_opacity_));
		this->setChecked(true);
	}
	else {
		track_brush_animation_->setStartValue(colorFromOpacity(this->style_.track_off_brush_, this->style_.track_off_opacity_));
		track_brush_animation_->setEndValue(colorFromOpacity(this->style_.track_off_brush_, this->style_.track_off_opacity_));
		thumb_brush_animation_->setStartValue(colorFromOpacity(this->style_.thumb_off_brush_, this->style_.thumb_off_opacity_));
		thumb_brush_animation_->setEndValue(colorFromOpacity(this->style_.thumb_off_brush_, this->style_.thumb_off_opacity_));
		this->setChecked(false);
	}
	/* set standard palettes */
	setSizePolicy(QSizePolicy(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Fixed));
}

QRect Switch::indicatorRect() {
	const auto w = this->style_.indicator_margin_.left() + this->style_.indicator_width_ + this->style_.indicator_margin_.right();
	return QRect(0, 0, w, this->style_.height_);
}

Switch::Switch(bool default_value, QWidget* parent) : SelectionControl(parent), value_(default_value) {
	init();
}

Switch::Switch(
	bool default_value, 
	QPushButton* link_button,
	QWidget* parent
) : 
	SelectionControl(parent), 
	link_button_(link_button), 
	value_(default_value), 
	mode_(Switch::LinkMode::LinkToButton)
{
	if (link_button == nullptr) {
		this->mode_ = Switch::LinkMode::NoLink;
	}
	init();
}

Switch::Switch(
	bool default_value, 
	QLabel* link_label, 
	QWidget* parent
) : 
	SelectionControl(parent), 
	link_label_(link_label), 
	value_(default_value), 
	mode_(Switch::LinkMode::LinkToLabel)
{
	if (link_label == nullptr) {
		this->mode_ = Switch::LinkMode::NoLink;
	}
	init();
}

QSize Switch::sizeHint() const {
	auto h = this->style_.height_;
	auto w = this->style_.indicator_margin_.left() + this->style_.indicator_width_ + this->style_.indicator_margin_.right();

	return QSize(w, h);
}

void Switch::paintEvent(QPaintEvent*) {
	/* for desktop usage we do not need Radial reaction */

	QPainter p(this);

	const auto _indicatorRect = indicatorRect();
	auto trackMargin = this->style_.indicator_margin_;
	QRectF trackRect = _indicatorRect.marginsRemoved(trackMargin);

	p.setOpacity(1.0);
	p.setPen(Qt::NoPen);
	/* draw track */
	p.setBrush(track_brush_animation_->currentValue().value<QColor>());
	p.setRenderHint(QPainter::Antialiasing, true);
	p.drawRoundedRect(trackRect, CORNER_RADIUS, CORNER_RADIUS);
	p.setRenderHint(QPainter::Antialiasing, false);
	trackRect.setLeft(thumb_position_animation_->currentValue().toInt() - THUMB_WIDTH / 2);
	trackRect.setRight(thumb_position_animation_->currentValue().toInt() + THUMB_WIDTH / 2);
	trackRect.setBottom(trackRect.bottom() + 2);
	trackRect.setTop(trackRect.top() - 2);
	auto thumbRect = trackRect;

	p.setBrush(thumb_brush_animation_->currentValue().value<QColor>());
	p.setRenderHint(QPainter::Antialiasing, true);
	p.drawRoundedRect(thumbRect, THUMB_RADIUS, THUMB_RADIUS);
	p.setRenderHint(QPainter::Antialiasing, false);
}

void Switch::resizeEvent(QResizeEvent* e) {
	SelectionControl::resizeEvent(e);
}

QString Switch::current_value() {
	return this->value_ ? SWITCH_ACCEPT : SWITCH_REJECT;
};

void Switch::toggle(Qt::CheckState state) {
	if (state == Qt::Checked) {
		this->value_ = true;
		const QVariant pos_end = (int)(indicatorRect().width() - THUMB_WIDTH / 2 - 2);
		const QVariant thumb_end = colorFromOpacity(this->style_.thumb_on_brush_, this->style_.thumb_on_opacity_);
		const QVariant track_end = colorFromOpacity(this->style_.track_on_brush_, this->style_.track_on_opacity_);

		if (!isVisible()) {
			thumb_position_animation_->setCurrentValue(pos_end);
			thumb_brush_animation_->setCurrentValue(thumb_end);
			track_brush_animation_->setCurrentValue(track_end);
		}
		else {
			thumb_position_animation_->interpolate(0, pos_end);
			thumb_brush_animation_->interpolate(colorFromOpacity(this->style_.thumb_off_brush_, this->style_.thumb_off_opacity_), thumb_end);
			track_brush_animation_->interpolate(colorFromOpacity(this->style_.track_off_brush_, this->style_.track_off_opacity_), track_end);
		}
		if (this->mode_ == Switch::LinkMode::LinkToButton) {
			this->link_button_->setStyleSheet("QPushButton{color:#4994c4;}");
		}
		else if (this->mode_ == Switch::LinkMode::LinkToLabel) {
			this->link_label_->setStyleSheet("QLabel{color:#4994c4;}");
		}
	}
	else { // Qt::Unchecked
		this->value_ = false;
		const QVariant pos_end = (int)(THUMB_WIDTH / 2);
		const QVariant thumb_end = colorFromOpacity(this->style_.thumb_off_brush_, this->style_.thumb_off_opacity_);
		const QVariant track_end = colorFromOpacity(this->style_.track_off_brush_, this->style_.track_off_opacity_);

		if (!isVisible()) {
			thumb_position_animation_->setCurrentValue(pos_end);
			thumb_brush_animation_->setCurrentValue(thumb_end);
			track_brush_animation_->setCurrentValue(track_end);
		}
		else {
			thumb_position_animation_->interpolate(thumb_position_animation_->currentValue().toInt(), pos_end);
			thumb_brush_animation_->interpolate(colorFromOpacity(this->style_.thumb_on_brush_, this->style_.thumb_on_opacity_), thumb_end);
			track_brush_animation_->interpolate(colorFromOpacity(this->style_.track_on_brush_, this->style_.track_on_opacity_), track_end);
		}
		if (this->mode_ == Switch::LinkMode::LinkToButton) {
			this->link_button_->setStyleSheet("QPushButton{color:#666;}");
		}
		else if (this->mode_ == Switch::LinkMode::LinkToLabel) {
			this->link_label_->setStyleSheet("QLabel{color:#666;}");
		}
	}
}
