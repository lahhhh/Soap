#pragma once

#include "Identifier.h"

#include <QtWidgets>
#include <QPainter>

#include <QtCore/qeasingcurve.h>

#define cyan500 QColor("#00bcd4")
#define gray50 QColor("#fafafa")
#define gray400 QColor("#bdbdbd")

namespace SwitchStyle {

	using Type = QEasingCurve::Type;

	struct Animation {
		Animation() = default;
		Animation(Type _easing, int _duration) :easing{ _easing }, duration{ _duration } {

		}

		Type easing;
		int duration;
	};

	struct Switch {
		Switch() :
			indicator_width_{ 60 },
			height_{ 30 },
			indicator_margin_{ QMargins(2, 8, 2, 8) },
			thumb_on_brush_{ cyan500 },
			thumb_on_opacity_{ 1 },
			track_on_brush_{ cyan500 },
			track_on_opacity_{ 0.5 },
			thumb_off_brush_{ gray50 },
			thumb_off_opacity_{ 1 },
			track_off_brush_{ QColor(0, 0, 0) },
			track_off_opacity_{ 0.38 },
			thumb_brush_animation_{ Animation(Type::Linear, 150) },
			track_brush_animation_{ Animation(Type::Linear, 150) },
			thumb_position_animation_{ Animation(Type::InOutQuad, 150) } {

		}

		int indicator_width_;
		int height_;
		QMargins indicator_margin_;
		QColor thumb_on_brush_;
		double thumb_on_opacity_;
		QColor track_on_brush_;
		double track_on_opacity_;
		QColor thumb_off_brush_;
		double thumb_off_opacity_;
		QColor track_off_brush_;
		double track_off_opacity_;
		Animation thumb_brush_animation_;
		Animation track_brush_animation_;
		Animation thumb_position_animation_;
	};

}

class Animator final : public QVariantAnimation {
	Q_OBJECT
		Q_PROPERTY(QObject* targetObject READ targetObject WRITE setTargetObject)

public:
	Animator(QObject* target, QObject* parent = nullptr);
	~Animator() override;

	QObject* targetObject() const;
	void setTargetObject(QObject* target);

	inline bool isRunning() const {
		return state() == Running;
	}

public slots:
	void setup(int duration, QEasingCurve easing = QEasingCurve::Linear);
	void interpolate(const QVariant& start, const QVariant& end);
	void setCurrentValue(const QVariant&);

protected:
	void updateCurrentValue(const QVariant& value) override final;
	void updateState(QAbstractAnimation::State newState, QAbstractAnimation::State oldState) override final;

private:
	QPointer<QObject> target;
};

class SelectionControl : public QAbstractButton {
	Q_OBJECT

public:
	explicit SelectionControl(QWidget* parent = nullptr);

	Qt::CheckState checkState() const;

Q_SIGNALS:
	void stateChanged(int);

protected:
	void enterEvent(QEnterEvent* e) override;

	void checkStateSet() override;
	void nextCheckState() override;
	virtual void toggle(Qt::CheckState state) = 0;
};

class Switch final : public SelectionControl {
	Q_OBJECT;
	static constexpr auto CORNER_RADIUS = 8.0;
	static constexpr auto THUMB_RADIUS = 8.0;
	static constexpr auto THUMB_WIDTH = 36.0;

	enum class LinkMode { NoLink, LinkToButton, LinkToLabel };

public:
	Switch(bool default_value, QWidget* parent = nullptr);
	Switch(bool default_value, QPushButton* link_button, QWidget* parent = nullptr);
	Switch(bool default_value, QLabel* link_label, QWidget* parent = nullptr);

	Switch::LinkMode mode_ = Switch::LinkMode::NoLink;

	QSize sizeHint() const override final;

	QString current_value();
	void set_status(bool checked);

	bool value_;

protected:
	void paintEvent(QPaintEvent*) override final;
	void resizeEvent(QResizeEvent*) override final;
	void toggle(Qt::CheckState) override final;

	void init();
	
	QRect indicatorRect();

	static inline QColor colorFromOpacity(const QColor& c, qreal opacity) {
		return QColor(c.red(), c.green(), c.blue(), qRound(opacity * 255.0));
	}

private:
	SwitchStyle::Switch style_;
	QPushButton* link_button_ = nullptr;
	QLabel* link_label_ = nullptr;

	QPointer<Animator> thumb_brush_animation_;
	QPointer<Animator> track_brush_animation_;
	QPointer<Animator> thumb_position_animation_;
};
