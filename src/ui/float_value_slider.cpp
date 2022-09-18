#include "float_value_slider.h"
#include "../crosshatching/util.hpp"

double ui::float_value_slider::range() const {
	return static_cast<double>(
		slider_->maximum()
		);
}

int ui::float_value_slider::value_to_position(double val) const {
	return val_to_pos_(val, range(), min_, max_);
}

double ui::float_value_slider::position_to_value(int pos) const {
	return pos_to_val_(pos, range(), min_, max_);
}

double ui::float_value_slider::value() const {
	return position_to_value(slider_->sliderPosition());
}

ui::float_value_slider::float_value_slider(const QString& txt, double min, double max, 
			double init_val, int range, val_to_pos_fn val_to_pos, pos_to_val_fn pos_to_val) :
		min_(min), max_(max),
		val_to_pos_(val_to_pos), pos_to_val_(pos_to_val) {
	QVBoxLayout* layout = new QVBoxLayout();
	QLabel* lbl_title = new QLabel(txt);

	slider_ = new QSlider();
	slider_->setOrientation(Qt::Vertical);
	slider_->setMinimum(0);
	slider_->setMaximum(range);
	slider_->setTickPosition(QSlider::TicksBothSides);
	slider_->setTickInterval(10);

	slider_->setSliderPosition(value_to_position(init_val));

	lbl_val_ = new QLabel(ch::to_string(init_val, 2).c_str());

	layout->addWidget(lbl_title);
	layout->addWidget(slider_);
	layout->addWidget(lbl_val_);
	this->setLayout(layout);

	connect(slider_, &QSlider::valueChanged, this, &float_value_slider::handle_position_change);
	connect(slider_, &QSlider::sliderReleased, this, &float_value_slider::handle_released);
}

void ui::float_value_slider::set(double value) {
	slider_->setSliderPosition(value_to_position(value));
}

void ui::float_value_slider::handle_position_change(int new_pos) {
	double new_value = position_to_value(new_pos);
	lbl_val_->setText(ch::to_string(new_value, 2).c_str());
	value_changed(new_value);
}

void ui::float_value_slider::handle_released() {
	slider_released();
}

int ui::float_value_slider::linear_value_to_position(double val, double range, double min, double max) {
	return static_cast<int>(range * ((val - min) / (max - min)));
}

double ui::float_value_slider::linear_position_to_value(int pos, double range, double min, double max) {
	return min + ((max - min) * pos) / range;
}

int ui::exp_val_to_pos(double val, double range, double min, double max) {
	double base = std::pow(max, 1.0 / range);
	return std::log(val) / std::log(base);
}

double ui::exp_pos_to_val(int pos, double range, double min, double max) {
	double base = std::pow(max, 1.0 / range);
	return std::pow(base, static_cast<double>(pos));
}

ui::int_value_slider::int_value_slider(const QString& txt, int min, int max, int init_val)
{
	QVBoxLayout* layout = new QVBoxLayout();
	QLabel* lbl_title = new QLabel(txt);

	slider_ = new QSlider();
	slider_->setOrientation(Qt::Vertical);
	slider_->setMinimum(min);
	slider_->setMaximum(max);
	slider_->setTickPosition(QSlider::TicksBothSides);
	slider_->setTickInterval(1);

	slider_->setSliderPosition(init_val);

	lbl_val_ = new QLabel(std::to_string(init_val).c_str());

	layout->addWidget(lbl_title);
	layout->addWidget(slider_);
	layout->addWidget(lbl_val_);
	this->setLayout(layout);

	connect(slider_, &QSlider::valueChanged, this, &int_value_slider::handle_position_change);
	connect(slider_, &QSlider::sliderReleased, this, &int_value_slider::handle_released);
}

void ui::int_value_slider::set(int val) {
	slider_->setSliderPosition(val);
}

int ui::int_value_slider::value() const {
	return slider_->sliderPosition();
}

void ui::int_value_slider::handle_position_change(int pos) {
	lbl_val_->setText(std::to_string(pos).c_str());
	value_changed(pos);
}

void ui::int_value_slider::handle_released() {
	slider_released();
}
