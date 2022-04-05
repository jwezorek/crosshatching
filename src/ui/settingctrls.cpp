#include "settingctrls.hpp"
#include "../crosshatching/util.hpp"

double ui::labeled_slider::range() const {
	return static_cast<double>(
		slider_->maximum()
	);
}

int ui::labeled_slider::value_to_position(double val) const {
	return static_cast<int>(range() * ((val - min_) / (max_ - min_)));
}

double ui::labeled_slider::position_to_value(int pos) const {
	return min_ + ((max_ - min_) * pos) / range();
}

double ui::labeled_slider::value() const {
	return position_to_value(slider_->sliderPosition());
}

ui::labeled_slider::labeled_slider(const QString& txt, double min, double max, double init_val, int range) :
	min_(min), max_(max)
{
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

	connect(slider_, &QSlider::valueChanged, this, &labeled_slider::handle_position_change);
	connect(slider_, &QSlider::sliderReleased, this, &labeled_slider::handle_released);
}

void ui::labeled_slider::set(double value) {
	slider_->setSliderPosition(value_to_position(value));
}

void ui::labeled_slider::handle_position_change(int new_pos) {
	double new_value = position_to_value(new_pos);
	lbl_val_->setText(ch::to_string(new_value, 2).c_str());
	value_changed(new_value);
}

void ui::labeled_slider::handle_released() {
	slider_released();
}

ui::preprocess_settings::preprocess_settings() :
	contrast_slider_(new ui::labeled_slider("beta", 1, 20, 1)),
	thresh_slider_(new ui::labeled_slider("sigma", 0, 1, 0.5)),
	scale_slider_(new ui::labeled_slider("Scale", 20, 100, 100))
{
	QHBoxLayout* layout = new QHBoxLayout();

	QGroupBox* contrast_box = new QGroupBox("Contrast");
	QHBoxLayout* cb_layout = new QHBoxLayout();
	cb_layout->addWidget(contrast_slider_);
	cb_layout->addWidget(thresh_slider_);
	contrast_box->setLayout(cb_layout);

	layout->addWidget(scale_slider_, 1);
	layout->addWidget(contrast_box, 2);
	setLayout(layout);

	connect(scale_slider_, &labeled_slider::slider_released, [this]() { scale_changed( this->scale() ); });
	connect(contrast_slider_, &labeled_slider::slider_released, [this]() { contrast_changed( this->contrast_params() ); });
	connect(thresh_slider_, &labeled_slider::slider_released, [this]() { contrast_changed( this->contrast_params() ); });
}

void ui::preprocess_settings::initialize() {
	contrast_slider_->set(1.0);
	thresh_slider_->set(0.5);
	scale_slider_->set(100.0);
}

double ui::preprocess_settings::scale() const {
	return scale_slider_->value();
}

double ui::preprocess_settings::beta() const {
	return contrast_slider_->value();
}

double ui::preprocess_settings::sigma() const {
	return thresh_slider_->value();
}

std::tuple<double, double> ui::preprocess_settings::contrast_params() const {
	return { contrast_slider_->value(), thresh_slider_->value() };
}