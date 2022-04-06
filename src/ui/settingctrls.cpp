#include "settingctrls.hpp"
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

ui::float_value_slider::float_value_slider(const QString& txt, double min, double max, double init_val, int range, 
		val_to_pos_fn val_to_pos, pos_to_val_fn pos_to_val) :
	min_(min), max_(max),
	val_to_pos_(val_to_pos), pos_to_val_(pos_to_val)
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

int exp_val_to_pos(double val, double range, double min, double max) {
	double base = std::pow(max, 1.0 / range);
	return std::log(val) / std::log(base);
}

double exp_pos_to_val(int pos, double range, double min, double max) {
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

ui::preprocess_settings::preprocess_settings() :
	contrast_slider_(new ui::float_value_slider("beta", 1, 20, 1, 1000, exp_val_to_pos, exp_pos_to_val)),
	thresh_slider_(new ui::float_value_slider("sigma", 0, 1, 0.5)),
	scale_slider_(new ui::float_value_slider("scale", 20, 100, 100))
{
	QHBoxLayout* layout = new QHBoxLayout();

	QGroupBox* contrast_box = new QGroupBox("contrast");
	QHBoxLayout* cb_layout = new QHBoxLayout();
	cb_layout->addWidget(contrast_slider_);
	cb_layout->addWidget(thresh_slider_);
	contrast_box->setLayout(cb_layout);

	layout->addWidget(scale_slider_, 1);
	layout->addWidget(contrast_box, 2);
	setLayout(layout);

	connect(scale_slider_, &float_value_slider::slider_released, [this]() { scale_changed( this->scale() ); });
	connect(contrast_slider_, &float_value_slider::slider_released, [this]() { contrast_changed( this->contrast_params() ); });
	connect(thresh_slider_, &float_value_slider::slider_released, [this]() { contrast_changed( this->contrast_params() ); });
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

ui::shock_filter_settings::shock_filter_settings() :
	sigma_slider_(new ui::int_value_slider("sigma 1", 1, 15, 5)),
	str_sigma_slider_(new ui::int_value_slider("sigma 2", 1, 15, 5)),
	blend_slider_(new ui::float_value_slider("blend", 0.0, 2.5, 0.5)),
	iter_slider_(new ui::int_value_slider("n", 0, 25, 0))
{
	QHBoxLayout* layout = new QHBoxLayout();

	layout->addWidget(sigma_slider_);
	layout->addWidget(str_sigma_slider_);
	layout->addWidget(blend_slider_);
	layout->addWidget(iter_slider_);
	setLayout(layout);

	connect(sigma_slider_, &int_value_slider::slider_released, [this]() { changed(state()); });
	connect(str_sigma_slider_, &int_value_slider::slider_released, [this]() { changed(state()); });
	connect(blend_slider_, &float_value_slider::slider_released, [this]() { changed(state()); });
	connect(iter_slider_, &int_value_slider::slider_released, [this]() { changed(state()); });
}

void ui::shock_filter_settings::initialize()
{
	sigma_slider_->set(11);
	str_sigma_slider_->set(11);
	blend_slider_->set(0.5);
	iter_slider_->set(0);
}


std::tuple<int, int, double, int>  ui::shock_filter_settings::state() const {
	return {
		sigma_slider_->value(),
		str_sigma_slider_->value(),
		blend_slider_->value(),
		iter_slider_->value()
	};
}