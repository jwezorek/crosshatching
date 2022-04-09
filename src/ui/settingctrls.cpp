#include "settingctrls.hpp"
#include "../crosshatching/util.hpp"

/*---------------------------------------------------------------------------------------------------------------------------*/

ui::image_processing_pipeline_item::image_processing_pipeline_item(QString title, int index) :
		content_area_(new QWidget()),
		index_(index) {
	QVBoxLayout* layout = new QVBoxLayout();
	QLabel* lbl_title = new QLabel(title);
	layout->addWidget(lbl_title);
	layout->addWidget(content_area_);
	setLayout(layout);
}

int ui::image_processing_pipeline_item::index() const {
	return index_;
}

cv::Mat ui::image_processing_pipeline_item::output() const {
	return output_;
}

void ui::image_processing_pipeline_item::clear_output() {
	output_ = {};
}

void ui::image_processing_pipeline_item::populate() {

}

void ui::image_processing_pipeline_item::initialize() {

}

cv::Mat ui::image_processing_pipeline_item::process_image(cv::Mat input) {
	return input;
}

bool ui::image_processing_pipeline_item::is_on() const {
	return false;
}

/*---------------------------------------------------------------------------------------------------------------------------*/

ui::scale_and_contrast::scale_and_contrast() :
	ui::image_processing_pipeline_item("Scale & Contrast", 0)
{}

void ui::scale_and_contrast::populate() {
	QHBoxLayout* layout = new QHBoxLayout();

	QGroupBox* contrast_box = new QGroupBox("contrast");
	QHBoxLayout* cb_layout = new QHBoxLayout();
	cb_layout->addWidget(contrast_slider_ = new ui::float_value_slider("beta", 1, 20, 1, 1000, ui::exp_val_to_pos, ui::exp_pos_to_val));
	cb_layout->addWidget(thresh_slider_ = new ui::float_value_slider("sigma", 0, 1, 0.5));
	contrast_box->setLayout(cb_layout);

	layout->addWidget(scale_slider_ = new ui::float_value_slider("scale", 20, 100, 100), 1);
	layout->addWidget(contrast_box, 2);
	content_area_->setLayout(layout);

	connect(scale_slider_, &float_value_slider::slider_released, [this]() { this->changed(this->index_); });
	connect(contrast_slider_, &float_value_slider::slider_released, [this]() { this->changed(this->index_); });
	connect(thresh_slider_, &float_value_slider::slider_released, [this]() { this->changed(this->index_); });
}

void ui::scale_and_contrast::initialize() {
	scale_slider_->set(100.0);
	contrast_slider_->set(1.0);
	thresh_slider_->set(0.5);
}

cv::Mat ui::scale_and_contrast::process_image(cv::Mat input) {

	auto scale = scale_slider_->value();

	scale /= 100.0;
	auto src_wd = input.cols;
	auto src_hgt = input.rows;
	int scaled_wd = static_cast<int>(scale * src_wd);
	int scaled_hgt = static_cast<int>(scale * src_hgt);

	cv::Mat scaled = ch::scale(input, scale);
	output_ = ch::apply_contrast(scaled, contrast_slider_->value(), thresh_slider_->value());
	return output_;
}

bool ui::scale_and_contrast::is_on() const {
	return scale_slider_->value() != 100.0 ||
		contrast_slider_->value() != 1.0 ||
		thresh_slider_->value() != 0.5;
}

/*---------------------------------------------------------------------------------------------------------------------------*/

ui::anisotropic_diffusion_filter::anisotropic_diffusion_filter() :
	ui::image_processing_pipeline_item("Anisotropic Diffusion", 1) {
}

void ui::anisotropic_diffusion_filter::populate() {
	QHBoxLayout* layout = new QHBoxLayout();

	layout->addWidget(alpha_slider_ = new float_value_slider("alpha", 0, 0.5, 0.15));
	layout->addWidget(k_slider_ = new float_value_slider("k", 0, 0.5, 0.15));
	layout->addWidget(iters_slider_ = new int_value_slider("iterations", 0, 25, 0));
	content_area_->setLayout(layout);

	connect(alpha_slider_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
	connect(k_slider_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
	connect(iters_slider_, &int_value_slider::slider_released, [this]() { changed(this->index()); });
}

void ui::anisotropic_diffusion_filter::initialize() {
	alpha_slider_->set(0.15);
	k_slider_->set(0.15);
	iters_slider_->set(0);
}

cv::Mat ui::anisotropic_diffusion_filter::process_image(cv::Mat input) {
	return output_ = ch::anisotropic_diffusion(
		input,
		alpha_slider_->value(),
		k_slider_->value(),
		iters_slider_->value()
	);
}

bool ui::anisotropic_diffusion_filter::is_on() const {
	return iters_slider_->value() > 0;
}

/*---------------------------------------------------------------------------------------------------------------------------*/

ui::shock_filter::shock_filter() :
	ui::image_processing_pipeline_item("Coherence Enhancing Shock Filter", 2) {
}

void ui::shock_filter::populate() {
		QHBoxLayout* layout = new QHBoxLayout();

		layout->addWidget(sigma_slider_ = new ui::int_value_slider("sigma 1", 1, 15, 5));
		layout->addWidget(str_sigma_slider_ = new ui::int_value_slider("sigma 2", 1, 15, 5));
		layout->addWidget(blend_slider_ = new ui::float_value_slider("blend", 0.0, 1.0, 0.5));
		layout->addWidget(iter_slider_ = new ui::int_value_slider("n", 0, 25, 0));
		content_area_->setLayout(layout);

		connect(sigma_slider_, &int_value_slider::slider_released, [this]() { changed(this->index()); });
		connect(str_sigma_slider_, &int_value_slider::slider_released, [this]() { changed(this->index()); });
		connect(blend_slider_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
		connect(iter_slider_, &int_value_slider::slider_released, [this]() { changed(this->index()); });
}

void ui::shock_filter::initialize() {
	sigma_slider_->set(5);
	str_sigma_slider_->set(5);
	blend_slider_->set(0.5);
	iter_slider_->set(0);
}

cv::Mat ui::shock_filter::process_image(cv::Mat input) {
	return output_ = ch::coherence_filter(
		input,
		2*sigma_slider_->value()+1,
		2*str_sigma_slider_->value()+1,
		blend_slider_->value(),
		iter_slider_->value()
	);
}

bool ui::shock_filter::is_on() const {
	return iter_slider_->value() > 0;
}

/*---------------------------------------------------------------------------------------------------------------------------*/

ui::mean_shift_segmentation::mean_shift_segmentation() :
	ui::image_processing_pipeline_item("Mean Shift Segmentation", 3) {
}

void ui::mean_shift_segmentation::populate() {
	QHBoxLayout* layout = new QHBoxLayout();

	layout->addWidget(sigmaS_slider_ = new ui::int_value_slider("sigma S", 0, 15, 0));
	layout->addWidget(sigmaR_slider_ = new ui::float_value_slider("sigma R", 0.0, 15, 3.0));
	layout->addWidget(min_size_slider_ = new ui::int_value_slider("min. size", 0, 100, 12));
	content_area_->setLayout(layout);

	connect(sigmaS_slider_, &int_value_slider::slider_released, [this]() { changed(this->index()); });
	connect(sigmaR_slider_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
	connect(min_size_slider_, &int_value_slider::slider_released, [this]() { changed(this->index()); });
}

void ui::mean_shift_segmentation::initialize() {
	sigmaS_slider_->set(0);
	sigmaR_slider_->set(3.0);
	min_size_slider_->set(12);
}

cv::Mat ui::mean_shift_segmentation::process_image(cv::Mat input) {
	return output_ = ch::do_segmentation(
		input,
		sigmaS_slider_->value(),
		sigmaR_slider_->value(),
		min_size_slider_->value()
	);
}

bool ui::mean_shift_segmentation::is_on() const {
	return sigmaS_slider_->value() != 0 && sigmaR_slider_->value() != 0.0;
}
