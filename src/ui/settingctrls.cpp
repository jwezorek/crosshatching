#include "settingctrls.hpp"
#include "../crosshatching/util.hpp"
#include "../crosshatching/raster_to_vector.hpp"

/*------------------------------------------------------------------------------------------------------*/

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

ui::pipeline_output ui::image_processing_pipeline_item::output() const {
	return output_;
}

void ui::image_processing_pipeline_item::clear_output() {
	output_ = {};
}

void ui::image_processing_pipeline_item::populate() {

}

void ui::image_processing_pipeline_item::initialize() {

}

ui::pipeline_output ui::image_processing_pipeline_item::process_image(pipeline_output input) {
	return input;
}

bool ui::image_processing_pipeline_item::is_on() const {
	return false;
}

ch::json ui::image_processing_pipeline_item::to_json() const {
    return {};
}

void ui::image_processing_pipeline_item::from_json(const ch::json& js) {
    
}

/*------------------------------------------------------------------------------------------------------*/

ui::scale_and_contrast::scale_and_contrast() :
	ui::image_processing_pipeline_item("Scale & Contrast", 0)
{}


std::tuple<float, float> ui::scale_and_contrast::bw_cutoff() const {
    float upper = static_cast<float>(black_white_cutoff_->GetUpperValue());
    float lower = static_cast<float>(black_white_cutoff_->GetLowerValue());
    return {
        1.0 - upper / 255.0,
        1.0 - lower / 255.0
    };
}

void ui::scale_and_contrast::populate() {
	QHBoxLayout* layout = new QHBoxLayout();

	QGroupBox* contrast_box = new QGroupBox("contrast");
	QHBoxLayout* cb_layout = new QHBoxLayout();
	cb_layout->addWidget(contrast_slider_ = 
		new ui::float_value_slider("beta", 1, 20, 1, 1000, ui::exp_val_to_pos, ui::exp_pos_to_val)
	);
	cb_layout->addWidget(thresh_slider_ = new ui::float_value_slider("sigma", 0, 1, 0.5));
	contrast_box->setLayout(cb_layout);

	layout->addWidget(scale_slider_ = new ui::float_value_slider("scale", 20, 100, 100), 1);
	layout->addWidget(contrast_box, 2);

    QVBoxLayout* bw_layout = new QVBoxLayout();
    bw_layout->addWidget(new QLabel("b/w cutoff"));
    bw_layout->addWidget(black_white_cutoff_ = new RangeSlider(Qt::Orientation::Vertical));
    bw_layout->addWidget(black_cutoff_ = new QLabel("1.0"));
    bw_layout->addWidget(white_cutoff_ = new QLabel("0.0"));
    black_white_cutoff_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding);
    black_white_cutoff_->SetRange(0, 255);
    layout->addLayout(bw_layout);

	content_area_->setLayout(layout);

	connect(scale_slider_, &float_value_slider::slider_released, 
		[this]() { this->changed(this->index_); });
	connect(contrast_slider_, &float_value_slider::slider_released, 
		[this]() { this->changed(this->index_); });
	connect(thresh_slider_, &float_value_slider::slider_released, 
		[this]() { this->changed(this->index_); });


    connect(black_white_cutoff_, &RangeSlider::upperHandleReleased,
        [this]() {
            float val = static_cast<float>(black_white_cutoff_->GetUpperValue());
            val = 1.0 - val / 255.0;
            white_cutoff_->setText(
                std::to_string(val).c_str()
            );
            this->changed(this->index_); 
        }
    );

    connect(black_white_cutoff_, &RangeSlider::lowerHandleReleased,
        [this]() {
            float val = static_cast<float>(black_white_cutoff_->GetLowerValue());
            val = 1.0 - val / 255.0;
            black_cutoff_->setText(
                std::to_string(val).c_str()
            );
            this->changed(this->index_);
        }
    );
}

void ui::scale_and_contrast::initialize() {
	scale_slider_->set(100.0);
	contrast_slider_->set(1.0);
	thresh_slider_->set(0.5);
}

ui::pipeline_output ui::scale_and_contrast::process_image(pipeline_output inp) {

	auto input = std::get<cv::Mat>(inp);
	auto scale = scale_slider_->value();

	scale /= 100.0;
	auto src_wd = input.cols;
	auto src_hgt = input.rows;
	int scaled_wd = static_cast<int>(scale * src_wd);
	int scaled_hgt = static_cast<int>(scale * src_hgt);

    auto [white, black] = bw_cutoff();

	cv::Mat scaled = ch::scale(input, scale);
	output_ = ch::apply_contrast(
        scaled, 
        contrast_slider_->value(), 
        thresh_slider_->value(),
        white, 
        black
    );
	return output_;
}

bool ui::scale_and_contrast::is_on() const {
    auto [white, black] = bw_cutoff();
    return scale_slider_->value() != 100.0 ||
        contrast_slider_->value() != 1.0 ||
        thresh_slider_->value() != 0.5 ||
        white != 0.0f ||
        black != 1.0f;
}

ch::json ui::scale_and_contrast::to_json() const {
    return ch::json{
        {"scale", scale_slider_->value()},
        {"contrast", contrast_slider_->value()},
        {"threshold", thresh_slider_->value()},
        {"bw_upper", black_white_cutoff_->GetUpperValue()},
        {"bw_lower", black_white_cutoff_->GetLowerValue()},
    };
}

void ui::scale_and_contrast::from_json(const ch::json& js) {
    scale_slider_->set(js["scale"].get<double>());
    contrast_slider_->set(js["contrast"].get<double>());
    thresh_slider_->set(js["threshold"].get<double>());
    black_white_cutoff_->SetUpperValue(js["bw_upper"].get<double>());
    black_white_cutoff_->SetLowerValue(js["bw_lower"].get<double>());
}


/*------------------------------------------------------------------------------------------------------*/
ui::edge_preserving_filter::edge_preserving_filter() :
    ui::image_processing_pipeline_item("Smoothing", 1) {
}

void ui::edge_preserving_filter::populate() {
    QVBoxLayout* v_layout = new QVBoxLayout();
    QHBoxLayout* h_layout = new QHBoxLayout();

    h_layout->addWidget(sigma_s_ = new float_value_slider("sigma S", 0, 200, 60));
    h_layout->addWidget(sigma_r_ = new float_value_slider("sigma R", 0, 1.0, 0.4));
    v_layout->addLayout(h_layout);
    v_layout->addWidget(filter_type_ = new QComboBox());

    filter_type_->addItem("Off");
    filter_type_->addItem("Recursive Filtering");
    filter_type_->addItem("Normalized Convolution Filtering");

    content_area_->setLayout(v_layout);

    connect(sigma_s_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
    connect(sigma_r_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
    connect(filter_type_, &QComboBox::currentIndexChanged, [this]() { changed(this->index()); });
}

void ui::edge_preserving_filter::initialize() {
    sigma_s_->set(60);
    sigma_r_->set(0.4);
    filter_type_->setCurrentIndex(0);
}

ui::pipeline_output ui::edge_preserving_filter::process_image(pipeline_output inp) {
    auto input = std::get<cv::Mat>(inp);
    return output_ = ch::edge_preserving_smoothing(
        input,
        filter_type_->currentIndex(),
        sigma_s_->value(),
        sigma_r_->value()
    );
}

bool ui::edge_preserving_filter::is_on() const {
    return filter_type_->currentIndex() > 0;
}

ch::json ui::edge_preserving_filter::to_json() const {
    return ch::json{
        {"sigma_s", sigma_s_->value()},
        {"sigma_r", sigma_r_->value()},
        {"filter_type", filter_type_->currentIndex()}
    };
}

void ui::edge_preserving_filter::from_json(const ch::json& js) {
    filter_type_->disconnect();

    sigma_s_->set(js["sigma_s"].get<double>());
    sigma_r_->set(js["sigma_r"].get<double>());
    filter_type_->setCurrentIndex(js["filter_type"].get<int>());

    connect(filter_type_, &QComboBox::currentIndexChanged, [this]() { changed(this->index()); });
}

/*------------------------------------------------------------------------------------------------------*/

ui::stylization_filter::stylization_filter() :
    ui::image_processing_pipeline_item("Stylization", 1) {
}

void ui::stylization_filter::populate() {
    QHBoxLayout* layout = new QHBoxLayout();

    layout->addWidget(sigma_s_ = new float_value_slider("sigma S", 0, 200, 60));
    layout->addWidget(sigma_r_ = new float_value_slider("sigma R", 0, 1.0, 0.0));
    content_area_->setLayout(layout);

    connect(sigma_s_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
    connect(sigma_r_, &float_value_slider::slider_released, [this]() { changed(this->index()); });
}

void ui::stylization_filter::initialize() {
    sigma_s_->set(60);
    sigma_r_->set(0);
}

ui::pipeline_output ui::stylization_filter::process_image(pipeline_output inp) {
    auto input = std::get<cv::Mat>(inp);
    return output_ = ch::stylize(
        input,
        sigma_s_->value(),
        sigma_r_->value()
    );
}

bool ui::stylization_filter::is_on() const {
    return sigma_r_->value() > 0;
}

/*------------------------------------------------------------------------------------------------------*/

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

ui::pipeline_output ui::anisotropic_diffusion_filter::process_image(pipeline_output inp) {
	auto input = std::get<cv::Mat>(inp);
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

/*------------------------------------------------------------------------------------------------------*/

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

		connect(sigma_slider_, &int_value_slider::slider_released, 
			[this]() { changed(this->index()); });
		connect(str_sigma_slider_, &int_value_slider::slider_released, 
			[this]() { changed(this->index()); });
		connect(blend_slider_, &float_value_slider::slider_released, 
			[this]() { changed(this->index()); });
		connect(iter_slider_, &int_value_slider::slider_released, 
			[this]() { changed(this->index()); });
}

void ui::shock_filter::initialize() {
	sigma_slider_->set(5);
	str_sigma_slider_->set(5);
	blend_slider_->set(0.5);
	iter_slider_->set(0);
}

ui::pipeline_output ui::shock_filter::process_image(pipeline_output inp) {
	auto input = std::get<cv::Mat>(inp);
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

ch::json ui::shock_filter::to_json() const {
    return ch::json{
       {"sigma", sigma_slider_->value()},
       {"str_sigma", str_sigma_slider_->value()},
       {"blend", blend_slider_->value()},
       {"iters", iter_slider_->value()}
    };
}

void ui::shock_filter::from_json(const ch::json& js) {
    sigma_slider_->set(js["sigma"].get<int>());
    str_sigma_slider_->set(js["str_sigma"].get<int>());
    blend_slider_->set(js["blend"].get<double>());
    iter_slider_->set(js["iters"].get<int>());
}

/*------------------------------------------------------------------------------------------------------*/

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

ui::pipeline_output ui::mean_shift_segmentation::process_image(pipeline_output inp) {
	auto input = std::get<cv::Mat>(inp);
    cv::Mat label_image;
	std::tie(output_, label_image) = ch::meanshift_segmentation(
		input,
		sigmaS_slider_->value(),
		sigmaR_slider_->value(),
		min_size_slider_->value()
	);
	return output_;
}

bool ui::mean_shift_segmentation::is_on() const {
	return sigmaS_slider_->value() != 0 && sigmaR_slider_->value() != 0.0;
}

ch::json ui::mean_shift_segmentation::to_json() const {
    return {
        { "sigma_s", sigmaS_slider_->value() },
        { "sigma_r", sigmaR_slider_->value()},
        { "min_size", min_size_slider_->value()},
    };
}

void ui::mean_shift_segmentation::from_json(const ch::json& js) {
    sigmaS_slider_->set(js["sigma_s"].get<int>());
    sigmaR_slider_->set(js["sigma_r"].get<double>());
    min_size_slider_->set(js["min_size"].get<int>());
}

/*------------------------------------------------------------------------------------------------------*/

ui::raster_to_vector::raster_to_vector() :
	ui::image_processing_pipeline_item("Raster to vector", 4) {
}

void ui::raster_to_vector::populate() {
	QHBoxLayout* layout = new QHBoxLayout();

	layout->addWidget(param_slider_ = new ui::float_value_slider("RDP param", 0.0, 8.0, 0.0, 100));
	content_area_->setLayout(layout);

	connect(param_slider_, &float_value_slider::slider_released, [this]() {
		changed(this->index()); }
	);

}

void ui::raster_to_vector::initialize() {
	param_slider_->set(0.0);
}

ui::pipeline_output ui::raster_to_vector::process_image(pipeline_output inp) {
	auto input = std::get<cv::Mat>(inp);
	auto param = param_slider_->value();
	if (param == 0.0) {
		return input;
	}
	auto colored_polys = ch::raster_to_vector(input, param);

	return output_ = std::make_shared<vector_graphics>(
		std::move(colored_polys),
		ch::mat_dimensions(input)
	);
}

bool ui::raster_to_vector::is_on() const {
	return param_slider_->value() > 0.0;
}

ch::json ui::raster_to_vector::to_json() const {
    return ch::json{ {"param", param_slider_->value()} };
}

void ui::raster_to_vector::from_json(const ch::json& js) {
    param_slider_->set(js["param"].get<double>());
}
