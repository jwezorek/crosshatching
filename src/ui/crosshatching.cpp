#include "crosshatching.h"
#include "settingctrls.hpp"
#include "../crosshatching/drawing.hpp"
#include "../crosshatching/util.hpp"
#include <QtWidgets>
#include <QSlider>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <numbers>
#include <stdexcept>
#include <array>

namespace {

	constexpr auto k_view_menu_index = 1;
	constexpr auto k_controls_width = 400;

	QMenu* create_view_menu(ui::crosshatching* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* view_menu = new QMenu(tr("&View"));

		QAction* action_bw = new QAction(tr("Black and white"), parent);
		action_bw->setCheckable(true);
		parent->connect(action_bw, &QAction::toggled, parent, &ui::crosshatching::handle_view_bw_change);

		QActionGroup* scale_actions = new QActionGroup(parent);
		std::array<int,10> scales = { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500 };
		for (int scale : scales) {
			std::string str = std::to_string(scale) + "%";
			auto scale_action = new QAction(QString::fromStdString(str), parent);
			scale_action->setCheckable(true);
			scale_actions->addAction(scale_action);
			double new_scale_value = static_cast<double>(scale) / 100.0;
			parent->connect(scale_action, &QAction::triggered,
				[parent, new_scale_value]() {
					parent->handle_view_scale_change(new_scale_value);
				}
			);
		}
		scale_actions->setExclusive(true);
		scale_actions->setEnabled(true);

		view_menu->addAction(action_bw);
		view_menu->addSeparator()->setText(tr("Scale"));
		view_menu->addActions(scale_actions->actions());

		return view_menu;
	}

	QMenu* create_file_menu(ui::crosshatching* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* file_menu = new QMenu(tr("&File"));

		QAction* actionOpen = new QAction(tr("&Open"), parent);
		parent->connect(actionOpen, &QAction::triggered, parent, &ui::crosshatching::open);
		file_menu->addAction(actionOpen);

		QAction* action_generate = new QAction(tr("&Generate"), parent);
		parent->connect(action_generate, &QAction::triggered, parent, &ui::crosshatching::generate);
		file_menu->addAction(action_generate);

		return file_menu;
	}
}

ui::crosshatching::crosshatching(QWidget *parent)
    : QMainWindow(parent)
{
	resize(QGuiApplication::primaryScreen()->availableGeometry().size() * 0.7);
	createMainMenu();
    setCentralWidget(createContent());
	handle_source_image_change(src_image_);
}

void ui::crosshatching::open()
{
	auto image_filename = QFileDialog::getOpenFileName(this,
		tr("Open image"), "/home", tr("PNG (*.png);; Bitmap (*.bmp);; JPeg (*.jpeg *.jpg)"));
	auto str = image_filename.toStdString();

	src_image_ = cv::imread(str.c_str());
	display(src_image_);
	change_source_image(src_image_);
}

void ui::crosshatching::generate() {
	//src_image_ = ch::do_segmentation(src_image_, 8, 3.0f, 12);
	//image_box_->setPixmap(QPixmap::fromImage(QImage(src_image_.data, src_image_.cols, src_image_.rows, src_image_.step, QImage::Format_BGR888)));
	try {
		auto make_pipeline_fn = [](double theta, double start) {
			return ch::make_run_pipeline_fn(
				ch::brush_pipeline{
					ch::make_ramp_fn(start, true,true),
					ch::make_linear_hatching_brush_fn(
						ch::make_lerped_normal_dist_fn(0, 50, 800, 100),
						ch::make_lerped_normal_dist_fn(200, 20, 0, 0.05),
						ch::make_lerped_normal_dist_fn(7, 0.5, 0.5, 0.05),
						ch::make_default_hatching_unit()
					),
					ch::make_one_param_brush_adaptor(ch::rotate, ch::make_constant_fn(theta)),
					ch::make_ramp_fn(0.20, false, true),
					ch::disintegrate
				}
			);
		};

		ch::brush br(
			ch::make_run_pipeline_fn(
				ch::brush_pipeline{
					ch::make_merge_fn({
						make_pipeline_fn(std::numbers::pi / 4.0, 0),
						make_pipeline_fn(-std::numbers::pi / 4.0, 0.65)
					}),
					ch::make_random_brush_adaptor(ch::jiggle, ch::normal_rnd_fn(0.0, 0.02))
				}
			),
			1
		);
		br.build_n(10);
		auto drawing = ch::generate_crosshatched_drawing("C:\\test\\paris.png", { 8, 3.0f, 12 }, 4.0, br);
		ch::to_svg("C:\\test\\paris-2.svg", drawing); 

	} catch (std::runtime_error err) {
		QMessageBox msg_box;
		msg_box.setText(err.what());
		msg_box.exec();
	}
}

void ui::crosshatching::createMainMenu()
{
	menuBar()->addMenu( create_file_menu(this) );
	menuBar()->addMenu( create_view_menu(this) );
	connect(this, &crosshatching::change_source_image, this, &crosshatching::handle_source_image_change);
}

void ui::crosshatching::handle_source_image_change(cv::Mat& img) {
	preprocess_settings_->initialize();
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	if (img.empty()) {
		view_menu->setEnabled(false);
		return;
	} 
	view_menu->setEnabled(true);
}

QWidget* ui::crosshatching::createContent()
{
	QTabWidget* tab = new QTabWidget();
	tab->setMaximumWidth(k_controls_width);

	tab->addTab(preprocess_settings_ = new ui::preprocess_settings(), tr("pre"));
	tab->addTab(new QWidget(), tr("filter"));
	tab->addTab(new QWidget(), tr("post"));

	QScrollArea* scroller = new QScrollArea();
	scroller->setWidget(image_box_ = new QLabel());
	scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

	QSplitter* splitter = new QSplitter();
	splitter->addWidget(tab);
	splitter->addWidget(scroller);
	splitter->setSizes(QList<int>({ INT_MAX, INT_MAX }));

	connect(preprocess_settings_, &preprocess_settings::scale_changed, this, &crosshatching::handle_scale_change);
	connect(preprocess_settings_, &preprocess_settings::contrast_changed, this, &crosshatching::handle_contrast_changed);

	return splitter;
}

void ui::crosshatching::display(cv::Mat img) {
	if (!img.empty()) {
		current_image_ = img;
	}

	cv::Mat mat = current_image_;
	if (view_state_.black_and_white) {
		mat = ch::convert_to_3channel_grayscale(mat);
	}
	if (view_state_.scale != 1.0) {
		mat = ch::scale(mat, view_state_.scale);
	}
	int wd = mat.cols;
	int hgt = mat.rows;
	int stride = mat.step;
	image_box_->setFixedHeight(hgt);
	image_box_->setFixedWidth(wd);
	image_box_->setPixmap(QPixmap::fromImage(QImage(mat.data, wd, hgt, stride, QImage::Format_BGR888)));
}


void ui::crosshatching::preprocess_image(double scale, double beta, double sigma) {
	if (src_image_.empty()) {
		return;
	}

	if (scale < 0) {
		scale = preprocess_settings_->scale();
	}

	if (beta < 0) {
		beta = preprocess_settings_->beta();
	}

	if (sigma < 0) {
		sigma = preprocess_settings_->sigma();
	}

	scale /= 100.0;
	auto [src_wd, src_hgt] = source_image_sz();
	int scaled_wd = static_cast<int>(scale * src_wd);
	int scaled_hgt = static_cast<int>(scale * src_hgt);

	preprocess_settings_->scaled_image = ch::scale(src_image_, scale);
	display(preprocess_settings_->contrast_applied_image = ch::apply_contrast(preprocess_settings_->scaled_image, beta, sigma));
}

void ui::crosshatching::handle_scale_change(double new_scale) {
	preprocess_image(new_scale, -1, -1);
	auto vs = view_state();
}

void ui::crosshatching::handle_contrast_changed(std::tuple<double, double> params) {
	auto [beta, sigma] = params;
	preprocess_image(-1, beta, sigma);
}

std::tuple<int, int> ui::crosshatching::source_image_sz() const {
	return { src_image_.cols, src_image_.rows };
}

void ui::crosshatching::handle_view_scale_change(double scale) {
	view_state_.scale = scale;
	display();
}

void ui::crosshatching::handle_view_bw_change(bool bw) {
	view_state_.black_and_white = bw;
	display();
}

ui::view_state ui::crosshatching::view_state() const {
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	auto items = view_menu->children();
	return {};
}