#include "crosshatching.h"
#include "settingctrls.hpp"
#include "../crosshatching/drawing.hpp"
#include "../crosshatching/util.hpp"
#include "../crosshatching/brush_language.hpp"
#include <QtWidgets>
#include <QSlider>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <numbers>
#include <stdexcept>
#include <array>
#include <chrono>

namespace {

	constexpr auto k_view_menu_index = 1;
	constexpr auto k_controls_width = 300;

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

	QTabWidget* tab_ctrl = new QTabWidget();
	tab_ctrl->setTabPosition(QTabWidget::TabPosition::East);
	tab_ctrl->addTab( createImageProcPipelineCtrls(), "image");
	tab_ctrl->addTab( createCrosshatchCtrls(), "crosshatch");

    setCentralWidget(tab_ctrl);
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
	auto start_time = std::chrono::system_clock::now();
	std::string script = R"(
        (pipe
            (merge
                (pipe 
                    (lin_brush
                        (norm_rnd (lerp 0 800) (lerp 50 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                	(rot 45)
                    (dis (ramp 0.20 false true))
                )
                (pipe 
                    (ramp 0.50 true true)
                    (lin_brush
                        (norm_rnd (lerp 0 800) (lerp 50 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                    (rot 315)
                    (dis (ramp 0.20 false true))
                )
            )
            (jiggle (norm_rnd 0.0 0.02))
        )
	)";
    auto result = ch::parse_brush_language(script);
	if (std::holds_alternative<std::string>(result)) {
		QMessageBox msg_box;
		msg_box.setText(std::get<std::string>(result).c_str());
		msg_box.exec();
		return;
	}
	ch::brush br(std::get<ch::brush_fn>(result));
	br.build_n(10);

	auto param = br.gray_value_to_param(1.0);
	qDebug() << param;
	cv::Mat mat = current_image_;
	auto drawing = ch::generate_crosshatched_drawing(mat, 5.0, br);
	ch::to_svg("C:\\test\\testo.svg", drawing);

	std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - start_time;
	QMessageBox msg_box;
	msg_box.setText(std::to_string(elapsed.count()).c_str());
	msg_box.exec();
}

/*
void ui::crosshatching::generate() {
	auto make_pipeline_fn = [](double gray, double theta, double start) {
		return ch::make_run_pipeline_fn(
			ch::brush_pipeline{
				ch::make_ramp_fn(start, true,true),
				ch::make_linear_hatching_brush_fn(
					ch::make_lerped_normal_dist_fn( gray * 790, 50, 800, 100),
					ch::make_lerped_normal_dist_fn( (1-gray)*200, 20, 0, 0.05),
					ch::make_lerped_normal_dist_fn(7, 0.5, 0.5, 0.05),
					ch::make_default_hatching_unit()
				),
				ch::make_random_brush_adaptor(ch::jiggle, ch::normal_rnd_fn(0.0, 0.02)),
				ch::make_one_param_brush_adaptor(ch::rotate, ch::make_constant_fn(theta)),
				ch::make_ramp_fn(0.20, false, true),
				ch::disintegrate
			}
		);
	};

	cv::Mat mat = current_image_;
	//mat = ch::do_segmentation(mat, 8, 3, 12);
	std::vector<ch::hierarchical_brush_component> brushes = {
		[make_pipeline_fn](double gray) { return make_pipeline_fn(gray, 0, 0); },
		[make_pipeline_fn](double gray) { return make_pipeline_fn(gray, std::numbers::pi / 4.0, 0); },
		[make_pipeline_fn](double gray) { return make_pipeline_fn(gray, -std::numbers::pi / 4.0, 0); }
	};
	auto drawing = ch::generate_hierarchical_drawing(mat, 5.0, brushes, { 0.3333, 0.6666 });
	ch::to_svg("C:\\test\\output.svg", drawing);
}
*/

void ui::crosshatching::createMainMenu()
{
	menuBar()->addMenu( create_file_menu(this) );
	menuBar()->addMenu( create_view_menu(this) );
	connect(this, &crosshatching::change_source_image, this, &crosshatching::handle_source_image_change);
}

void ui::crosshatching::handle_source_image_change(cv::Mat& img) {
	for (auto* pipeline_stage : imgproc_pipeline_) {
		pipeline_stage->initialize();
	}
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	if (img.empty()) {
		view_menu->setEnabled(false);
		return;
	} 
	view_menu->setEnabled(true);
}

QWidget* ui::crosshatching::createCrosshatchCtrls() {
	return new QWidget();
}

QWidget* ui::crosshatching::createImageProcPipelineCtrls()
{
	QTabWidget* tab_ctrl = new QTabWidget();
	tab_ctrl->setMaximumWidth(k_controls_width);

	imgproc_pipeline_ = std::vector<ui::image_processing_pipeline_item*>{
		new ui::scale_and_contrast(),
		new ui::anisotropic_diffusion_filter(),
		new ui::shock_filter(),
		new ui::mean_shift_segmentation()
	};

	for (auto* stage: imgproc_pipeline_) {
		stage->populate();
		tab_ctrl->addTab(stage, (std::string("stage ") + std::to_string(stage->index() + 1)).c_str());
		connect(stage, &ui::image_processing_pipeline_item::changed, this, &crosshatching::handle_pipeline_change);
	}

	QScrollArea* scroller = new QScrollArea();
	scroller->setWidget(image_box_ = new QLabel());
	scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

	QSplitter* splitter = new QSplitter();
	splitter->addWidget(tab_ctrl);
	splitter->addWidget(scroller);
	splitter->setSizes(QList<int>({ INT_MAX, INT_MAX }));

	return splitter;
}

cv::Mat ui::crosshatching::input_to_nth_stage(int index) const {
	if (index == 0) {
		return src_image_;
	} 
	cv::Mat input;
	int j = index - 1;
	do {
		input = imgproc_pipeline_.at(j--)->output();
	} while (input.empty() && j > 1);
	if (input.empty()) {
		input = src_image_;
	}
	return input;
}

void ui::crosshatching::handle_pipeline_change(int index) {
	std::for_each(imgproc_pipeline_.begin() + index, imgproc_pipeline_.end(),
		[](ui::image_processing_pipeline_item* stage) {
			stage->clear_output();
		}
	);

	cv::Mat mat = input_to_nth_stage(index);
	std::for_each(imgproc_pipeline_.begin() + index, imgproc_pipeline_.end(),
		[&mat](ui::image_processing_pipeline_item* stage) {
			if (stage->is_on()) {
				mat = stage->process_image(mat);
			}
		}
	);

	display(mat);
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