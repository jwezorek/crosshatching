#include "crosshatching.h"
#include "../crosshatching/drawing.hpp"
#include <QtWidgets>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <numbers>
#include <stdexcept>
#include <array>

namespace {

	QMenu* create_view_menu(crosshatching* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* view_menu = new QMenu(tr("&View"));

		QAction* action_bw = new QAction(tr("Black and white"), parent);
		action_bw->setCheckable(true);

		QActionGroup* scale_actions = new QActionGroup(parent);
		std::array<int, 5> scales = { 100, 125, 150, 175, 200 };
		for (int scale : scales) {
			std::string str = std::to_string(scale) + "%";
			auto scale_action = new QAction(QString::fromStdString(str), parent);
			scale_action->setCheckable(true);
			scale_actions->addAction(scale_action);
		}
		scale_actions->setExclusive(true);
		scale_actions->setEnabled(true);

		view_menu->addAction(action_bw);
		view_menu->addSeparator()->setText(tr("Scale"));
		view_menu->addActions(scale_actions->actions());

		return view_menu;
	}

	QMenu* create_file_menu(crosshatching* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* file_menu = new QMenu(tr("&File"));

		QAction* actionOpen = new QAction(tr("&Open"), parent);
		parent->connect(actionOpen, &QAction::triggered, parent, &crosshatching::open);
		file_menu->addAction(actionOpen);

		QAction* action_generate = new QAction(tr("&Generate"), parent);
		parent->connect(action_generate, &QAction::triggered, parent, &crosshatching::generate);
		file_menu->addAction(action_generate);

		return file_menu;
	}
}

crosshatching::crosshatching(QWidget *parent)
    : QMainWindow(parent)
{
	resize(QGuiApplication::primaryScreen()->availableGeometry().size() * 0.7);
	createMainMenu();
    setCentralWidget(createContent());
}

void crosshatching::open()
{
	auto image_filename = QFileDialog::getOpenFileName(this,
		tr("Open image"), "", tr("PNG Files (*.png)"));
	auto str = image_filename.toStdString();
	src_image_ = cv::imread(str.c_str());
	

	int wd = src_image_.cols;
	int hgt = src_image_.rows;
	int stride = src_image_.step;

	image_box_->setFixedHeight(hgt);
	image_box_->setFixedWidth(wd);
	image_box_->setPixmap(QPixmap::fromImage(QImage(src_image_.data, wd,  hgt, stride, QImage::Format_BGR888)));
}

void crosshatching::generate() {
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

void crosshatching::createMainMenu()
{
	menuBar()->addMenu(create_file_menu(this));
	menuBar()->addMenu(create_view_menu(this));
}

QWidget* crosshatching::createContent()
{
	QTabWidget* tab = new QTabWidget();
	tab->setMaximumWidth(300);

	tab->addTab(new QWidget(), tr("Contrast"));
	tab->addTab(new QWidget(), tr("filter"));
	tab->addTab(new QWidget(), tr("segmentation"));

	QScrollArea* scroller = new QScrollArea();
	scroller->setWidget(image_box_ = new QLabel());
	scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

	QSplitter* splitter = new QSplitter();
	splitter->addWidget(tab);
	splitter->addWidget(scroller);
	splitter->setSizes(QList<int>({ INT_MAX, INT_MAX }));

	return splitter;
}
