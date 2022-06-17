#include "crosshatching.h"
#include "settingctrls.hpp"
#include "treepanel.h"
#include "dialogs.h"
#include "../crosshatching/drawing.hpp"
#include "../crosshatching/util.hpp"
#include "../crosshatching/brush_language.hpp"
#include <QtWidgets>
#include <QSlider>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <range/v3/all.hpp>
#include <numbers>
#include <stdexcept>
#include <array>
#include <chrono>
#include <fstream>
#include <memory>
#include <map>
#include <array>

namespace r = ranges;
namespace rv = ranges::views;

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

		QAction* action_debug = new QAction(tr("&Debug"), parent);
		parent->connect(action_debug, &QAction::triggered, parent, &ui::crosshatching::debug);
		file_menu->addAction(action_debug);

		return file_menu;
	}

	class layer_panel : public ui::list_panel {

	public:
		layer_panel() :
			ui::list_panel("layers", 3, [&]() { this->add_layer(); }, [&]() { this->delete_layer(); }) {

			//layers_->horizontalHeader()->setStretchLastSection(true);
			list()->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
			list()->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
			list()->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Fixed);
			list()->setColumnWidth(1, 75);
			list()->setColumnWidth(2, 75);
			//list()->setSelectionMode(QAbstractItemView::NoSelection);
			list()->setSelectionBehavior(QAbstractItemView::SelectRows);
			list()->connect(list(), &QTableWidget::cellDoubleClicked, this, &layer_panel::cell_double_clicked);
			setEnabled(false);
		}

		void set_brush_names(const std::vector<std::string>& brush_names) {
			brush_names_ = brush_names;
			setEnabled(!brush_names_.empty());
		}

	private:

		std::vector<std::string> brush_names_;
		std::map<double, std::string> layers_;


		std::tuple<std::string, double> row(int n) const {
			auto brush = list()->item(n, 0)->text().toStdString();
			auto val = list()->item(n, 2)->text().toDouble();
			return { brush,val };
		}

		void add_layer() {
			auto result = ui::layer_dialog::create_layer_item(brush_names_, list()->rowCount() == 0);
			if (result) {
				auto [brush, end_of_range] = *result;
				insert_layer(brush, end_of_range);
			}
		}

		void cell_double_clicked(int r, int col) {
			auto [brush, val] = row(r);
			auto result = ui::layer_dialog::edit_layer_item(brush_names_, brush, val);
			if (result) {
				layers_.erase(val);
				auto [brush, end_of_range] = *result;
				insert_layer(brush, end_of_range);
			}
		}

		static void delete_layer() {

		}

		void insert_layer(const std::string& brush, double end_of_range) {
			layers_[end_of_range] = brush;
			sync_layers_to_ui();
		}

		void setRowText(int row, const std::string& brush, const std::string& from, const std::string& to) {
			std::array<QTableWidgetItem*, 3> items = {
				new QTableWidgetItem(brush.c_str()),
				new QTableWidgetItem(from.c_str()),
				new QTableWidgetItem(to.c_str())
			};
			for (auto [col, item] : rv::enumerate(items)) {
				item->setFlags(item->flags() & ~Qt::ItemIsEditable);
				list()->setItem(row, col, item );
			}
		}

		void sync_layers_to_ui() {
			list()->setRowCount(layers_.size());
			int row = 0;
			std::string prev = "0.0";
			for (const auto& [val, name] : layers_) {
				std::string curr = std::to_string(val);
				setRowText(row, name, prev, curr);
				prev = curr;
				++row;
			}
		}
	};

	class brush_panel : public ui::tree_panel {

	public:

		brush_panel( layer_panel& layers ) : 
				tree_panel("brushes", [&]() {this->add_brush_node(); }, [&]() {this->delete_brush_node(); }),
				layer_panel_(layers) {

		}

		std::vector<std::string> brush_names() const {
			std::vector<std::string> brushes(tree()->topLevelItemCount());
			for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
				QTreeWidgetItem* item = tree()->topLevelItem(i);
				brushes[i] = item->text(0).toStdString();
			}
			return brushes;
		}

	private:

		layer_panel& layer_panel_;

		class brush_item : public QTreeWidgetItem {
		public:
			brush_item(const std::string& name, ch::brush_expr_ptr expr) :
				QTreeWidgetItem(static_cast<QTreeWidget*>(nullptr), QStringList(QString(name.c_str()))),
				brush_expr(expr)
			{}

			brush_item(ch::brush_expr_ptr expr) :
				QTreeWidgetItem(static_cast<QTreeWidget*>(nullptr), QStringList(QString(expr->to_short_string().c_str()))),
				brush_expr(expr)
			{}

			ch::brush_expr_ptr brush_expr;
		};

		static void insert_brush_item(brush_item* parent, brush_item* item) {
			parent->addChild(item);
			if (item->brush_expr->is_expression()) {
				auto expr = item->brush_expr;
				auto children = expr->children();
				if (!children) {
					return;
				}
				for (auto child : *children) {
					insert_brush_item(item, new brush_item(child));
				}
			}
		}

		static void insert_toplevel_item(QTreeWidget* tree, const std::string& name, ch::brush_expr_ptr expr) {
			auto toplevel_item = new brush_item(name, expr);
			tree->addTopLevelItem(toplevel_item);
			auto children = expr->children();

			if (!children) {
				return;
			}

			for (auto child : *children) {
				insert_brush_item(toplevel_item, new brush_item(child));
			}
		}

	   void add_brush_node() {
			if (tree()->selectedItems().empty()) {
				auto result = ui::brush_dialog::create_brush();
				if (result) {
					const auto& [name, brush] = *result;
					insert_toplevel_item(tree(), name, brush);
					sync_layer_panel();
				}
			} else {
				// TODO: add child item
			}
		}

		void delete_brush_node() {
			if (!tree()->selectedItems().empty()) {
				QMessageBox mb;
				mb.setText( tree()->selectedItems().first()->text(0) + " delete");
				mb.exec();
			}
		}

		void sync_layer_panel() {
			layer_panel_.set_brush_names(brush_names());
		}
	};

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
    auto result = ch::brush_language_to_func(script);
	if (std::holds_alternative<std::string>(result)) {
		QMessageBox msg_box;
		msg_box.setText(std::get<std::string>(result).c_str());
		msg_box.exec();
		return;
	}
	//ch::brush br(std::get<ch::brush_fn>(result));
	//br.build_n(10);

	cv::Mat mat = current_image_;
	auto drawing = ch::generate_crosshatched_drawing(mat, { {std::get<ch::brush_fn>(result), 1.0} }, ch::crosshatching_params(5.0));
	ch::to_svg("C:\\test\\testo.svg", drawing);

	std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - start_time;
	QMessageBox msg_box;
	msg_box.setText(std::to_string(elapsed.count()).c_str());
	msg_box.exec();
}

void ui::crosshatching::debug() {
	ch::debug(current_image_, segmentation());
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
	QSplitter* vert_splitter = new QSplitter();
	vert_splitter->setOrientation(Qt::Orientation::Vertical);
	layers_ = new layer_panel();
	vert_splitter->addWidget(brushes_ = new brush_panel(*static_cast<layer_panel*>(layers_)));
	vert_splitter->addWidget(layers_);

	vert_splitter->setMaximumWidth(100);
	QSplitter* splitter = new QSplitter();
	splitter->addWidget(vert_splitter);
	splitter->addWidget(new QWidget());
	splitter->setSizes(QList<int>({ 120,500 }));

	return splitter;
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

cv::Mat ui::crosshatching::segmentation() const {
	return static_cast<ui::mean_shift_segmentation*>(imgproc_pipeline_.back())->labels();
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