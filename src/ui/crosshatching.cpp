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
#include <unordered_map>
#include <filesystem>

namespace r = ranges;
namespace rv = ranges::views;
namespace fs = std::filesystem;

namespace {

	constexpr auto k_view_menu_index = 1;
	constexpr auto k_controls_width = 300;
	constexpr auto k_swatch_scale = 4.0;

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

		QAction* action_open = new QAction(tr("&Open ..."), parent);
		parent->connect(action_open, &QAction::triggered, parent, &ui::crosshatching::open);
		file_menu->addAction(action_open);

		QAction* action_save = new QAction(tr("Save processed image ..."), parent);
		parent->connect(action_save, &QAction::triggered, parent, &ui::crosshatching::save_processed_image);
		file_menu->addAction(action_save);

		file_menu->addSeparator();

		QAction* exit = new QAction(tr("&Exit"), parent);
		parent->connect(exit, &QAction::triggered, []() {QApplication::quit(); });
		file_menu->addAction(exit);

		return file_menu;
	}

	QMenu* create_crosshatching_menu(ui::crosshatching* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* ch_menu = new QMenu(tr("&Crosshatching"));

		QAction* actionOpen = new QAction(tr("&Settings"), parent);
		parent->connect(actionOpen, &QAction::triggered, parent, &ui::crosshatching::edit_settings);
		ch_menu->addAction(actionOpen);

		ch_menu->addSeparator();

		QAction* action_test = new QAction(tr("&Test"), parent);
		parent->connect(action_test, &QAction::triggered, parent, &ui::crosshatching::test);
		ch_menu->addAction(action_test);

		QAction* redo_test_swatch = new QAction(tr("Change test swatch"), parent);
		parent->connect(redo_test_swatch, &QAction::triggered, parent, &ui::crosshatching::redo_test_swatch);
		ch_menu->addAction(redo_test_swatch);

		ch_menu->addSeparator();

		QAction* action_generate = new QAction(tr("&Generate"), parent);
		parent->connect(action_generate, &QAction::triggered, parent, &ui::crosshatching::generate);
		ch_menu->addAction(action_generate);

		return ch_menu;
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
			list()->setHorizontalHeaderLabels(QStringList{ "brush", "from value", "to value" });

			setEnabled(false);
		}

		void set_brush_names(const std::vector<std::string>& brush_names) {
			brush_names_ = brush_names;
			setEnabled(!brush_names_.empty());
		}

		std::vector<std::tuple<std::string, double>> layers() const {
			return layers_ |
				rv::transform(
					[](const auto& p)->std::tuple<std::string, double> {
						const auto& [val, name] = p;
						return { name, val };
					}
			) | r::to_vector;
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

	    void delete_layer() {
			auto current_row = list()->currentRow();
			if (current_row >= 0) {
				auto [brush, val] = row(current_row);
				layers_.erase(val);
				sync_layers_to_ui();
			}
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

		std::unordered_map<std::string, ch::brush_fn> brush_dictionary() const {
			std::unordered_map<std::string, ch::brush_fn> dictionary;
			for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
				QTreeWidgetItem* item = tree()->topLevelItem(i);
				std::string name = item->text(0).toStdString();
				brush_item* bi = static_cast<brush_item*>(item);
				dictionary[name] = std::get<ch::brush_fn>( bi->brush_expr->eval() );
			}
			return dictionary;
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
	handle_source_image_change( img_proc_ctrls_.src );
}

void ui::crosshatching::open()
{
	auto image_filename = QFileDialog::getOpenFileName(this,
		tr("Open image"), "/home", tr("PNG (*.png);; Bitmap (*.bmp);; JPeg (*.jpeg *.jpg)"));
	auto str = image_filename.toStdString();

	img_proc_ctrls_.src = cv::imread(str.c_str());
	display( img_proc_ctrls_.src );
	change_source_image( img_proc_ctrls_.src );
	img_proc_ctrls_.src_filename = fs::path(str).filename().string();
}

void ui::crosshatching::test() { 
	if (crosshatching_.swatch.empty()) {
		redo_test_swatch();
	} 
	auto box = std::make_unique<ui::progress>();
	auto layers = ch::generate_ink_layer_images(crosshatching_.swatch, {}, brush_per_intervals());
	auto result = box->run(
		"drawing test swatch...",
		[&]()->std::any {
			auto drawing = ch::generate_crosshatched_drawing(
				{ {}, layers, drawing_params() },
				{ box->progress_func(), {}, {} }
			);
			return { drawing };
		}
	);
	//TODO: error handling
	auto test = std::any_cast<ch::drawing>(result);
	test = ch::scale(test, k_swatch_scale);
	cv::Mat test_swatch = paint_drawing(test);
	crosshatching_.drawing_swatch->set_image(test_swatch); 
}

void ui::crosshatching::redo_test_swatch() {

	crosshatching_.swatch = ui::test_swatch_picker::get_test_swatch(img_proc_ctrls_.current);
	crosshatching_.img_swatch->set_image(crosshatching_.swatch, k_swatch_scale);
	crosshatching_.drawing_swatch->set_image(crosshatching_.swatch, k_swatch_scale);
}

/*
void ui::crosshatching::test() {
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
*/

void ui::crosshatching::generate() {
	auto progress_box = std::make_unique<ui::drawing_progress>( drawing_job() );
	progress_box->exec();
}

void ui::crosshatching::save_processed_image() {

}

void ui::crosshatching::edit_settings() {

}

void ui::crosshatching::createMainMenu()
{
	menuBar()->addMenu( create_file_menu(this) );
	menuBar()->addMenu( create_view_menu(this) );
	menuBar()->addMenu( create_crosshatching_menu(this));
	connect(this, &crosshatching::change_source_image, this, &crosshatching::handle_source_image_change);
}

void ui::crosshatching::handle_source_image_change(cv::Mat& img) {
	for (auto* pipeline_stage : img_proc_ctrls_.pipeline) {
		pipeline_stage->initialize();
	}
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	auto crosshatch_menu = menuBar()->actions()[k_view_menu_index+1];
	auto enable_image_dependent_menus = !img.empty();
	view_menu->setEnabled(enable_image_dependent_menus);
	crosshatch_menu->setEnabled(enable_image_dependent_menus);
}

QWidget* ui::crosshatching::createCrosshatchCtrls() {
	QSplitter* vert_splitter = new QSplitter();
	vert_splitter->setOrientation(Qt::Orientation::Vertical);
	crosshatching_.layers = new layer_panel();
	vert_splitter->addWidget(crosshatching_.brushes = new brush_panel(*static_cast<layer_panel*>( crosshatching_.layers )));
	vert_splitter->addWidget(crosshatching_.layers);

	vert_splitter->setMaximumWidth(100);
	QSplitter* splitter = new QSplitter();
	splitter->addWidget(vert_splitter);
	
	QWidget* picture_panel = new QWidget();
	auto layout = new QHBoxLayout(picture_panel);
	layout->addWidget( crosshatching_.img_swatch = new ui::image_box() );
	layout->addWidget( crosshatching_.drawing_swatch = new ui::image_box());

	int swatch_box_sz = test_swatch_picker::swatch_sz() * k_swatch_scale;
	crosshatching_.img_swatch->setFixedSize(QSize(swatch_box_sz, swatch_box_sz));
	crosshatching_.drawing_swatch->setFixedSize(QSize(swatch_box_sz, swatch_box_sz));

	QScrollArea* scroller = new QScrollArea();
	scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	scroller->setWidget(picture_panel);

	splitter->addWidget(scroller);
	splitter->setSizes(QList<int>({ 180,500 }));

	return splitter;
}

QWidget* ui::crosshatching::createImageProcPipelineCtrls()
{
	QTabWidget* tab_ctrl = new QTabWidget();
	tab_ctrl->setMaximumWidth(k_controls_width);

	img_proc_ctrls_.pipeline = std::vector<ui::image_processing_pipeline_item*>{
		new ui::scale_and_contrast(),
		new ui::anisotropic_diffusion_filter(),
		new ui::shock_filter(),
		new ui::mean_shift_segmentation()
	};

	for (auto* stage: img_proc_ctrls_.pipeline) {
		stage->populate();
		tab_ctrl->addTab(stage, (std::string("stage ") + std::to_string(stage->index() + 1)).c_str());
		connect(stage, &ui::image_processing_pipeline_item::changed, this, &crosshatching::handle_pipeline_change);
	}

	QScrollArea* scroller = new QScrollArea();
	scroller->setWidget(img_proc_ctrls_.img_box = new image_box());
	scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

	QSplitter* splitter = new QSplitter();
	splitter->addWidget(tab_ctrl);
	splitter->addWidget(scroller);
	splitter->setSizes(QList<int>({ INT_MAX, INT_MAX }));

	return splitter;
}

cv::Mat ui::crosshatching::input_to_nth_stage(int index) const {
	if (index == 0) {
		return img_proc_ctrls_.src;
	} 

	cv::Mat input;
	int j = index;
	do {
		input = img_proc_ctrls_.pipeline.at(--j)->output();
	} while (input.empty() && j > 0);

	if (input.empty()) {
		input = img_proc_ctrls_.src;
	}

	return input;
}

cv::Mat ui::crosshatching::segmentation() const {
	return static_cast<ui::mean_shift_segmentation*>(img_proc_ctrls_.pipeline.back())->labels();
}

ch::crosshatching_job ui::crosshatching::drawing_job() const {
	return {
		image_src_filename(),
		layers(),
		drawing_params()
	};
}

std::string ui::crosshatching::image_src_filename() const {
	return img_proc_ctrls_.src_filename;
}

cv::Mat ui::crosshatching::processed_image() const {
	return img_proc_ctrls_.current;
}

std::vector<std::tuple<ch::brush_fn, double>> ui::crosshatching::brush_per_intervals() const {
	auto brush_dict = static_cast<brush_panel*>(crosshatching_.brushes)->brush_dictionary();
	auto layer_name_vals = static_cast<layer_panel*>(crosshatching_.layers)->layers();
	return layer_name_vals |
			rv::transform(
				[&](const auto& p)->std::tuple<ch::brush_fn, double> {
					const auto& [brush_name, val] = p;
					return { brush_dict.at(brush_name), val };
				}
		) | r::to_vector;
}

std::vector<std::tuple<ch::brush_fn, cv::Mat>> ui::crosshatching::layers() const {
	auto brush_and_threshold = brush_per_intervals();
	return generate_ink_layer_images( processed_image(), segmentation(), brush_and_threshold );
}

std::vector<cv::Mat> ui::crosshatching::layer_images() const {
	auto layer_pairs = layers();
	return layer_pairs |
		rv::transform(
			[](const auto& tup)->cv::Mat {
				return std::get<1>(tup);
			}
		) | r::to_vector;
}

ch::crosshatching_params ui::crosshatching::drawing_params() const {
	//TODO
	return {};
}

void ui::crosshatching::handle_pipeline_change(int index) {
	std::for_each(img_proc_ctrls_.pipeline.begin() + index, img_proc_ctrls_.pipeline.end(),
		[](ui::image_processing_pipeline_item* stage) {
			stage->clear_output();
		}
	);

	cv::Mat mat = input_to_nth_stage(index);
	std::for_each(img_proc_ctrls_.pipeline.begin() + index, img_proc_ctrls_.pipeline.end(),
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
		img_proc_ctrls_.current = img;
	}

	cv::Mat mat = img_proc_ctrls_.current;
	if (img_proc_ctrls_.view_state.black_and_white) {
		mat = ch::convert_to_3channel_grayscale(mat);
	}
	if (img_proc_ctrls_.view_state.scale != 1.0) {
		mat = ch::scale(mat, img_proc_ctrls_.view_state.scale);
	}
	int wd = mat.cols;
	int hgt = mat.rows;
	int stride = mat.step;
	img_proc_ctrls_.img_box->setFixedHeight(hgt);
	img_proc_ctrls_.img_box->setFixedWidth(wd);
	img_proc_ctrls_.img_box->set_image(mat);
}

std::tuple<int, int> ui::crosshatching::source_image_sz() const {
	return { img_proc_ctrls_.src.cols, img_proc_ctrls_.src.rows };
}

void ui::crosshatching::handle_view_scale_change(double scale) {
	img_proc_ctrls_.view_state.scale = scale;
	display();
}

void ui::crosshatching::handle_view_bw_change(bool bw) {
	img_proc_ctrls_.view_state.black_and_white = bw;
	display();
}

ui::view_state ui::crosshatching::view_state() const {
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	auto items = view_menu->children();
	return {};
}