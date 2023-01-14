#include "main_window.h"
#include "settingctrls.hpp"
#include "treepanel.h"
#include "dialogs.h"
#include "../crosshatching/drawing.hpp"
#include "../crosshatching/util.hpp"
#include "../crosshatching/ink_layers.hpp"
#include "../crosshatching/brush_lang.hpp"
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
	constexpr auto k_controls_width = 400;
	constexpr auto k_swatch_scale = 4.0;

	QMenu* create_view_menu(ui::main_window* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* view_menu = new QMenu(tr("&View"));

		QAction* action_bw = new QAction(tr("Black and white"), parent);
		action_bw->setCheckable(true);
		parent->connect(action_bw, &QAction::toggled, parent, &ui::main_window::handle_view_bw_change);

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

	QMenu* create_file_menu(ui::main_window* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* file_menu = new QMenu(tr("&File"));

		QAction* action_open = new QAction(tr("&Open ..."), parent);
		parent->connect(action_open, &QAction::triggered, parent, &ui::main_window::open);
		file_menu->addAction(action_open);

		QAction* action_save = new QAction(tr("Save processed image ..."), parent);
		parent->connect(action_save, &QAction::triggered, parent, &ui::main_window::save_processed_image);
		file_menu->addAction(action_save);

		QAction* action_debug = new QAction(tr("Debug"), parent);
		parent->connect(action_debug, &QAction::triggered, parent, &ui::main_window::debug);
		file_menu->addAction(action_debug);

		file_menu->addSeparator();

		QAction* exit = new QAction(tr("&Exit"), parent);
		parent->connect(exit, &QAction::triggered, []() {QApplication::quit(); });
		file_menu->addAction(exit);

		return file_menu;
	}

	QMenu* create_crosshatching_menu(ui::main_window* parent) {
		auto tr = [parent](const char* str) {return parent->tr(str); };
		QMenu* ch_menu = new QMenu(tr("&Crosshatching"));

		QAction* actionOpen = new QAction(tr("&Settings"), parent);
		parent->connect(actionOpen, &QAction::triggered, parent, &ui::main_window::edit_settings);
		ch_menu->addAction(actionOpen);

		ch_menu->addSeparator();

		QAction* action_test = new QAction(tr("&Test"), parent);
		parent->connect(action_test, 
			&QAction::triggered, parent, &ui::main_window::test);
		ch_menu->addAction(action_test);

		QAction* redo_test_swatch = new QAction(tr("Change test swatch"), parent);
		parent->connect(redo_test_swatch, 
			&QAction::triggered, parent, &ui::main_window::redo_test_swatch);
		ch_menu->addAction(redo_test_swatch);

		ch_menu->addSeparator();

		QAction* action_generate = new QAction(tr("&Generate"), parent);
		parent->connect(action_generate, 
			&QAction::triggered, parent, &ui::main_window::generate);
		ch_menu->addAction(action_generate);

		return ch_menu;
	}

}

ui::main_window::main_window(QWidget *parent)
    : QMainWindow(parent)
{
	resize(QGuiApplication::primaryScreen()->availableGeometry().size() * 0.7);
	create_main_menu();

	QTabWidget* tab_ctrl = new QTabWidget();
	tab_ctrl->setTabPosition(QTabWidget::TabPosition::East);
	tab_ctrl->addTab( create_image_processing_pipeline_tools(), "image");
	tab_ctrl->addTab( create_drawing_tools(), "crosshatch");

    setCentralWidget(tab_ctrl);
	handle_source_image_change( img_proc_ctrls_.src );
}

void ui::main_window::open()
{
	auto image_filename = QFileDialog::getOpenFileName(this,
		tr("Open image"), "/home", tr("PNG (*.png);; Bitmap (*.bmp);; JPeg (*.jpeg *.jpg)"));
	auto str = image_filename.toStdString();
	if (!str.empty()) {
		img_proc_ctrls_.src = cv::imread(str.c_str());
		int n = img_proc_ctrls_.src.channels();
		display(img_proc_ctrls_.src);
		change_source_image(img_proc_ctrls_.src);
		img_proc_ctrls_.src_filename = fs::path(str).filename().string();
	}
}

void ui::main_window::test() { /*\
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

	int swatch_box_sz = test_swatch_picker::swatch_sz() * k_swatch_scale;
	if (swatch_box_sz != test_swatch.rows) {
		double scale = static_cast<double>(swatch_box_sz) / test_swatch.rows;
		test_swatch = ch::convert_to_3channel_grayscale(test_swatch);
		cv::resize(test_swatch, test_swatch, cv::Size(), scale, scale, cv::INTER_AREA);
	}

	set_swatch_view(test_swatch, false); */
}

void ui::main_window::set_swatch_view(cv::Mat swatch, bool left) {
	if (left) {
		crosshatching_.swatch_ = swatch;
		crosshatching_.img_swatch_->set_image(swatch, k_swatch_scale);
		crosshatching_.drawing_swatch_->set_image(swatch, k_swatch_scale);
	} else {
		crosshatching_.drawing_swatch_->set_image(swatch);
	}
	crosshatching_.set_view(drawing_tools::view::swatch);
}

void ui::main_window::set_layer_view() { 
    crosshatching_.layers_.clear();
	auto layers = layer_images();
	int n = static_cast<int>(layers.size());
	auto names = rv::concat(
			rv::iota(0, n) |
			rv::transform(
				[](int i)->std::string {
					return "layer " + std::to_string(i + 1);
				}
			), 
			rv::single(std::string("image"))
		);
	auto images = layers; // rv::concat(layers, rv::single(processed_image()));
	auto tab_content = rv::zip(names, images) |
		r::to<std::vector<std::tuple<std::string, cv::Mat>>>();
	crosshatching_.layer_viewer_->set_content(tab_content);
	crosshatching_.set_view(drawing_tools::view::layers);
}

void ui::main_window::set_drawing_view(cv::Mat drawing) {
	crosshatching_.drawing_->setFixedSize(drawing.cols, drawing.rows);
	crosshatching_.drawing_->set_image(drawing);
	crosshatching_.set_view(drawing_tools::view::drawing);
}

void ui::main_window::redo_test_swatch() {
	//set_swatch_view(ui::test_swatch_picker::get_test_swatch(img_proc_ctrls_.current), true);
}

void ui::main_window::generate() {
	if (!crosshatching_.layers_panel_->has_content()) {
		return;
	}
	auto progress_box = std::make_unique<ui::drawing_progress>( this, drawing_job() );
	progress_box->exec();
}

void ui::main_window::save_processed_image() {
	
}

void ui::main_window::edit_settings() {
	auto result = ui::settings::edit_settings(crosshatching_.params_);
	if (result) {
		crosshatching_.params_ = result.value();
	}
}

void ui::main_window::debug() {
    auto perlin = ch::perlin_noise({ 512,512 }, 12345, 8, 8.0);
    auto gray_perlin = ch::float_noise_to_grayscale(perlin);
    display(gray_perlin);

    ch::debug_brushes();
}

void ui::main_window::create_main_menu()
{
	menuBar()->addMenu( create_file_menu(this) );
	menuBar()->addMenu( create_view_menu(this) );
	menuBar()->addMenu( create_crosshatching_menu(this));
	connect(this, &main_window::change_source_image, this, &main_window::handle_source_image_change);
}

void ui::main_window::handle_source_image_change(cv::Mat& img) {
	for (auto* pipeline_stage : img_proc_ctrls_.pipeline) {
		pipeline_stage->initialize();
	}
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	auto crosshatch_menu = menuBar()->actions()[k_view_menu_index+1];
	auto enable_image_dependent_menus = !img.empty();
	view_menu->setEnabled(enable_image_dependent_menus);
	crosshatch_menu->setEnabled(enable_image_dependent_menus);
}

QWidget* ui::main_window::create_drawing_tools() {
	QSplitter* vert_splitter = new QSplitter();
	vert_splitter->setOrientation(Qt::Orientation::Vertical);
	crosshatching_.layers_panel_ = new layer_panel();
	vert_splitter->addWidget(
		crosshatching_.brushes_ = new brush_panel(
			*static_cast<layer_panel*>( crosshatching_.layers_panel_)
		)
	);
	vert_splitter->addWidget(crosshatching_.layers_panel_);

	vert_splitter->setMaximumWidth(100);
	QSplitter* splitter = new QSplitter();
	splitter->addWidget(vert_splitter);
	
	QWidget* picture_panel = new QWidget();
	auto layout = new QHBoxLayout(picture_panel);
	layout->addWidget( crosshatching_.img_swatch_ = new ui::image_box() );
	layout->addWidget( crosshatching_.drawing_swatch_ = new ui::image_box());

	int swatch_box_sz = test_swatch_picker::swatch_sz() * k_swatch_scale;
	crosshatching_.img_swatch_->setFixedSize(QSize(swatch_box_sz, swatch_box_sz));
	crosshatching_.drawing_swatch_->setFixedSize(QSize(swatch_box_sz, swatch_box_sz));

	QScrollArea* swatch_viewer = new QScrollArea();
	swatch_viewer->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	swatch_viewer->setWidget(picture_panel);

	crosshatching_.layer_viewer_ = new image_tab_ctrl();

	QScrollArea* drawing_scroller = new QScrollArea();
	drawing_scroller->setWidget(crosshatching_.drawing_ = new image_box());
	drawing_scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

	crosshatching_.viewer_stack_ = new QStackedWidget();

	crosshatching_.viewer_stack_->addWidget( swatch_viewer );
	crosshatching_.viewer_stack_->addWidget( crosshatching_.layer_viewer_ );
	crosshatching_.viewer_stack_->addWidget( drawing_scroller );

	splitter->addWidget(crosshatching_.viewer_stack_);
	splitter->setSizes(QList<int>({ 180,500 }));

	connect(crosshatching_.layers_panel_, &layer_panel::layers_changed,
		this, &main_window::set_layer_view);

	return splitter;
}

QWidget* ui::main_window::create_image_processing_pipeline_tools()
{
	QTabWidget* tab_ctrl = new QTabWidget();
	tab_ctrl->setMaximumWidth(300);

	img_proc_ctrls_.pipeline = std::vector<ui::image_processing_pipeline_item*>{
		new ui::scale_and_contrast(),
		new ui::edge_preserving_filter(),
		new ui::shock_filter(),
		new ui::mean_shift_segmentation(),
		new ui::raster_to_vector()
	};

	for (auto* stage: img_proc_ctrls_.pipeline) {
		stage->populate();
		tab_ctrl->addTab(stage, ("<" + std::to_string(stage->index() + 1) + ">").c_str());
		connect(stage, &ui::image_processing_pipeline_item::changed, 
			this, &main_window::handle_pipeline_change);
	}

	QScrollArea* scroller = new QScrollArea();
	scroller->setWidget(img_proc_ctrls_.img_box = new image_box());
	scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

	QWidget* ctrl = new QWidget();
	QHBoxLayout* horz_layout = new QHBoxLayout(ctrl);
	horz_layout->addWidget(tab_ctrl);
	horz_layout->addWidget(scroller);

	return ctrl;
}

bool is_pipeline_output_empty(ui::pipeline_output po) {
    if (std::holds_alternative<std::monostate>(po)) {
        return true;
    }
    if (std::holds_alternative<cv::Mat>(po) && std::get<cv::Mat>(po).empty()) {
        return true;
    }
    return false;
}

ui::pipeline_output ui::main_window::input_to_nth_stage(int index) const {
	if (index == 0) {
		return img_proc_ctrls_.src;
	} 

	pipeline_output input;
	int j = index;
	do {
		input = img_proc_ctrls_.pipeline.at(--j)->output();
	} while (is_pipeline_output_empty(input) && j > 0);

	if (is_pipeline_output_empty(input)) {
		input = img_proc_ctrls_.src;
	}

	return input;
}

cv::Mat ui::main_window::segmentation() const {
	return static_cast<ui::mean_shift_segmentation*>(img_proc_ctrls_.pipeline.back())->labels();
}

ch::crosshatching_job ui::main_window::drawing_job() const {
	return {
		image_src_filename(),
		layers(),
		drawing_params()
	};
}

std::string ui::main_window::image_src_filename() const {
	return img_proc_ctrls_.src_filename;
}

ui::vector_graphics_ptr ui::main_window::vector_output() const {
	if (std::holds_alternative<ui::vector_graphics_ptr>(img_proc_ctrls_.current)) {
		return std::get<ui::vector_graphics_ptr>(img_proc_ctrls_.current) ;
	} 
	return nullptr;
}

std::tuple<std::vector<ch::brush_expr_ptr>, std::vector<double>> 
			ui::main_window::brush_per_intervals() const {
	auto brush_dict = static_cast<brush_panel*>(crosshatching_.brushes_)->brush_dictionary();
	auto layer_name_vals = static_cast<layer_panel*>(crosshatching_.layers_panel_)->layers();
	return {
		layer_name_vals | rv::transform(
			[&](const auto& p)->ch::brush_expr_ptr {
				return brush_dict.at(std::get<0>(p));
			}
		) | r::to_vector,
		layer_name_vals | rv::transform(
			[&](const auto& p)->double {
				return std::get<1>(p);
			}
		) | r::to_vector,
	}; 
}

ch::ink_layers ui::main_window::layers() const {
	if (!vector_output()) {
		return {};
	}
    if (!crosshatching_.layers_.empty()) {
        return crosshatching_.layers_;
    }
	auto vo = *vector_output();
	auto [brushes, intervals] = brush_per_intervals();
	auto gray_value_polygons = ch::to_monochrome(vo.polygons, true);

	auto gray_value_ranges = intervals |
			rv::transform(
				[](double val) {
					return static_cast<uchar>(val * 255.0);
				}
			) | r::to_vector;

	return crosshatching_.layers_ = ch::split_into_layers(
        vo.sz, gray_value_polygons, brushes, gray_value_ranges
    );
}


ch::dimensions<int> ui::main_window::dimensions() const {
	auto curr_pipeline_output = img_proc_ctrls_.current;
	if (std::holds_alternative<cv::Mat>(curr_pipeline_output)) {
		return ch::mat_dimensions(std::get<cv::Mat>(curr_pipeline_output));
	} else if (std::holds_alternative<vector_graphics_ptr>(curr_pipeline_output)) {
		return std::get<vector_graphics_ptr>(curr_pipeline_output)->sz;
	} else {
		throw std::runtime_error("this shouldnt happen");
	}
}

std::vector<cv::Mat> ui::main_window::layer_images() const {
	
	double view_scale = img_proc_ctrls_.view_state.scale;
	auto sz = view_scale * dimensions();
	auto layers = this->layers().content;
	r::reverse(layers);

	auto all_polygons = layers |
		rv::transform(
			[view_scale](const auto& il) {
				return il.content |
					rv::transform(
						[view_scale](const ch::ink_layer_item& ili)->ch::gray_polygon {
							return { 
								ili.value, 
								ch::scale(ili.poly, view_scale) 
							};
						}
					);
			}
		) |
		rv::join |
		r::to_vector;
	auto composite = ch::paint_polygons(all_polygons, sz, true );

	return rv::concat(
		layers | 
			rv::transform(
				[view_scale, sz](const ch::ink_layer& il)->cv::Mat {
					return ch::paint_polygons(
						ch::scale(ch::to_polygons(il), view_scale),
						sz, 
						true
					);
				}
			), 
		rv::single(composite)
	) | r::to_vector;
}

ch::parameters ui::main_window::drawing_params() const {
	return crosshatching_.params_;
}

void ui::main_window::handle_pipeline_change(int index) {
	std::for_each(img_proc_ctrls_.pipeline.begin() + index, img_proc_ctrls_.pipeline.end(),
		[](ui::image_processing_pipeline_item* stage) {
			stage->clear_output();
		}
	);

	auto work_in_prog = input_to_nth_stage(index);
	if (is_pipeline_output_empty(work_in_prog)) {
		return;
	}
	std::for_each(img_proc_ctrls_.pipeline.begin() + index, img_proc_ctrls_.pipeline.end(),
		[&work_in_prog](ui::image_processing_pipeline_item* stage) {
			if (stage->is_on()) {
				work_in_prog = stage->process_image(work_in_prog);
			}
		}
	);

	display(work_in_prog);
}

void ui::main_window::display(pipeline_output work_in_prog) {
	if (! is_pipeline_output_empty(work_in_prog)) {
		img_proc_ctrls_.current = work_in_prog;
	}
	const auto& current = img_proc_ctrls_.current;
	double view_scale = img_proc_ctrls_.view_state.scale;
	cv::Mat mat;
	if (std::holds_alternative<cv::Mat>(current)) {
		mat = std::get<cv::Mat>(current);
		if (img_proc_ctrls_.view_state.black_and_white) {
			mat = ch::convert_to_3channel_grayscale(mat);
		}
		if (view_scale != 1.0) {
			mat = ch::scale(mat, view_scale);
		}
	} else if (std::holds_alternative<vector_graphics_ptr>(current)) {
		auto ptr_to_vg = std::get<vector_graphics_ptr>(current);
		std::vector<ch::colored_polygon> polys = (view_scale != 1.0) ?
			ch::scale(ptr_to_vg->polygons, view_scale) :
			ptr_to_vg->polygons;
		mat = ch::paint_polygons(polys, view_scale * ptr_to_vg->sz);
	}
	int wd = mat.cols;
	int hgt = mat.rows;
	img_proc_ctrls_.img_box->setFixedHeight(hgt);
	img_proc_ctrls_.img_box->setFixedWidth(wd);
	img_proc_ctrls_.img_box->set_image(mat);
}

std::tuple<int, int> ui::main_window::source_image_sz() const {
	return { img_proc_ctrls_.src.cols, img_proc_ctrls_.src.rows };
}

void ui::main_window::handle_view_scale_change(double scale) {
	img_proc_ctrls_.view_state.scale = scale;
	display();
	if (crosshatching_.layer_viewer_->has_images()) {
		set_layer_view();
	}
}

void ui::main_window::handle_view_bw_change(bool bw) {
	img_proc_ctrls_.view_state.black_and_white = bw;
	display();
}

ui::view_state ui::main_window::view_state() const {
	auto view_menu = menuBar()->actions()[k_view_menu_index];
	auto items = view_menu->children();
	return {};
}

void ui::drawing_tools::set_view(view v)
{
	viewer_stack_->setCurrentIndex( static_cast<int>(v) );
}

/*------------------------------------------------------------------------------------------------------*/

ui::layer_panel::layer_panel() :
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

void ui::layer_panel::set_brush_names(const std::vector<std::string>& brush_names) {
	brush_names_ = brush_names;
	setEnabled(!brush_names_.empty());
}

std::vector<std::tuple<std::string, double>> ui::layer_panel::layers() const {
	return layers_ |
		rv::transform(
			[](const auto& p)->std::tuple<std::string, double> {
				const auto& [val, name] = p;
				return { name, val };
			}
	) | r::to_vector;
}

bool ui::layer_panel::has_content() const {
	return !layers_.empty();
}

std::tuple<std::string, double> ui::layer_panel::row(int n) const {
	auto brush = list()->item(n, 0)->text().toStdString();
	auto val = list()->item(n, 2)->text().toDouble();
	return { brush,val };
}

void ui::layer_panel::add_layer() {
	auto result = ui::layer_dialog::create_layer_item(brush_names_, list()->rowCount() == 0);
    if (!result) {
        return;
    }
	auto [brush, end_of_range] = *result;
	insert_layer(brush, end_of_range);
	emit layers_changed();
}

void ui::layer_panel::cell_double_clicked(int r, int col) {
	auto [brush, val] = row(r);
	auto result = ui::layer_dialog::edit_layer_item(brush_names_, brush, val);
	if (result) {
		layers_.erase(val);
		auto [brush, end_of_range] = *result;
		insert_layer(brush, end_of_range);
		emit layers_changed();
	}
}

void ui::layer_panel::delete_layer() {
	auto current_row = list()->currentRow();
	if (current_row >= 0) {
		auto [brush, val] = row(current_row);
		layers_.erase(val);
		sync_layers_to_ui();
		emit layers_changed();
	}
}

void ui::layer_panel::insert_layer(const std::string& brush, double end_of_range) {
	layers_[end_of_range] = brush;
	sync_layers_to_ui();
}

void ui::layer_panel::setRowText(int row, const std::string& brush, 
		const std::string& from, const std::string& to) {
	std::array<QTableWidgetItem*, 3> items = {
		new QTableWidgetItem(brush.c_str()),
		new QTableWidgetItem(from.c_str()),
		new QTableWidgetItem(to.c_str())
	};
	for (auto [col, item] : rv::enumerate(items)) {
		item->setFlags(item->flags() & ~Qt::ItemIsEditable);
		list()->setItem(row, col, item);
	}
}

void ui::layer_panel::sync_layers_to_ui() {
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

/*------------------------------------------------------------------------------------------*/

ui::brush_panel::brush_panel(layer_panel& layers) :
		tree_panel("brushes", [&]() {this->add_brush_node(); }, 
			[&]() {this->delete_brush_node(); }),
		layer_panel_(layers) {
	connect(tree(), &QTreeWidget::itemDoubleClicked, this, &brush_panel::handle_double_click);

    auto brush1 = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                    ( 
                    brush
                        (define lines 
                            (quote
                                (brush 
                                    (horz_strokes
                                        (norm_rnd (lerp 0 800) (lerp 50 100))
                                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                                    )
                                    (dis (ramp 0.20 false true))
                                    (jiggle (norm_rnd 0.0 0.5))
                                )
                            )
                        )
                        (rot 45)
                        lines
                        (rot -45)
                        lines
                    )
            )"
        )
    );
    auto brush2 = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                (brush
                    (horz_strokes
                        (norm_rnd (lerp 0 800) (lerp 50 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                    (dis (ramp 0.2 0 1))
                    (jiggle (norm_rnd 0 0.5))
                )
            )"
        )
        );
    insert_toplevel_item(tree(), "main", brush1);
    insert_toplevel_item(tree(), "horz", brush2); 
    sync_layer_panel();
}

std::vector<std::string> ui::brush_panel::brush_names() const {
	std::vector<std::string> brushes(tree()->topLevelItemCount());
	for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
		QTreeWidgetItem* item = tree()->topLevelItem(i);
		brushes[i] = item->text(0).toStdString();
	}
	return brushes;
}

std::unordered_map<std::string, ch::brush_expr_ptr> ui::brush_panel::brush_dictionary() const {
	std::unordered_map<std::string, ch::brush_expr_ptr> dictionary;
	for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
		QTreeWidgetItem* item = tree()->topLevelItem(i);
		std::string name = item->text(0).toStdString();
		brush_item* bi = static_cast<brush_item*>(item);
		dictionary[name] = bi->brush_expression;
	}
	return dictionary;
}

ui::brush_panel::brush_item::brush_item(const std::string& name, ch::brush_expr_ptr expr) :
	QTreeWidgetItem(static_cast<QTreeWidget*>(nullptr), QStringList(QString(name.c_str()))),
	brush_expression(expr), 
	is_toplevel(true)
{}

ui::brush_panel::brush_item::brush_item(ch::brush_expr_ptr expr) :
	QTreeWidgetItem(
		static_cast<QTreeWidget*>(nullptr),
		QStringList(QString(expr->short_string().c_str()))),
	brush_expression(expr), 
	is_toplevel(false)
{}

bool ui::brush_panel::brush_item::is_leaf() const {
	return brush_expression->children().empty();
}

void ui::brush_panel::insert_brush_item(brush_item* parent, brush_item* item) {
	parent->addChild(item);
	auto expr = item->brush_expression;
	for (auto child : expr->children()) {
		insert_brush_item(item, new brush_item(child));
	}
}

void ui::brush_panel::insert_toplevel_item(QTreeWidget* tree, 
		const std::string& name, ch::brush_expr_ptr expr) {
	auto toplevel_item = new brush_item(name, expr);
	tree->addTopLevelItem(toplevel_item);

	for (auto child : expr->children()) {
		insert_brush_item(toplevel_item, new brush_item(child));
	}
}

void ui::brush_panel::add_brush_node() {
	auto result = ui::brush_dialog::create_brush();
	if (result) {
		const auto& [name, brush] = *result;
		insert_toplevel_item(tree(), name, brush);
		sync_layer_panel();
	}
}

void ui::brush_panel::delete_brush_node() {
	if (!tree()->selectedItems().empty()) {
		auto item = tree()->selectedItems().first();
		if (!item->parent()) {
			tree()->removeItemWidget(item, 0);
			delete item;
		}
	}
}

void ui::brush_panel::sync_layer_panel() {
	layer_panel_.set_brush_names(brush_names());
}

QTreeWidgetItem* toplevel_parent(QTreeWidgetItem* twi) {
	while (twi->parent()) {
		twi = twi->parent();
	}
	return twi;
}

void ui::brush_panel::handle_double_click(QTreeWidgetItem* item, int column) {
	/*
	brush_item* bi = static_cast<brush_item*>(item);
	if (bi->is_toplevel) {
		auto result = brush_dialog::edit_brush(bi->text(0).toStdString(), 
			bi->brush_expression->to_formatted_string());
		if (result) {
			auto [new_name, expr] = result.value();
			tree()->removeItemWidget(item, 0);
			delete item;
			insert_toplevel_item(tree(), new_name, expr);
		} 
	} else {
		if (bi->is_leaf()) {
			return; //TODO
		}
		auto result = brush_dialog::edit_brush_expr(bi->brush_expression->to_formatted_string());
		if (result) {
			auto toplevel_item = toplevel_parent(item);
			auto toplevel_brush_item = static_cast<brush_item*>(toplevel_item);
			auto parent = static_cast<brush_item*>(item->parent());
			parent->as_expr()->replace_child(bi->brush_expression, result);

			std::string test1 = parent->brush_expression->to_string();
			std::string toplevel = toplevel_brush_item->brush_expression->to_formatted_string();

			auto toplevel_parent_body = toplevel_brush_item->brush_expression;
			auto toplevel_parent_name = toplevel_brush_item->text(0).toStdString();
			tree()->removeItemWidget(toplevel_item, 0);
			delete toplevel_item;
			insert_toplevel_item(tree(), toplevel_parent_name, toplevel_parent_body);
		}
	}
	*/
}