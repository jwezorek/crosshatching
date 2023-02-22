#include "main_window.h"
#include "settingctrls.hpp"
#include "treepanel.h"
#include "dialogs.h"
#include "../crosshatching/drawing.hpp"
#include "../crosshatching/util.hpp"
#include "../crosshatching/ink_layers.hpp"
#include "../crosshatching/brush_lang.hpp"
#include "../crosshatching/json.hpp"
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

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;
namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {

	constexpr auto k_view_menu_index = 1;
	constexpr auto k_controls_width = 400;
	constexpr auto k_swatch_scale = 4.0;

    template<typename T>
    concept byte_collection = requires(const T & t) {
        std::same_as<decltype(t.size()), typename T::size_type>;
        std::same_as<decltype(t.data()), const unsigned char*>;
    };

    void write_data(std::ofstream& ofs, const byte_collection auto& bytes) {
        auto sz = bytes.size();
        ofs.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
        ofs.write(reinterpret_cast<const char*>(bytes.data()), sz);
    }

    template<byte_collection collection>
    void read_data(std::ifstream& ifs, collection& bytes) {
        typename collection::size_type sz;
        ifs.read(reinterpret_cast<char*>(&sz), sizeof(sz));
        bytes.resize(sz);
        ifs.read(reinterpret_cast<char*>(bytes.data()), sz);
    }

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

		QAction* action_open = new QAction(tr("&Open image ..."), parent);
		parent->connect(action_open, &QAction::triggered, parent, &ui::main_window::open);
		file_menu->addAction(action_open);

		QAction* action_save = new QAction(tr("Save processed image ..."), parent);
		parent->connect(action_save, &QAction::triggered, parent, &ui::main_window::save_processed_image);
		file_menu->addAction(action_save);
        file_menu->addSeparator();

        QAction* action_open_ch = new QAction(tr("Open project ..."), parent);
        parent->connect(action_open_ch, &QAction::triggered, parent, &ui::main_window::open_ch_file);
        file_menu->addAction(action_open_ch);

        QAction* action_save_ch = new QAction(tr("Save project ..."), parent);
        parent->connect(action_save_ch, &QAction::triggered, parent, &ui::main_window::save_ch_file);
        file_menu->addAction(action_save_ch);

        file_menu->addSeparator();

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
    tab_ctrl->addTab( create_region_map_tools(), "region map");

    connect(tab_ctrl, &QTabWidget::currentChanged, this, &main_window::tab_changed);

    setCentralWidget(tab_ctrl);
	handle_source_image_change( img_proc_ctrls_.src );
}

void ui::main_window::tab_changed(int index) {
    if (index == 2) {
        rgn_map_.rgn_props->repopulate_ctrls();
    }
}

const ui::brush_panel& ui::main_window::brush_panel() const {
    return *crosshatching_.brushes_;
}

std::string ui::main_window::state_to_json() const {
    json pipeline_settings = json::array();
    for (auto js_obj : img_proc_ctrls_.pipeline | 
            rv::transform([](auto ptr) {return ptr->to_json(); })) {
        pipeline_settings.push_back(js_obj);
    }

    auto brushes = crosshatching_.brushes_->to_json();
    auto layers = crosshatching_.layers_panel_->to_json();

    auto brush_dict = crosshatching_.brushes_->brush_name_dictionary();
    auto content = crosshatching_.layers_.to_json(brush_dict);

    json ch_state = {
        {"pipeline" , pipeline_settings},
        {"brushes" , brushes},
        {"layers", layers},
        {"content", std::move(content)}
    };

    return ch_state.dump();
}

void  ui::main_window::json_to_state(const std::string& str) {
    json ch_state = json::parse(str);

    const auto& pipeline_json = ch_state["pipeline"];
    for (auto item : img_proc_ctrls_.pipeline) {
        item->from_json(pipeline_json[item->index()]);
        emit item->changed(item->index());
    }

    crosshatching_.layers_panel_->from_json(ch_state["layers"]);
    crosshatching_.brushes_->from_json( ch_state["brushes"] );

    auto brush_dict = crosshatching_.brushes_->brush_dictionary();
    crosshatching_.layers_.from_json(ch_state["content"], brush_dict);

    sync_layers_to_ui();
}

void ui::main_window::save_project(const std::string& fname) const {
    std::vector<uchar> src_img_as_png;
    cv::imencode(".png", img_proc_ctrls_.src, src_img_as_png);
    auto json_str = state_to_json();

    std::ofstream ofs(fname, std::ios::out | std::ios::binary);
    write_data(ofs, img_proc_ctrls_.src_filename);
    write_data(ofs, src_img_as_png);
    write_data(ofs, json_str);
    ofs.close();
}

void ui::main_window::open_project(const std::string& fname) {
    std::ifstream ifs(fname, std::ios::in | std::ios::binary);
    std::vector<uchar> src_img_as_png;
    std::string src_fname;
    std::string json;

    read_data(ifs, src_fname);
    read_data(ifs, src_img_as_png);
    read_data(ifs, json);

    img_proc_ctrls_.src = cv::imdecode(src_img_as_png, cv::IMREAD_COLOR);
    int n = img_proc_ctrls_.src.channels();
    display(img_proc_ctrls_.src);
    change_source_image(img_proc_ctrls_.src);
    img_proc_ctrls_.src_filename = src_fname;

    json_to_state(json);
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

void ui::main_window::sync_layers_to_ui() {
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

    if (!rgn_map_.has_rgn_maps()) {
        auto ink_layers = this->layers();
        rgn_map_.rgn_props->set_layers(
            img_proc_ctrls_.view_state.scale,
            ink_layers
        );
    } else {
        rgn_map_.rgn_props->set_scale(img_proc_ctrls_.view_state.scale);
    }
}

void ui::main_window::rebuild_layers() {
    crosshatching_.layers_.clear();
    sync_layers_to_ui();
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
	//TODO
}

void ui::main_window::open_ch_file() {
    auto documentsPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    auto ch_filename = QFileDialog::getOpenFileName(this,
        tr("Open crosshatching project"), documentsPath, tr("Crosshatching projects Files (*.chp)"));
    auto str = ch_filename.toStdString();
    if (!str.empty()) {
        open_project(str);
    }
}

void ui::main_window::save_ch_file() {
    auto documentsPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save File"), documentsPath, tr("Crosshatching projects Files (*.chp)")
    );
    if (!filename.isNull()) {
        save_project(filename.toStdString());
    }

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
		crosshatching_.brushes_ = new ui::brush_panel(
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
		this, &main_window::rebuild_layers);

	return splitter;
}

QWidget* ui::main_window::create_region_map_tools() {
    QSplitter* vert_splitter = new QSplitter();
    vert_splitter->setOrientation(Qt::Orientation::Horizontal);
  
    //QScrollArea* scroller = new QScrollArea();
    //scroller->setWidget(rgn_map_.rgn_map_stack = new QStackedWidget());
    //scroller->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

    rgn_map_.rgn_map_stack = new QStackedWidget();
    rgn_map_.rgn_props = new rgn_map_panel(this, rgn_map_.rgn_map_stack);
    vert_splitter->addWidget(rgn_map_.rgn_props);
    vert_splitter->addWidget(rgn_map_.rgn_map_stack);
    vert_splitter->setSizes({100, 800});

    return vert_splitter;
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

ch::crosshatching_job ui::main_window::drawing_job() const {
	return {
		image_src_filename(),
		*layers(),
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
	auto brush_dict = static_cast<ui::brush_panel*>(crosshatching_.brushes_)->brush_dictionary();
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

ch::ink_layers* ui::main_window::layers() {
    if (crosshatching_.layers_.empty()) {
        if (!vector_output()) {
            return nullptr;
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

        crosshatching_.layers_ = ch::split_into_layers_simple(
            vo.sz, gray_value_polygons, brushes, gray_value_ranges
        );
    }
    return &crosshatching_.layers_;
}

bool ui::main_window::has_layers() const {
    return !crosshatching_.layers_.empty();
}

std::vector<std::string> ui::main_window::brush_names() const {
    return crosshatching_.brushes_->brush_names();
}

QStackedWidget* ui::main_window::regions_stack() const {
    return rgn_map_.rgn_map_stack;
}

const ch::ink_layers* ui::main_window::layers() const {
    auto unconst = const_cast<main_window*>(this);
    return unconst->layers();
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
	auto layers = this->layers()->content;

	auto all_polygons = layers |
		rv::transform(
			[view_scale](const auto& il) {
				return il |
					rv::transform(
						[view_scale](const ch::ink_layer_item& ili)->ch::gray_polygon {
							return { 
								ili.value, 
								ch::scale(ili.poly, view_scale) 
							};
						}
					);
			}
		) | rv::join | r::to_vector;
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
        sync_layers_to_ui();
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
	emit layers_changed();
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

json ui::layer_panel::to_json() const {
    json js_ary = json::array();
    for (const auto& [val, brush_name] : layers_) {
        json pair = json::array();
        pair.push_back(val);
        pair.push_back(brush_name);
        js_ary.push_back(pair);
    }
    return js_ary;
}

void ui::layer_panel::from_json(const json& json) {
    layers_.clear();
    for (const auto& pair : json) {
        auto thresh = pair[0].get<double>();
        auto brush_name = pair[1].get<std::string>();
        layers_[thresh] = brush_name;
        sync_layers_to_ui();
    }
}

/*------------------------------------------------------------------------------------------*/

ui::brush_panel::brush_panel(layer_panel& layers) :
		tree_panel("brushes", [&]() {this->add_brush_node(); }, 
			[&]() {this->delete_brush_node(); }),
		layer_panel_(layers) {
	connect(tree(), &QTreeWidget::itemDoubleClicked, this, &brush_panel::handle_double_click);

    auto interlocking_diagonal = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                    ( 
                    brush
                        (define lines 
                            (quote
                                (brush 
                                    (horz_strokes
                                        (norm_rnd (lerp 100 1600) (lerp 10 100))
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
    auto horz = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                (brush
                    (horz_strokes
                        (norm_rnd (lerp 100 1600) (lerp 10 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                    (dis (ramp 0.2 0 1))
                    (jiggle (norm_rnd 0 0.5))
                )
            )"
        )
    );
    auto vert = std::get<ch::brush_expr_ptr>(
        ch::parse(
            R"(
                (brush
                    (rot 90)
                    (horz_strokes
                        (norm_rnd (lerp 100 1600) (lerp 10 100))
                        (norm_rnd (lerp 200 0) (lerp 20 0.05))
                        (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
                    )
                    (dis (ramp 0.2 0 1))
                    (jiggle (norm_rnd 0 0.5))
                )
            )"
        )
    );
    auto solid = std::get<ch::brush_expr_ptr>(ch::parse("(solid value)"));
        
    insert_toplevel_item(tree(), "interlocking-diagonal", interlocking_diagonal);
    insert_toplevel_item(tree(), "horizontal", horz);
    insert_toplevel_item(tree(), "vertical", vert);
    insert_toplevel_item(tree(), "solid", solid);

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


std::vector<ch::brush_expr_ptr> ui::brush_panel::brushes() const {
    std::vector<ch::brush_expr_ptr> brushes(tree()->topLevelItemCount());
    for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
        brush_item* bi = static_cast<brush_item*>(tree()->topLevelItem(i));
        brushes[i] = bi->brush_expression;
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

std::unordered_map<ch::brush_expr*, std::string> ui::brush_panel::brush_name_dictionary() const {
    std::unordered_map<ch::brush_expr*, std::string> dict;
    for (int i = 0; i < tree()->topLevelItemCount(); ++i) {
        QTreeWidgetItem* item = tree()->topLevelItem(i);
        std::string name = item->text(0).toStdString();
        brush_item* bi = static_cast<brush_item*>(item);
        dict[bi->brush_expression.get()] = name;
    }
    return dict;
}

json ui::brush_panel::to_json() const {
    auto brushes = brush_dictionary();
    json brushes_json = json::array();
    for (const auto& [name, expr] : brushes) {
        json pair = json::array();
        pair.push_back(name);
        pair.push_back(ch::pretty_print(*expr));
        brushes_json.push_back(pair);
    }
    return brushes_json;
}

void ui::brush_panel::from_json(const ch::json& json) {
    tree()->clear();
    for (const auto& pair : json) {
        auto brush_name = pair[0].get<std::string>();
        auto brush_expr = std::get<ch::brush_expr_ptr>(ch::parse(pair[1].get<std::string>()));
        insert_toplevel_item(tree(), brush_name, brush_expr);
    }
    sync_layer_panel();
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
    if (expr) {
        for (auto child : expr->children()) {
            insert_brush_item(toplevel_item, new brush_item(child));
        }
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
            sync_layer_panel();
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
	brush_item* bi = static_cast<brush_item*>(item);
	if (bi->is_toplevel) {
        auto result = brush_dialog::edit_brush(
            bi->text(0).toStdString(),
            ch::pretty_print(*bi->brush_expression)
        );
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
		auto result = brush_dialog::edit_brush_expr(ch::pretty_print(*bi->brush_expression));
		if (result) {
			auto toplevel_item = toplevel_parent(item);
			auto toplevel_brush_item = static_cast<brush_item*>(toplevel_item);
			auto parent = static_cast<brush_item*>(item->parent());
			parent->brush_expression->replace_child(bi->brush_expression, result);

			auto toplevel_parent_body = toplevel_brush_item->brush_expression;
			auto toplevel_parent_name = toplevel_brush_item->text(0).toStdString();
			tree()->removeItemWidget(toplevel_item, 0);
			delete toplevel_item;
			insert_toplevel_item(tree(), toplevel_parent_name, toplevel_parent_body);
		}
	}
}