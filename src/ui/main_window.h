#pragma once

#include "../crosshatching/drawing.hpp"
#include "../crosshatching/brush_lang.hpp"
#include "../crosshatching/ink_layers.hpp"
#include "../crosshatching/util.hpp"
#include "layer_panel.hpp"
#include "brush_panel.hpp"
#include "rgn_tool_panel.hpp"
#include "settingctrls.hpp"
#include "treepanel.h"
#include "image_box.h"
#include "image_tab_ctrl.h"
#include <QtWidgets>
#include <QtWidgets/QMainWindow>
#include <opencv2/core.hpp>
#include <tuple>

/*------------------------------------------------------------------------------------------------*/

namespace ui {

    struct view_state {
        double scale;
        bool black_and_white;
        view_state() : scale(1.0), black_and_white(false)
        {}
    };

    struct image_processing_tools {
        std::string src_filename;
        std::vector<image_processing_pipeline_item*> pipeline;
        image_box* img_box;
        cv::Mat src;
        pipeline_output current;
        view_state view_state;
    };

    struct drawing_tools {

        enum class view : int {
            swatch = 0,
            layers = 1,
            drawing = 2
        };

        brush_panel* brushes_;
        layer_panel* layers_panel_;
        image_box* img_swatch_;
        image_box* drawing_swatch_;
        image_box* drawing_;
        image_tab_ctrl* layer_viewer_;
        QStackedWidget* viewer_stack_;
        cv::Mat swatch_;
        ch::parameters params_;
        mutable ch::ink_layers layers_;

        void set_view(view v);
    };

    using layer_tuple = std::tuple<ch::brush_expr_ptr, std::vector<ch::gray_polygon>>;

    class main_window : public QMainWindow
    {

        Q_OBJECT

    public:
        main_window(QWidget* parent = Q_NULLPTR);

        void open();
        void save_processed_image();
        void open_ch_file();
        void save_ch_file();
        void generate();
        void test();
        void redo_test_swatch();
        void edit_settings();
        void debug();

        void set_swatch_view(cv::Mat swatch, bool left);

        void sync_layers_to_ui();
        void rebuild_layers();
        void set_drawing_view(cv::Mat drawing);
        std::vector<std::string> brush_names() const;
        QStackedWidget* regions_stack() const;
        const brush_panel& brush_panel() const;
        std::tuple<std::vector<ch::brush_expr_ptr>, std::vector<double>>
            brush_per_intervals() const;
        bool has_layers() const;
        ch::ink_layers* layers();

    signals:
        void change_source_image(cv::Mat& img);

    public slots:
        void handle_source_image_change(cv::Mat& img);
        void handle_view_scale_change(double scale);
        void handle_view_bw_change(bool bw);

    private:
        void create_main_menu();
        QWidget* create_image_processing_pipeline_tools();
        QWidget* create_drawing_tools();
        QWidget* create_region_map_tools();

        void display(pipeline_output work_in_prog = {});
        void handle_pipeline_change(int index);
        pipeline_output input_to_nth_stage(int index) const;
        ch::crosshatching_job drawing_job() const;
        const ch::ink_layers* layers() const;
        std::vector<cv::Mat> layer_images() const;
        ch::parameters drawing_params() const;
        std::string image_src_filename() const;
        vector_graphics_ptr vector_output() const;
        ch::dimensions<int> dimensions() const;
        std::tuple<int, int> source_image_sz() const;
        ui::view_state view_state() const;
        void save_project(const std::string& fname) const;
        void open_project(const std::string& fname);
        std::string state_to_json() const;
        void json_to_state(const std::string& str);

        image_processing_tools img_proc_ctrls_;
        drawing_tools crosshatching_;
        rgn_map_tools rgn_map_;
    };

};