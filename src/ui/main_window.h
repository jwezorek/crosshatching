#pragma once

#include "../crosshatching/drawing.hpp"
#include "../crosshatching/brush_language.hpp"
#include "settingctrls.hpp"
#include "treepanel.h"
#include "image_box.h"
#include "image_tab_ctrl.h"
#include <QtWidgets>
#include <QtWidgets/QMainWindow>
#include <opencv2/core.hpp>
#include <tuple>

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
        cv::Mat current;
        view_state view_state;
    };

	class layer_panel : public ui::list_panel {

        Q_OBJECT

	public:
		layer_panel();
		void set_brush_names(const std::vector<std::string>& brush_names);
        std::vector<std::tuple<std::string, double>> layers() const;

    signals:
        void layers_changed();

	private:

		std::vector<std::string> brush_names_;
		std::map<double, std::string> layers_;

        std::tuple<std::string, double> row(int n) const;
        void add_layer();
        void cell_double_clicked(int r, int col);
        void delete_layer();
        void insert_layer(const std::string& brush, double end_of_range);
        void setRowText(int row, const std::string& brush, const std::string& from, const std::string& to);
        void sync_layers_to_ui();
	};

	class brush_panel : public ui::tree_panel {

        Q_OBJECT

	public:

		brush_panel(layer_panel& layers);
		std::vector<std::string> brush_names() const;
        std::unordered_map<std::string, ch::brush_fn> brush_dictionary() const;

	private:

		layer_panel& layer_panel_;

		class brush_item : public QTreeWidgetItem {
		public:
            brush_item(const std::string& name, ch::brush_expr_ptr expr);
            brush_item(ch::brush_expr_ptr expr);
			ch::brush_expr_ptr brush_expr;
		};

        static void insert_brush_item(brush_item* parent, brush_item* item);
        static void insert_toplevel_item(QTreeWidget* tree, const std::string& name, ch::brush_expr_ptr expr);
        void add_brush_node();
        void delete_brush_node();
        void sync_layer_panel();
	};

    struct drawing_tools {

        enum class view : int {
            swatch = 0,
            layers = 1,
            drawing = 2
        };

        brush_panel* brushes;
        layer_panel* layers;
        image_box* img_swatch;
        image_box* drawing_swatch;
        image_box* drawing;
        image_tab_ctrl* layer_viewer;
        QStackedWidget* viewer_stack;
        cv::Mat swatch;

        void set_view(view v);
    };

    class main_window : public QMainWindow
    {
        Q_OBJECT

    public:
        main_window(QWidget* parent = Q_NULLPTR);

        void open();
        void save_processed_image();
        void generate();
        void test();
        void redo_test_swatch();
        void edit_settings();

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

        void set_swatch_view(cv::Mat swatch, bool left);
        void set_layer_view();
        void set_drawing_view(cv::Mat drawing);

        void display(cv::Mat mat = {});
        void handle_pipeline_change(int index);
        cv::Mat input_to_nth_stage(int index) const;
        cv::Mat segmentation() const;
        ch::crosshatching_job drawing_job() const;
        std::vector<std::tuple<ch::brush_fn, double>> brush_per_intervals() const;
        std::vector<std::tuple<ch::brush_fn, cv::Mat>> layers() const;
        std::vector<cv::Mat> layer_images() const;
        ch::crosshatching_params drawing_params() const;
        std::string image_src_filename() const;
        cv::Mat processed_image() const;

        std::tuple<int, int> source_image_sz() const;
        ui::view_state view_state() const;

        image_processing_tools img_proc_ctrls_;
        drawing_tools crosshatching_;
    };

};