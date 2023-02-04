#pragma once

#include "../crosshatching/drawing.hpp"
#include "../crosshatching/brush_lang.hpp"
#include "../crosshatching/ink_layers.hpp"
#include "rgn_map_ctrls.hpp"
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

	class layer_panel : public ui::list_panel {

        Q_OBJECT

	public:
		layer_panel();
		void set_brush_names(const std::vector<std::string>& brush_names);
        std::vector<std::tuple<std::string, double>> layers() const;
        bool has_content() const;

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
        std::unordered_map<std::string, ch::brush_expr_ptr> brush_dictionary() const;

	private:

		layer_panel& layer_panel_;

		class brush_item : public QTreeWidgetItem {
		public:
            brush_item(const std::string& name, ch::brush_expr_ptr expr);
            brush_item(ch::brush_expr_ptr expr);
            
			ch::brush_expr_ptr brush_expression;
            bool is_toplevel;

            bool is_leaf() const;
		};

        static void insert_brush_item(brush_item* parent, brush_item* item);
        static void insert_toplevel_item(QTreeWidget* tree, const std::string& name, ch::brush_expr_ptr expr);
        void add_brush_node();
        void delete_brush_node();
        void sync_layer_panel();
        void handle_double_click(QTreeWidgetItem* item, int column);
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

    struct rgn_map_tools {
        rgn_map_ctrl* rgn_map_;
        rgn_map_panel* rgn_props_;
    };

    using layer_tuple = std::tuple<ch::brush_expr_ptr, std::vector<ch::gray_polygon>>;

    class main_window : public QMainWindow
    {
        friend class rgn_map_panel;

        Q_OBJECT

    public:
        main_window(QWidget* parent = Q_NULLPTR);

        void open();
        void save_processed_image();
        void generate();
        void test();
        void redo_test_swatch();
        void edit_settings();
        void debug();

        void set_swatch_view(cv::Mat swatch, bool left);
        void set_layer_view();
        void set_drawing_view(cv::Mat drawing);
        std::vector<std::string> brush_names() const;

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
        cv::Mat segmentation() const;
        ch::crosshatching_job drawing_job() const;
        std::tuple<std::vector<ch::brush_expr_ptr>, std::vector<double>> 
            brush_per_intervals() const;
        const ch::ink_layers* layers() const;
        ch::ink_layers* layers();
        bool has_layers() const;
        std::vector<cv::Mat> layer_images() const;
        ch::parameters drawing_params() const;
        std::string image_src_filename() const;
        vector_graphics_ptr vector_output() const;
        ch::dimensions<int> dimensions() const;
        std::tuple<int, int> source_image_sz() const;
        ui::view_state view_state() const;
        void tab_changed(int index);

        image_processing_tools img_proc_ctrls_;
        drawing_tools crosshatching_;
        rgn_map_tools rgn_map_;
    };

};