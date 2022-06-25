#pragma once

#include "../crosshatching/drawing.hpp"
#include "settingctrls.hpp"
#include "treepanel.h"
#include "cv_image_box.h"
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

    struct image_processing_ctrls {
        std::string src_filename;
        std::vector<image_processing_pipeline_item*> pipeline;
        cv_image_box* image_box;
        cv::Mat src;
        cv::Mat current;
        view_state view_state;
    };

    struct crosshatching_ctrls {
        tree_panel* brushes;
        list_panel* layers;
        cv_image_box* img_swatch;
        cv_image_box* drawing_swatch;
        cv::Mat swatch;
    };

    class crosshatching : public QMainWindow
    {
        Q_OBJECT

    public:
        crosshatching(QWidget* parent = Q_NULLPTR);

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
        void createMainMenu();
        QWidget* createImageProcPipelineCtrls();
        QWidget* createCrosshatchCtrls();
        void display(cv::Mat mat = {});
        void handle_pipeline_change(int index);
        cv::Mat input_to_nth_stage(int index) const;
        cv::Mat segmentation() const;
        ch::crosshatching_job drawing_job() const;
        std::vector<std::tuple<ch::brush_fn, double>> layers() const;
        ch::crosshatching_params drawing_params() const;
        std::string image_src_filename() const;
        cv::Mat processed_image() const;

        std::tuple<int, int> source_image_sz() const;
        ui::view_state view_state() const;

        image_processing_ctrls img_proc_ctrls_;
        crosshatching_ctrls crosshatching_;
    };

};