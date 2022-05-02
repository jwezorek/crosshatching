#pragma once

#include "settingctrls.hpp"
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

    class crosshatching : public QMainWindow
    {
        Q_OBJECT

    public:
        crosshatching(QWidget* parent = Q_NULLPTR);

        void open();
        void generate();
        void debug();

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

        std::tuple<int, int> source_image_sz() const;
        ui::view_state view_state() const;

        std::vector<image_processing_pipeline_item*> imgproc_pipeline_;

        QLabel* image_box_;
        cv::Mat src_image_;
        cv::Mat current_image_;
        ui::view_state view_state_;
    };

};