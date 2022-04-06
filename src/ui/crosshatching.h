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

    signals:
        void change_source_image(cv::Mat& img);

    public slots:
        void handle_source_image_change(cv::Mat& img);
        void handle_view_scale_change(double scale);
        void handle_view_bw_change(bool bw);

    private:
        void createMainMenu();
        QWidget* createContent();
        void display(cv::Mat mat = {});
        void preprocess_image(double scale, double beta, double sigma);
        cv::Mat preprocess_output() const;
        cv::Mat filter_output() const;

        void handle_scale_change(double new_scale);
        void handle_contrast_changed(std::tuple<double, double> params);
        void handle_shock_filter_changed(std::tuple<int, int, double, int> params);
        void handle_segmentation_changed(std::tuple<int, double, int> params);

        std::tuple<int, int> source_image_sz() const;
        ui::view_state view_state() const;

        preprocess_settings* preprocess_settings_;
        shock_filter_settings* shock_filter_settings_;
        segmentation_settings* segmentation_settings_;

        QLabel* image_box_;
        cv::Mat src_image_;
        cv::Mat current_image_;
        ui::view_state view_state_;
    };

};