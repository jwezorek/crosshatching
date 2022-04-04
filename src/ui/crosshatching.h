#pragma once

#include "settingctrls.hpp"
#include <QtWidgets>
#include <QtWidgets/QMainWindow>
#include <opencv2/core.hpp>
#include <tuple>

namespace ui {

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

    private:
        void createMainMenu();
        QWidget* createContent();

        void handle_scale_change(double new_scale);
        void handle_contrast_changed(std::tuple<double, double> params);

        preprocess_settings* preprocess_settings_;
        QLabel* image_box_;
        cv::Mat src_image_;
    };

};