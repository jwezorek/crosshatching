#pragma once

#include <QLabel>
#include <opencv2/core.hpp>

namespace ui {
    class image_box : public QLabel
    {
        Q_OBJECT

    public:
        image_box(cv::Mat img = {});
        void set_image(cv::Mat img);
        void set_image(QImage img);
        void set_image(cv::Mat img, float scale);
        ~image_box();
    };
}