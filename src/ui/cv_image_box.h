#pragma once

#include <QLabel>
#include <opencv2/core.hpp>

namespace ui {
    class cv_image_box : public QLabel
    {
        Q_OBJECT

    public:
        cv_image_box(cv::Mat img = {});
        void set_image(cv::Mat img);
        void set_image(cv::Mat img, float scale);
        ~cv_image_box();
    };
}