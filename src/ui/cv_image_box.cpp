#include "cv_image_box.h"
#include "../crosshatching/util.hpp"

ui::cv_image_box::cv_image_box(cv::Mat img)
    : QLabel(nullptr)
{
    if (!img.empty()) {
        set_image(img);
    }
}

void ui::cv_image_box::set_image(cv::Mat img) {
    this->setPixmap(QPixmap::fromImage(QImage(img.data, img.cols, img.rows, img.step, QImage::Format_BGR888)));
}

void ui::cv_image_box::set_image(cv::Mat img, float scale) {
    auto mat = ch::scale(img, scale);
    set_image(mat);
}

ui::cv_image_box::~cv_image_box()
{
}
