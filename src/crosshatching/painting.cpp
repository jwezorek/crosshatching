#include "painting.hpp"

namespace {

}

cv::Mat ch::qimage_to_mat(QImage img) {
    return {};
}

QImage ch::mat_to_qimage(cv::Mat mat) {
    return {};
}

void ch::paint_polygon(QPainter& g, const polygon& poly, color col) {

}

void ch::paint_strokes(QPainter& g, strokes str) {

}

std::tuple<uchar, ch::polygon> ch::raster_to_vector_grayscale(cv::Mat mat, double param) {
    return {};
}

std::tuple<ch::color, ch::polygon> ch::raster_to_vector(cv::Mat mat, double param) {
    return {};
}