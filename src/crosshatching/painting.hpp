#pragma once

#include "geometry.hpp"
#include "util.hpp"
#include <QPainter>
#include <QImage>
#include <opencv2/core.hpp>
#include <vector>
#include <tuple>

namespace ch {

    cv::Mat qimage_to_mat(QImage img);
    QImage mat_to_qimage(cv::Mat mat);
    void paint_polygon(QPainter& g, const polygon& poly, color col);
    void paint_strokes(QPainter& g, strokes str);

    std::tuple<uchar, polygon> raster_to_vector_grayscale(cv::Mat mat, double param);
    std::tuple<color, polygon> raster_to_vector(cv::Mat mat, double param);
}