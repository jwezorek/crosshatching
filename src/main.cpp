#include "geometry.hpp"
#include "crosshatch.hpp"
#include "brush.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <iostream>
#include <functional>
#include <tuple>
#include <vector>
#include <random>
#include <numbers>


int main()
{
    auto brush = ch::make_linear_hatching_brush_fn(
        ch::make_lerped_normal_dist_fn(10, 25, 100, 20),
        ch::make_lerped_normal_dist_fn(100, 20, 0, 0.05),
        ch::make_lerped_normal_dist_fn(7, 0.5, 0.75, 0.05),
        ch::make_default_hatching_unit()
    );

    int n = 100;
    for (int i = 0; i <= n; i++) {
        auto mat = ch::paint_cross_hatching(2, brush( i * (1.0/n)));
        cv::imshow("crosshatch", mat);
        int k = cv::waitKey(0);
    }

    return 0;
}
