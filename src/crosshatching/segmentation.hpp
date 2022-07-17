#pragma once

#include "util.hpp"
#include <opencv2/core.hpp>
#include <unordered_set>
#include <vector>

namespace ch {

    namespace detail {
        class edge_set;
    }

    class segmentation {
    private:
        
        struct connected_component;

        cv::Mat label_mat_;
        std::vector<connected_component> connected_components_;

        void get_neighbors(int x, int y, std::array<int, 4>& ary);
        void insert_edges(int label, int x, int y, detail::edge_set& edges);
        void init_new_connected_component(int label, uchar val, int x, int y, detail::edge_set& edges);
        void update_connected_component(int label, int x, int y, detail::edge_set& edges);
        bool is_uninitialized(int index);

    public:
        segmentation(cv::Mat img, cv::Mat label_image);
        segmentation(cv::Mat mask);
        ~segmentation();

        uchar value(int index) const;
        std::vector<int> neighbor_indices(int index) const;
        cv::Rect bounds(int index) const;
        cv::Size dimensions() const;
        int count() const;

        template<typename T>
        void paint_component(int index, cv::Mat dst, T v) const {
            auto cc_bounds = this->bounds(index);
            cv::Mat roi = cv::Mat(label_mat_, cc_bounds);
            cv::Mat mask;
            cv::inRange(roi, index, index, mask);
            cv::Mat dst_roi(dst, cc_bounds);
            dst_roi.setTo(v, mask);
        }

        void paint_component(int index, cv::Mat dst) const;
        cv::Point2i component_to_point(int index) const;
        int point_to_component(const cv::Point2i& pt) const;

        struct cc_properties {
            int index;
            uchar color;
            int area;
        };

        cc_properties properties(int index) const;
        cv::Mat paint_components() const;
    };

    cv::Mat grayscale_to_label_image(cv::Mat input);

}