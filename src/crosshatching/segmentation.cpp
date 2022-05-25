#include "segmentation.hpp"
#include "util.hpp"
#include <boost/functional/hash.hpp>
#include <unordered_set>
#include <tuple>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

namespace r = ranges;
namespace rv = ranges::views;

/*--------------------------------------------------------------------------------------------------------------*/
// edge_set

class ch::segmentation::edge_set {
private:

    struct hasher {
        size_t operator()(const std::tuple<int, int>& edge) const {
            const auto& [i, j] = edge;
            std::size_t seed = 0;
            boost::hash_combine(seed, i);
            boost::hash_combine(seed, j);
            return seed;
        }
    };

    std::unordered_set<std::tuple<int, int>, hasher> set_;

public:
    void insert(int i, int j) {
        set_.insert(std::tuple<int, int>(i, j));
    }

    bool contains(int i, int j) {
        return set_.find(std::tuple<int, int>(i, j)) != set_.end();
    }
};

struct ch::segmentation::connected_component {
    int index;
    uchar value;
    cv::Rect bounds;
    cv::Point2i point;
    int area;
    std::vector<connected_component*> neighbors;

    connected_component(int i = -1, cv::Rect b = {}, int a = 0, uchar v = 255, ch::point p = { -1,-1 }) :
        index(i), point(p), area(a), value(v)
    {}
};

/*--------------------------------------------------------------------------------------------------------------*/
// segmentation

void ch::segmentation::get_neighbors(int x, int y, std::array<int, 4>& ary) {
    int label = label_mat_.at<int>(y, x);
    std::fill(ary.begin(), ary.end(), -1);
    int i = 0;
    int neigh;
    if ((y > 0) && ((neigh = label_mat_.at<int>(y - 1, x)) != label)) {
        ary[i++] = neigh;
    }
    if ((y < label_mat_.rows - 1) && ((neigh = label_mat_.at<int>(y + 1, x)) != label)) {
        ary[i++] = neigh;
    }
    if ((x > 0) && ((neigh = label_mat_.at<int>(y, x - 1)) != label)) {
        ary[i++] = neigh;
    }
    if ((x < label_mat_.cols - 1) && ((neigh = label_mat_.at<int>(y, x + 1)) != label)) {
        ary[i++] = neigh;
    }
}

void ch::segmentation::insert_edges(int label, int x, int y, edge_set& edges) {
    std::array<int, 4> neighbors;
    get_neighbors(x, y, neighbors);
    for (auto neigh : neighbors) {
        if (neigh == -1) {
            break;
        }
        if (edges.contains(label, neigh)) {
            continue;
        }
        edges.insert(label, neigh);
        connected_components_[label].neighbors.push_back(&connected_components_[neigh]);
    }
}

void ch::segmentation::init_new_connected_component(int label, uchar val, int x, int y, edge_set& edges) {
    auto& c = connected_components_[label];
    c.index = label;
    c.value = val;
    c.bounds = cv::Rect(x, y, 1, 1);
    c.area = 1;
    c.point = { x,y };
    insert_edges(label, x, y, edges);
}

void ch::segmentation::update_connected_component(int label, int x, int y, edge_set& edges) {
    auto& c = connected_components_[label];
    c.index = label;
    c.area++;
    c.bounds = union_rect_and_pt(c.bounds, { x,y });
    ch::segmentation::insert_edges(label, x, y, edges);
}

bool ch::segmentation::is_uninitialized(int index) {
    return connected_components_[index].index == -1;
}

ch::segmentation::segmentation(cv::Mat img, cv::Mat label_image) :
    label_mat_(label_image) {
    edge_set edges;
    int n = max_val_in_mat(label_image) + 1;
    connected_components_.resize(n);
    for (int y = 0; y < img.rows; y++) {
        for (int x = 0; x < img.cols; x++) {
            int label = label_image.at<int>(y, x);
            if (is_uninitialized(label)) {
                uchar value = img.at<uchar>(y, x);
                init_new_connected_component(label, value, x, y, edges);
            } else {
                update_connected_component(label, x, y, edges);
            }
        }
    }
}

ch::segmentation::segmentation(cv::Mat mask) {
    cv::Mat stats;
    cv::Mat centroids;
    cv::connectedComponentsWithStats(mask, label_mat_, stats, centroids);

    int n = stats.rows - 1;
    connected_components_.resize(n);
    for (int row = 1; row < n + 1; ++row) {
        int index = row - 1;
        auto& cc = connected_components_[index];
        cc.index = index;
        cc.bounds = cv::Rect(
            stats.at<int>(cv::Point(0, row)),
            stats.at<int>(cv::Point(1, row)),
            stats.at<int>(cv::Point(2, row)),
            stats.at<int>(cv::Point(3, row))
        );
        cc.area = stats.at<int>(cv::Point(4, row));
    }

    label_mat_.forEach<int>(
        [this, &mask](int& lbl, const int* position) -> void {
            lbl--;
            if (lbl >= 0) {
                int x = position[1];
                int y = position[0];
                auto& cc = this->connected_components_[lbl];
                if (cc.point.x < 0) {
                    cc.point = { x,y };
                    cc.value = mask.at<uchar>(y, x);
                }
            }
        }
    );
}

uchar ch::segmentation::value(int index) const {
    return connected_components_[index].value;
}

std::vector<int> ch::segmentation::neighbor_indices(int index) const {
    const auto& neighbors = connected_components_[index].neighbors;
    return neighbors |
        rv::transform(
            [](const auto* cc)->int {
                return cc->index;
            }
        ) |
        r::to_vector;
}

cv::Rect ch::segmentation::bounds(int index) const {
    return connected_components_[index].bounds;
}

cv::Size ch::segmentation::dimensions() const {
    return { label_mat_.cols, label_mat_.rows };
}

int ch::segmentation::count() const {
    return static_cast<int>(connected_components_.size());
}

void ch::segmentation::paint_component(int index, cv::Mat dst) const {
    const auto& c = connected_components_.at(index);
    paint_component(index, dst, c.value);
}

cv::Point2i ch::segmentation::component_to_point(int index) const {
    return connected_components_.at(index).point;
}

int ch::segmentation::point_to_component(const cv::Point2i& pt) const {
    return label_mat_.at<int>(pt.y, pt.x);
}

ch::segmentation::cc_properties ch::segmentation::properties(int index) const {
    const auto& cc = connected_components_.at(index);
    return {
        cc.index,
        cc.value,
        cc.area
    };
}

cv::Mat ch::segmentation::paint_components() const {
    cv::Mat img(dimensions(), CV_8UC1, cv::Scalar(255));
    for (int i = 0; i < count(); ++i) {
        paint_component(i, img);
    }
    return img;
}

ch::segmentation::~segmentation() = default;

/*--------------------------------------------------------------------------------------------------------------*/

cv::Mat ch::grayscale_to_label_image(cv::Mat input) {
    auto values = ch::unique_gray_values(input);
    cv::Mat output(input.size(), CV_32SC1, cv::Scalar(-1));
    int label = 0;
    for (auto val : values) {
        cv::Mat mask;
        cv::inRange(input, val, val, mask);
        ch::segmentation seg(mask);
        for (int i = 0; i < seg.count(); ++i) {
            seg.paint_component(i, output, label++);
        }
    }
    return output;
}