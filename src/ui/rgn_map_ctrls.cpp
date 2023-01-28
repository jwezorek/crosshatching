#include "rgn_map_ctrls.hpp"
#include <QtWidgets>
#include <range/v3/all.hpp>
#include "../crosshatching/util.hpp"
#include <numbers>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

ui::flow_direction_panel::flow_direction_panel() {
    set_image(ch::blank_monochrome_bitmap(100));
}

ui::rgn_properties_panel::rgn_properties_panel() {
    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("brush"));
    layout->addWidget(new QComboBox());
    layout->addWidget(new QLabel("flow"));
    layout->addWidget(new flow_direction_panel());
    layout->addStretch();
}

ui::rgn_map_ctrl::rgn_map_ctrl() :
        base_sz_(-1),
        scale_(1.0),
        selecting_(false),
        cursor_radius_(1)
{
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);
}

void ui::rgn_map_ctrl::set_regions(vector_graphics_ptr gfx) {
    auto n = gfx->polygons.size();
    scaled_regions_.clear();
    colors_.clear();
    scaled_regions_.reserve(n);
    colors_.reserve(n);

    for (const auto& [color, poly] : gfx->polygons) {
        scaled_regions_.push_back(
            (scale_ != 1.0) ? ch::scale(poly, scale_) : poly
        );
        colors_.push_back(ch::color_to_monochrome(color));
    }

    tree_ = {};
    for (const auto& [i, poly] : rv::enumerate(scaled_regions_)) {
        int index = static_cast<int>(i);
        tree_.insert({ index, &poly });
    }
    base_sz_ = gfx->sz;
    auto sz = scale_ * base_sz_;
    setFixedSize({ sz.wd, sz.hgt });
    display();
}

int ui::rgn_map_ctrl::poly_at_point(const ch::point& pt) const {
    auto pair = tree_.find(pt);
    if (pair) {
        return pair->first;
    }
    return -1;
}

ch::polygon ui::rgn_map_ctrl::cursor_poly(float r) {
    int n = 20;
    ch::polygon output;
    auto& verts = output.outer();
    auto angle = (2.0 * std::numbers::pi) / n;
    for (int i = 0; i < n; i++) {
        float theta = angle * static_cast<float>(i);
        verts.emplace_back( r * std::cos(theta), r * std::sin(theta));
    }
    return output;
}

void ui::rgn_map_ctrl::display() {
    if (scaled_regions_.empty()) {
        return;
    }
    auto sz = scale_ * base_sz_;
    cv::Mat mat(sz.hgt, sz.wd, CV_8UC3);
    mat.setTo(cv::Scalar(255, 255, 255));
    QImage img = ch::mat_to_qimage(mat, false);
    QPainter g(&img);
    g.setRenderHint(QPainter::Antialiasing, true);

    for (const auto& [index, poly] : rv::enumerate(scaled_regions_)) {
        uchar gray = colors_[index];
        ch::paint_polygon(g, poly, ch::rgb(gray, gray, gray));
        if (selected_.contains(index)) {
            g.setOpacity(0.5);
            ch::paint_polygon(g, poly, ch::color(255, 0, 0));
            g.setOpacity(1.0);
        }
        if (curr_loc_ && cursor_radius_ > 1) {
            static std::vector<ch::polygon> cursor_polys;
            if (cursor_polys.empty()) {
                cursor_polys = rv::iota(0, 11) |
                    rv::transform(
                        [](int index)->ch::polygon {
                            if (index == 0) {
                                index = 10;
                            }
                            double radius = (5.0 * index * index) / 10.0 + 5.0;
                            return cursor_poly(radius);
                        }
                ) | r::to_vector;
            }
            auto cursor_poly = ch::transform(
                cursor_polys[cursor_radius_], ch::translation_matrix(*curr_loc_)
            );
            ch::paint_polygon(g, cursor_poly, ch::color(0, 0, 0), false, 1);
        }
    }
    set_image(mat);
}

void ui::rgn_map_ctrl::set_scale(double sc) {
    scale_ = sc;
    display();
}

void ui::rgn_map_ctrl::keyPressEvent(QKeyEvent* event) {
    auto key = event->key();
    if (key >= Qt::Key_0 && key <= Qt::Key_9) {
        cursor_radius_ = (key != Qt::Key_0) ? static_cast<int>(key - Qt::Key_0) : 10;
        display();
    }
}

void ui::rgn_map_ctrl::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::MouseButton::LeftButton) {
        selecting_ = true;
    }
}

void ui::rgn_map_ctrl::mouseMoveEvent(QMouseEvent* event) {
    setFocus();
    bool repaint = false;
    auto qpt = event->position();
    auto new_loc = ch::point{ static_cast<float>(qpt.x()), static_cast<float>(qpt.y()) };
    if (!curr_loc_ || new_loc != *curr_loc_) {
        curr_loc_ = new_loc;
        repaint = true;
    }
    if ((event->buttons() & Qt::LeftButton) && selecting_ && curr_loc_) {
        auto poly_index = poly_at_point( *curr_loc_ );
        if (poly_index >= 0) {
            selected_.insert(poly_index);
            repaint = true;
        }
    }
    if (repaint) {
        display();
    }
}

void ui::rgn_map_ctrl::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::MouseButton::LeftButton) {
        selecting_ = false;
    }
}