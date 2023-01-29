#include "rgn_map_ctrls.hpp"
#include <QtWidgets>
#include <range/v3/all.hpp>
#include "../crosshatching/util.hpp"
#include <numbers>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    ch::polygon cursor_poly_aux(float r) {
        int n = 20;
        ch::polygon output;
        auto& verts = output.outer();
        auto angle = (2.0 * std::numbers::pi) / n;
        for (int i = 0; i < n; i++) {
            float theta = angle * static_cast<float>(i);
            verts.emplace_back(r * std::cos(theta), r * std::sin(theta));
        }
        return output;
    }

    ch::box qrect_to_box(const QRect& rect) {
        ch::point bottomLeft(rect.left(), rect.top());
        ch::point topRight(rect.right(), rect.bottom());
        return ch::box(bottomLeft, topRight);
    }

    ch::point qpoint_to_pt(const QPoint& pt) {
        return {
            static_cast<float>(pt.x()),
            static_cast<float>(pt.y())
        };
    }
}

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
    update();
}

std::vector<int> ui::rgn_map_ctrl::polys_at_point(const ch::point& pt) const {
    return (cursor_radius_ == 1) ?
        find_intersecting_polys(pt) :
        find_intersecting_polys(cursor_poly(scale_,pt, cursor_radius_));
}

ch::polygon ui::rgn_map_ctrl::cursor_poly(double scale, const ch::point& pt, int sz_index) const {
    static std::unordered_map<int,std::vector<ch::polygon>> cursor_polys;
    auto scale_key = static_cast<int>(scale * 100.0);
    if (cursor_polys[scale_key].empty()) {
        cursor_polys[scale_key] = rv::iota(0, 11) |
            rv::transform(
                [scale](int index)->ch::polygon {
                    if (index == 0) {
                        index = 10;
                    }
                    double radius = (5.0 * index * index) / 10.0 + 5.0;
                    return cursor_poly_aux(scale * radius);
                }
        ) | r::to_vector;
    }
    return ch::transform(
        cursor_polys[scale_key][cursor_radius_], ch::translation_matrix(*curr_loc_)
    );
}

void ui::rgn_map_ctrl::set_scale(double sc) {
    scale_ = sc;
    update();
}

void ui::rgn_map_ctrl::paintEvent(QPaintEvent* event) {
    auto dirty_rect = event->rect();

    QPainter painter(this);
    painter.fillRect(dirty_rect, Qt::white);
    painter.setRenderHint(QPainter::Antialiasing, true);
    if (scaled_regions_.empty()) {
        return;
    }

    auto visible_polys = find_intersecting_polys(qrect_to_box(dirty_rect));
    for (auto poly_index : visible_polys) {
        uchar gray = colors_[poly_index];
        const auto& poly = scaled_regions_[poly_index];
        ch::paint_polygon(painter, poly, ch::rgb(gray, gray, gray));
        if (selected_.contains(poly_index)) {
            painter.setOpacity(0.5);
            ch::paint_polygon(painter, poly, ch::color(255, 0, 0));
            painter.setOpacity(1.0);
        }
    }

    if (cursor_radius_ > 1 && rect().contains(mapFromGlobal(QCursor::pos()))) {
        auto pt = qpoint_to_pt(mapFromGlobal(QCursor::pos()));
        ch::paint_polygon(
            painter,
            cursor_poly(scale_, pt, cursor_radius_),
            ch::color(0, 0, 0),
            false
        );
    }

}

void ui::rgn_map_ctrl::keyPressEvent(QKeyEvent* event) {
    auto key = event->key();
    if (key >= Qt::Key_0 && key <= Qt::Key_9) {
        cursor_radius_ = (key != Qt::Key_0) ? static_cast<int>(key - Qt::Key_0) : 10;
        update();
    }
}

void ui::rgn_map_ctrl::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::MouseButton::LeftButton) {
        selecting_ = true;
        mouseMoveEvent(event);
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
        auto poly_indices = polys_at_point( *curr_loc_ );
        if (! poly_indices.empty()) {
            bool adding_to_selection = !(QApplication::keyboardModifiers() & Qt::AltModifier);
            if (adding_to_selection) {
                selected_.insert(poly_indices.begin(), poly_indices.end());
            } else {
                for (auto removee : poly_indices) {
                    selected_.erase(removee);
                }
            }
            repaint = true;
        }
    }
    if (repaint) {
        update();
    }
}

void ui::rgn_map_ctrl::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::MouseButton::LeftButton) {
        selecting_ = false;
        mouseMoveEvent(event);
    }
}