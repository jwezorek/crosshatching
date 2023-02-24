#include "rgn_map_ctrl.h"
#include "../crosshatching/util.hpp"
#include <numbers>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    std::vector<std::tuple<ch::point, ch::point>> transform(
        const std::vector<std::tuple<ch::point, ch::point>>& lines,
        const ch::matrix& mat) {
        return lines |
            rv::transform(
                [&mat](auto&& tup)->std::tuple<ch::point, ch::point> {
                    const auto& [a, b] = tup;
                    return {
                        ch::transform(a, mat),
                        ch::transform(b, mat)
                    };
                }
        ) | r::to_vector;
    }

    QPixmap createPatternedPixmap(qreal theta)
    {
        float dim = 256;
        static std::vector<std::tuple<ch::point, ch::point>> lines;

        if (lines.empty()) {
            auto big_dim = 2 * dim;
            for (auto y = -dim; y <= dim; y += 8) {
                lines.push_back(
                    { {-dim,y},{dim,y} }
                );
            }
        }
        auto rotated = transform(lines, ch::rotation_matrix(theta));
        rotated = transform(rotated, ch::translation_matrix(dim / 2, dim / 2));

        QPixmap pixmap(256, 256);

        QPainter painter(&pixmap);
        painter.setRenderHint(QPainter::Antialiasing, false);
        painter.setPen(QPen(Qt::black, 1));
        pixmap.fill(Qt::transparent);
        for (const auto& [u, v] : rotated) {
            painter.drawLine(u.x, u.y, v.x, v.y);
        }
        painter.end();

        return pixmap;
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

    QPointF pt_to_qpointf(const ch::point& pt) {
        return QPointF(
            pt.x,
            pt.y
        );
    }

    ch::color get_nth_color(int n) {
        static std::vector<ch::color> colors = {
            ch::rgb(180,0,0),
            ch::rgb(255,97,3),
            ch::rgb(255,255,0),
            ch::rgb(0,255,0),
            ch::rgb(0,0,255),
            ch::rgb(138,43,226),
            ch::rgb(255, 192, 203),
            ch::rgb(255, 160, 122),
            ch::rgb(255, 250, 205),
            ch::rgb(152, 251, 152),
            ch::rgb(173, 216, 230),
            ch::rgb(230, 230, 250)
        };
        if (n >= 0 && n < colors.size()) {
            return colors[n];
        }
        return ch::rgb(255, 255, 255);
    }

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
}


ui::rgn_selection::rgn_selection()
{
}

const std::unordered_set<int>& ui::rgn_selection::selected_ids() const {
    return selection_;
}

bool ui::rgn_selection::has_selection() const {
    return !selection_.empty();
}

void ui::rgn_selection::clear() {
    old_selection_ = std::move(selection_);
    selection_.clear();
}

void ui::rgn_selection::reselect() {
    selection_ = std::move(old_selection_);
    old_selection_.clear();
}

bool ui::rgn_selection::is_selected_id(int id) const {
    return selection_.contains(id);
}

void ui::rgn_selection::insert(std::span<int> ids) {
    selection_.insert(ids.begin(), ids.end());
}

void ui::rgn_selection::remove(std::span<int> ids) {
    for (auto id : ids) {
        selection_.erase(id);
    }
}

ui::rgn_map_ctrl::rgn_map_ctrl(std::span<ch::brush_expr_ptr> brushes) :
        base_sz_(-1),
        scale_(1.0),
        is_showing_brushes_(false),
        is_showing_flow_(false),
        selecting_(false),
        cursor_radius_(1),
        set_brushes_(brushes.size()),
        brush_info_(
            brushes |
            rv::transform(
                [](auto br_ptr)->std::pair<const ch::brush_expr*, brush_info> {
                    return { br_ptr.get(), {} };
                }
            ) | r::to< std::unordered_map<const ch::brush_expr*, brush_info>>()
        ) {
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);

    for (auto [i, br] : rv::enumerate(brushes)) {
        brush_info_[br.get()].color = get_nth_color(i);
    }
}

const ui::rgn_selection& ui::rgn_map_ctrl::current_selection() const {
    return selection_;
}

void ui::rgn_map_ctrl::set_scale(double scale) {
    scale_ = scale;
    auto n = layer_->size();
    scaled_regions_.clear();
    scaled_regions_.reserve(n);

    for (auto& ili : *layer_) {
        scaled_regions_.push_back(
            (scale_ != 1.0) ? ch::scale(ili.poly, scale_) : ili.poly
        );
    }

    tree_ = {};
    for (const auto& [i, poly] : rv::enumerate(scaled_regions_)) {
        int index = static_cast<int>(i);
        tree_.insert({ index, &poly });
    }
    auto scaled_sz = scale_ * base_sz_;
    setFixedSize({ scaled_sz.wd, scaled_sz.hgt });
    update();
}

void ui::rgn_map_ctrl::set(ch::brush_expr_ptr default_brush, double sc,
    const ch::dimensions<int>& sz, ch::ink_layer* layer) {
    def_brush_ = default_brush;
    layer_ = layer;
    auto n = layer->size();
    colors_.clear();
    colors_.reserve(n);

    for (auto& ili : *layer) {
        colors_.push_back(ch::ink_shade_to_color(ili.value));
    }
    base_sz_ = sz;
    set_scale(sc);
}


std::vector<int> ui::rgn_map_ctrl::polys_at_point(const ch::point& pt) const {
    return (cursor_radius_ == 1) ?
        find_intersecting_polys(pt) :
        find_intersecting_polys(cursor_poly(scale_, pt, cursor_radius_));
}

ch::polygon ui::rgn_map_ctrl::cursor_poly(
    double scale, const ch::point& pt, int sz_index) const {
    static std::unordered_map<int, std::vector<ch::polygon>> cursor_polys;
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

void ui::rgn_map_ctrl::show_brushes(bool v) {
    is_showing_brushes_ = v;
}

void ui::rgn_map_ctrl::show_flow(bool v) {
    is_showing_flow_ = v;
}

bool ui::rgn_map_ctrl::is_showing_brushes() const {
    return is_showing_brushes_;
}

bool ui::rgn_map_ctrl::is_showing_flow() const {
    return is_showing_flow_;
}

ch::ink_layer* ui::rgn_map_ctrl::layer() const {
    return layer_;
}

void ui::rgn_map_ctrl::set_brush_of_selection(ch::brush_expr_ptr br) {
    if (!selection_.has_selection()) {
        return;
    }
    for (auto [index, ili] : rv::enumerate(selection_.selected_layer_items(layer_))) {
        brush_info_[ili->brush.get()].items.erase(index);
        ili->brush = br;
        brush_info_[br.get()].items.insert(index);
    }
    selection_.clear();
    emit selection_changed();
    update();
}


void ui::rgn_map_ctrl::set_flow_of_selection(double flow) {
    if (!selection_.has_selection()) {
        return;
    }
    for (auto [index, ili] : rv::enumerate(selection_.selected_layer_items(layer_))) {
        ili->flow_dir = flow;
    }
    selection_.clear();
    emit selection_changed();
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
        const auto& poly = scaled_regions_[poly_index];
        ch::paint_polygon(painter, poly, colors_[poly_index]);
        if (selection_.is_selected_id(poly_index)) {
            painter.setOpacity(0.5);
            ch::paint_polygon(painter, poly, ch::color(255, 0, 0));
            painter.setOpacity(1.0);
        }
        auto brush = layer_->at(poly_index).brush;
        if (brush != def_brush_ && is_showing_brushes_) {
            painter.setOpacity(0.3);
            ch::paint_polygon(painter, poly, brush_info_[brush.get()].color);
            painter.setOpacity(1.0);
        }
        auto flow = layer_->at(poly_index).flow_dir;
        if (flow != 0.0 && is_showing_flow_) {
            int key = static_cast<int>(flow * 1000);
            if (!patterns_.contains(key)) {
                patterns_[key] = createPatternedPixmap(flow);
            }
            QBrush patternBrush(patterns_[key]);
            patternBrush.setStyle(Qt::TexturePattern);
            painter.setCompositionMode(QPainter::CompositionMode_DestinationIn);
            ch::fill_polygon(painter, poly, patternBrush);
            painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
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

    if (flow_pt1_ && flow_pt2_) {
        painter.setPen(Qt::black);
        painter.drawLine(pt_to_qpointf(*flow_pt1_), pt_to_qpointf(*flow_pt2_));
    }
}

void ui::rgn_map_ctrl::keyPressEvent(QKeyEvent* event) {
    auto key = event->key();
    if (key >= Qt::Key_0 && key <= Qt::Key_9) {
        cursor_radius_ = (key != Qt::Key_0) ? static_cast<int>(key - Qt::Key_0) : 10;
        update();
    }
    if ((key == Qt::Key_R) && (event->modifiers() & Qt::ControlModifier)) {
        selection_.reselect();
        emit selection_changed();
        update();
    }
}

void ui::rgn_map_ctrl::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::MouseButton::LeftButton) {
        selecting_ = true;
        mouseMoveEvent(event);
    }
}

void ui::rgn_map_ctrl::handle_flow_drag(bool is_flow_dragging, const ch::point& pt) {

    if (!is_flow_dragging) {
        auto y_diff = pt.y - flow_pt1_->y;
        auto x_diff = pt.x - flow_pt1_->x;
        flow_pt1_ = flow_pt2_ = {};
        emit flow_assigned(std::atan2(y_diff, x_diff));
        update();
        return;
    }
    if (!flow_pt1_) {
        flow_pt1_ = pt;
        update();
        return;
    }
    flow_pt2_ = pt;
    update();
}

void ui::rgn_map_ctrl::mouseMoveEvent(QMouseEvent* event) {
    setFocus();
    bool selection_state_changed = false;
    bool repaint = false;
    auto qpt = event->position();
    auto new_loc = ch::point{ static_cast<float>(qpt.x()), static_cast<float>(qpt.y()) };
    if (!curr_loc_ || new_loc != *curr_loc_) {
        curr_loc_ = new_loc;
        repaint = true;
    }

    bool left_mouse_btn_down = event->buttons() & Qt::LeftButton;
    auto mods = QApplication::keyboardModifiers();
    bool is_flow_dragging = left_mouse_btn_down && (mods & Qt::ControlModifier);
    if (is_flow_dragging || flow_pt1_) {
        handle_flow_drag(is_flow_dragging, new_loc);
        return;
    }

    if (left_mouse_btn_down && selecting_ && curr_loc_) {
        auto poly_indices = polys_at_point(*curr_loc_);
        if (!poly_indices.empty()) {
            auto debug_layer_item = layer_->at(poly_indices.front());
            bool removing_from_selection = mods & Qt::AltModifier;
            bool adding_to_selection = (!removing_from_selection) && (mods & Qt::ShiftModifier);

            if (!adding_to_selection && !removing_from_selection) {
                selection_.clear();
                selection_.insert(poly_indices);
            } else if (adding_to_selection) {
                selection_.insert(poly_indices);
            } else if (removing_from_selection) {
                selection_.remove(poly_indices);
            }

            selection_state_changed = true;
            repaint = true;
        }
    }
    if (selection_state_changed) {
        emit selection_changed();
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

