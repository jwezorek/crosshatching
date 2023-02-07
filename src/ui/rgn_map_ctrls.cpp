#include "rgn_map_ctrls.hpp"
#include "main_window.h"
#include <QtWidgets>
#include <range/v3/all.hpp>
#include "../crosshatching/util.hpp"
#include <numbers>
#include <qdebug.h>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    constexpr auto k_no_selection = "<no selection>";
    constexpr auto k_heterogeneous = "<heterogeneous selection>";

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

    auto get_region_maps(QStackedWidget* stack) {
        int n = stack->count();
        return rv::iota(0,n) |
            rv::transform(
                [stack](auto i)->ui::rgn_map_ctrl* {
                    auto scroller = static_cast<QScrollArea*>(stack->widget(i));
                    return static_cast<ui::rgn_map_ctrl*>( scroller->widget());
                }
            );
    }
}

ui::flow_direction_panel::flow_direction_panel() {
    set_image(ch::blank_monochrome_bitmap(100));
}

std::vector<ch::brush_expr_ptr> ui::rgn_map_panel::get_brush_defaults() const {
    return std::get<0>(parent_->brush_per_intervals());
}

ui::rgn_map_panel::rgn_map_panel(main_window* parent, QStackedWidget* stack) :
        parent_(parent) {

    auto layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("layer"));
    layout->addWidget(layer_cbo_ = new QComboBox());
    layout->addSpacerItem(new QSpacerItem(0, 30));
    layout->addWidget(brush_cbx_ = new QCheckBox("show set brushes"));
    layout->addWidget(flow_cbx_ = new QCheckBox("show set flow"));
    layout->addSpacerItem(new QSpacerItem(0, 30));

    auto box = new QGroupBox("properties");
    auto box_layout = new QVBoxLayout(box);
    box_layout->addWidget(new QLabel("brush"));
    box_layout->addWidget(curr_brush_cbo_ = new QComboBox());
    box_layout->addWidget(new QLabel("flow"));
    box_layout->addWidget(flow_ctrl_ = new flow_direction_panel());

    layout->addWidget(box);
    layout->addStretch();

    curr_brush_cbo_->setEditable(true);
    curr_brush_cbo_->lineEdit()->setReadOnly(true);
    stack_ = stack;

    connect(layer_cbo_, &QComboBox::currentIndexChanged, 
        stack_, &QStackedWidget::setCurrentIndex);
    connect(curr_brush_cbo_, &QComboBox::currentIndexChanged,
        this, &rgn_map_panel::handle_brush_change);
    connect(curr_brush_cbo_->lineEdit(), &QLineEdit::textChanged,
        [this](const QString& txt) {
            if (txt != k_no_selection && txt != k_heterogeneous) {
                int index = curr_brush_cbo_->findText(txt);
                if (index >= 0) {
                    handle_brush_change(index);
                }
            }
        }
    );
}

void ui::rgn_map_panel::repopulate_ctrls() {
    if (!parent_->has_layers()) {
        return;
    }
    int num_layers = static_cast<int>(parent_->layers()->count());
    layer_cbo_->clear();
    for (int i = 0; i < num_layers; ++i) {
        std::string lbl = "layer " + std::to_string(i);
        layer_cbo_->insertItem(i, lbl.c_str());
    }

    curr_brush_cbo_->clear();
    auto brushes = parent_->brush_names();
    for (const auto& [i,brush] : rv::enumerate(brushes)) {
        curr_brush_cbo_->insertItem(i, brush.c_str());
    }

    default_brushes_ = get_brush_defaults();
    name_to_brush_ = parent_->brush_panel().brush_dictionary();
    brush_to_name_ = name_to_brush_ |
        rv::transform(
            [](auto&& v)->std::unordered_map<ch::brush_expr*, std::string>::value_type {
                auto [name, br] = v;
                return { br.get(), name };
            }
        ) | r::to<std::unordered_map<ch::brush_expr*, std::string>>();

    curr_brush_cbo_->lineEdit()->setText(k_no_selection);
}

void ui::rgn_map_panel::set_layers(double scale, ch::ink_layers* ink_layers) {
    auto def_brushes = get_brush_defaults();
    int old_num_layers = stack_->count();
    if (old_num_layers > 0) {
        for (int i = old_num_layers - 1; i >= 0; --i) {
            auto rgn_map = stack_->widget(i);
            stack_->removeWidget(rgn_map);
            delete rgn_map;
        }
    }
    layer_cbo_->clear();

    if (!ink_layers || ink_layers->empty()) {
        return;
    }

    int n = static_cast<int>(ink_layers->content.size());
    for (int i = 0; i < n; ++i) {
        std::string lbl = "layer " + std::to_string(i);
        layer_cbo_->insertItem(i, lbl.c_str());
    }

    for (int i = 0; i < n; ++i) {
        auto rgn_map = new rgn_map_ctrl();
        rgn_map->set(def_brushes[i], scale, ink_layers->sz, &(ink_layers->content[i]));
        auto scroller = new QScrollArea();
        scroller->setWidget(rgn_map);
        stack_->addWidget(scroller);
        connect(
            rgn_map, &rgn_map_ctrl::selection_changed,
            this, &rgn_map_panel::handle_selection_change
        );
        rgn_map->setVisible(true);
    }
    int test = stack_->count();
    stack_->setCurrentIndex(0);
    stack_->setVisible(true);
}

ui::rgn_map_ctrl* ui::rgn_map_panel::current_rgn_map() const {
    return static_cast<rgn_map_ctrl*>(
        static_cast<QScrollArea*>(stack_->currentWidget())->widget()
    );
}

void ui::rgn_map_panel::handle_brush_change(int brush_index) {
    if (name_to_brush_.empty()) {
        return;
    }
    auto brush_name = curr_brush_cbo_->itemText(brush_index).toStdString();
    if (!name_to_brush_.contains(brush_name)) {
        qDebug() << "bad rgn_map_panel::handle_brush_change: " << brush_name.c_str();
        return;
    }
    auto brush = name_to_brush_[brush_name];
    for (auto ili_ptr : current_rgn_map()->selected_layer_items()) {
        ili_ptr->brush = brush;
    }
}

void ui::rgn_map_panel::handle_selection_change() {
    std::optional<ch::brush_expr_ptr> selected_brush = {};
    auto rgn_map = current_rgn_map();
    if (!rgn_map->has_selection()) {
        curr_brush_cbo_->lineEdit()->setText(k_no_selection);
        return;
    }
    auto selected = rgn_map->selected_layer_items() | r::to_vector;
    for ( auto sel : selected) {
        if (!selected_brush) {
            selected_brush = sel->brush;
        } else {
            if (selected_brush != sel->brush) {
                selected_brush = {};
                break;
            }
        }
    }
    if (!selected_brush) {
        curr_brush_cbo_->lineEdit()->setText(k_heterogeneous);
        return;
    }
    auto brush_name = brush_to_name_[selected_brush->get()];
    curr_brush_cbo_->setCurrentText(brush_name.c_str());
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

void ui::rgn_map_ctrl::set(ch::brush_expr_ptr default_brush, double scale, 
        const ch::dimensions<int>& sz, ch::ink_layer* layer) {
    def_brush_ = default_brush;
    scale_ = scale;
    layer_ = layer;
    auto n = layer->size();
    scaled_regions_.clear();
    colors_.clear();
    scaled_regions_.reserve(n);
    colors_.reserve(n);

    for ( auto& ili : *layer ) {
        scaled_regions_.push_back(
            (scale_ != 1.0) ? ch::scale(ili.poly, scale_) : ili.poly
        );
        colors_.push_back(ch::ink_shade_to_color(ili.value));
    }

    tree_ = {};
    for (const auto& [i, poly] : rv::enumerate(scaled_regions_)) {
        int index = static_cast<int>(i);
        tree_.insert({ index, &poly });
    }
    base_sz_ = sz;
    auto scaled_sz = scale_ * base_sz_;
    setFixedSize({ scaled_sz.wd, scaled_sz.hgt });
    update();
}

std::vector<int> ui::rgn_map_ctrl::polys_at_point(const ch::point& pt) const {
    return (cursor_radius_ == 1) ?
        find_intersecting_polys(pt) :
        find_intersecting_polys(cursor_poly(scale_,pt, cursor_radius_));
}

ch::polygon ui::rgn_map_ctrl::cursor_poly(
        double scale, const ch::point& pt, int sz_index) const {
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

const std::unordered_set<int>& ui::rgn_map_ctrl::selected() const {
    return selected_;
}

bool ui::rgn_map_ctrl::has_selection() const {
    return !selected_.empty();
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
    bool selection_state_changed = false;
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
                selection_state_changed = true;
            } else {
                for (auto removee : poly_indices) {
                    selected_.erase(removee);
                    selection_state_changed = true;
                }
            }
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

