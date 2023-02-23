#include "rgn_map_ctrls.hpp"
#include "main_window.h"
#include <QtWidgets>
#include <range/v3/all.hpp>
#include "../crosshatching/util.hpp"
#include <numbers>
#include <qdebug.h>
#include <fstream>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    constexpr auto k_no_selection = "<no selection>";
    constexpr auto k_heterogeneous = "<heterogeneous selection>";

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

/*------------------------------------------------------------------------------------------------*/

ui::flow_direction_panel::flow_direction_panel() : direction_(0) {
    setFixedSize(QSize(130, 130));
}

void ui::flow_direction_panel::set_direction(double theta) {
    direction_ = theta;
    update();
}

std::optional<double> ui::flow_direction_panel::direction() const {
    return direction_;
}

float ui::flow_direction_panel::radius() const {
    return static_cast<float>(width() - 2) / 2.0f;
}

QPointF ui::flow_direction_panel::center() const {
    auto r = static_cast<float>(width()) / 2.0f;
    return QPointF(r, r);
}

void ui::flow_direction_panel::clear() {
    direction_ = {};
    update();
}

void ui::flow_direction_panel::paintEvent(QPaintEvent* event) {
    QPainter g(this);
    g.setRenderHint(QPainter::Antialiasing, true);
    auto r = radius();
    auto c = center();
    g.setPen(Qt::white);
    g.drawEllipse(c, r, r);

    if (direction_) {
        auto theta = *direction_;
        g.setPen(QPen(Qt::white, 3.0f));
        QPointF end_pt = {
            r * std::cos(theta) + c.x(),
            r * std::sin(theta) + c.y()
        };

        g.drawLine(c, end_pt);
    }
}

/*------------------------------------------------------------------------------------------------*/

ui::select_button::select_button(const std::string& txt)
{
    setText(txt.c_str());
    connect(this, &QPushButton::clicked, this, &select_button::show_popup);

    list_ = new QListWidget;
    connect(list_, &QListWidget::itemClicked, this, &select_button::item_selected);

    popup_ = new QWidget;
    popup_->setWindowFlags(Qt::Popup);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(list_);
    popup_->setLayout(layout);
}

void ui::select_button::set_items(std::span<const std::string> items) {
    list_->clear();
    for (const auto& item_txt : items) {
        list_->addItem(item_txt.c_str());
    }
}

void ui::select_button::show_popup() {
    QPoint buttonPos = mapToGlobal(QPoint(0, height()));
    popup_->move(buttonPos);
    popup_->show();
}

void ui::select_button::item_selected(QListWidgetItem* item) {
    popup_->hide();
    emit item_clicked(item->text());
}

/*------------------------------------------------------------------------------------------------*/

ui::rgn_map_ctrl::rgn_map_ctrl(const std::vector<ch::brush_expr_ptr>* brushes) :
        base_sz_(-1),
        scale_(1.0),
        is_showing_brushes_(false),
        is_showing_flow_(false),
        selecting_(false),
        cursor_radius_(1),
        brushes_(brushes),
        set_brushes_(brushes->size()),
        brush_info_(
            *brushes |
               rv::transform(
                   [](auto br_ptr)->std::pair<const ch::brush_expr*, brush_info> {
                       return { br_ptr.get(), {} };
                   }
               ) | r::to< std::unordered_map<const ch::brush_expr*, brush_info>>()
        )   
{
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);

    for (auto [i,br] : rv::enumerate(*brushes)) {
        brush_info_[br.get()].color = get_nth_color(i);
    }
}

const ui::rgn_map_ctrl::selection& ui::rgn_map_ctrl::current_selection() const {
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

    for ( auto& ili : *layer ) {
        colors_.push_back(ch::ink_shade_to_color(ili.value));
    }
    base_sz_ = sz;
    set_scale(sc);
}

ui::rgn_map_ctrl::selection::selection()
{
}

const std::unordered_set<int>& ui::rgn_map_ctrl::selection::selected_ids() const {
    return selection_;
}

bool ui::rgn_map_ctrl::selection::has_selection() const {
    return !selection_.empty();
}

void ui::rgn_map_ctrl::selection::clear() {
    old_selection_ = std::move(selection_);
    selection_.clear();
}

bool ui::rgn_map_ctrl::selection::is_selected_id(int id) const {
    return selection_.contains(id);
}

void ui::rgn_map_ctrl::selection::insert(std::span<int> ids) {
    selection_.insert(ids.begin(), ids.end());
}

void ui::rgn_map_ctrl::selection::remove(std::span<int> ids) {
    for (auto id : ids) {
        selection_.erase(id);
    }
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
    if (! selection_.has_selection()) {
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
    if (! selection_.has_selection()) {
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
        auto poly_indices = polys_at_point( *curr_loc_ );
        if (! poly_indices.empty()) {
            auto debug_layer_item = layer_->at(poly_indices.front());
            bool removing_from_selection = mods & Qt::AltModifier;
            bool adding_to_selection = (!removing_from_selection) && (mods & Qt::ShiftModifier);

            if (!adding_to_selection && !removing_from_selection) {
                selection_.clear();
                selection_.insert(poly_indices);
            } else if (adding_to_selection) {
                selection_.insert(poly_indices);
            } else if (removing_from_selection){
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

/*------------------------------------------------------------------------------------------------*/

std::vector<ch::brush_expr_ptr> ui::rgn_tool_panel::get_brush_defaults() const {
    return std::get<0>(parent_->brush_per_intervals());
}

ui::rgn_tool_panel::rgn_tool_panel(main_window* parent, QStackedWidget* stack) :
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
    box_layout->addWidget(brush_lbl_ = new QLabel("brush: <>"));
    box_layout->addWidget(select_brush_btn_ = new select_button("select brush"));
    box_layout->addWidget(new QLabel("flow"));
    box_layout->addWidget(flow_ctrl_ = new flow_direction_panel());

    layout->addWidget(box);
    layout->addStretch();

    stack_ = stack;

    connect(layer_cbo_, &QComboBox::currentIndexChanged,
        stack_, &QStackedWidget::setCurrentIndex);
    connect(select_brush_btn_, &select_button::item_clicked,
        this, &rgn_tool_panel::handle_brush_change);
    connect(brush_cbx_, &QCheckBox::stateChanged,
        [this](bool checked) {
            auto rm = this->current_rgn_map();
            rm->show_brushes(checked);
            rm->update();
        }
    );
    connect(flow_cbx_, &QCheckBox::stateChanged,
        [this](bool checked) {
            auto rm = this->current_rgn_map();
            rm->show_flow(checked);
            rm->update();
        }
    );
}

void ui::rgn_tool_panel::repopulate_ctrls() {
    if (!parent_->has_layers()) {
        return;
    }
    int num_layers = static_cast<int>(parent_->layers()->count());
    layer_cbo_->clear();
    for (int i = 0; i < num_layers; ++i) {
        std::string lbl = "layer " + std::to_string(i);
        layer_cbo_->insertItem(i, lbl.c_str());
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

    brush_lbl_->setText(k_no_selection);
    select_brush_btn_->set_items(parent_->brush_names());

}

void layer_debug(const std::string& output_file, ch::ink_layer& layer, double scale) {

    auto [x1, y1, wd, hgt] = ch::bounding_rectangle(
        layer |
        rv::transform([](const auto& ili) {return ili.poly; }) |
        r::to_vector
    );

    std::ofstream outfile(output_file);
    outfile << ch::svg_header(static_cast<int>(scale * wd), static_cast<int>(scale * hgt));
    for (const auto& ili : layer) {
        outfile << ch::polygon_to_svg(ili.poly, ch::to_svg_color(get_nth_color(ili.id)), scale) << std::endl;
    }
    for (const auto& ili : layer) {
        ch::point pt = ch::representative_point(ch::scale(ili.poly, scale));
        auto x = pt.x;
        auto y = pt.y;
        std::string text = std::to_string(ili.id);
        outfile << "<text x=\"" << x << "\" y=\"" << y << "\" text-anchor=\"middle\">" << text << "</text>"
            << std::endl;
    }
    outfile << "</svg>" << std::endl;
    outfile.close();
}

void debug(ch::ink_layers& ink_layers) {
    if (ink_layers.content.size() != 2) {
        return;
    }
    auto& layers = ink_layers.content;
    for (int i = 0; i < layers.size(); ++i) {
        std::string fname = "C:\\test\\aaaa-debug-lay-" + std::to_string(i) + ".svg";
        layer_debug(fname, layers[i], 20.0);
    }
}

void ui::rgn_tool_panel::set_layers(double scale, ch::ink_layers* ink_layers) {
    //debug(*ink_layers);
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
    brushes_ = parent_->brush_panel().brushes();
    for (int i = 0; i < n; ++i) {
        auto rgn_map = new rgn_map_ctrl(&brushes_);
        rgn_map->set(def_brushes[i], scale, ink_layers->sz, &(ink_layers->content[i]));
        auto scroller = new QScrollArea();
        scroller->setWidget(rgn_map);
        stack_->addWidget(scroller);
        connect(
            rgn_map, &rgn_map_ctrl::selection_changed,
            this, &rgn_tool_panel::handle_selection_change
        );
        connect(
            rgn_map, &rgn_map_ctrl::flow_assigned,
            this, &rgn_tool_panel::handle_flow_assigned
        );
        rgn_map->setVisible(true);
    }
    int test = stack_->count();
    stack_->setCurrentIndex(0);
    stack_->setVisible(true);
}

void  ui::rgn_tool_panel::set_scale(double scale) {
    auto rgn_maps = stack_->findChildren<rgn_map_ctrl*>();
    for (auto rm : rgn_maps) {
        rm->set_scale(scale);
    }
}

/*------------------------------------------------------------------------------------------------*/

ui::rgn_map_ctrl* ui::rgn_tool_panel::current_rgn_map() const {
    return static_cast<rgn_map_ctrl*>(
        static_cast<QScrollArea*>(stack_->currentWidget())->widget()
    );
}

void ui::rgn_tool_panel::handle_brush_change(QString qbrush_name) {
    if (name_to_brush_.empty()) {
        return;
    }
    auto brush_name = qbrush_name.toStdString();
    if (!name_to_brush_.contains(brush_name)) {
        qDebug() << "bad rgn_map_panel::handle_brush_change: " << brush_name.c_str();
        return;
    }
    auto brush = name_to_brush_[brush_name];
    current_rgn_map()->set_brush_of_selection(brush);
    brush_lbl_->setText(std::string("brush: " + brush_name).c_str());
}


void ui::rgn_tool_panel::handle_flow_assigned(double theta) {
    current_rgn_map()->set_flow_of_selection(theta);
    //flow_ctrl_->set_direction(theta);
}

void ui::rgn_tool_panel::handle_selection_change_brush() {
    std::optional<ch::brush_expr_ptr> selected_brush = {};

    auto rgn_map = current_rgn_map();
    const auto& sel = rgn_map->current_selection();

    if (!sel.has_selection()) {
        brush_lbl_->setText(k_no_selection);
        return;
    }
    auto selected = sel.selected_layer_items(rgn_map->layer()) | r::to_vector;
    for (auto sel : selected) {
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
        brush_lbl_->setText(k_heterogeneous);
        return;
    }
    auto brush_name = brush_to_name_[selected_brush->get()];
    brush_lbl_->setText(("brush: " + brush_name).c_str());
}

void ui::rgn_tool_panel::handle_selection_change_flow() {
    std::optional<float> selected_flow = {};

    auto rgn_map = current_rgn_map();
    const auto& sel = rgn_map->current_selection();

    if (! sel.has_selection()) {
        flow_ctrl_->clear();
        return;
    }
    auto selected = sel.selected_layer_items(rgn_map->layer()) | r::to_vector;
    for (auto sel : selected) {
        if (!selected_flow) {
            selected_flow = sel->flow_dir;
        } else {
            if (selected_flow != sel->flow_dir) {
                selected_flow = {};
                break;
            }
        }
    }

    if (!selected_flow) {
        flow_ctrl_->clear();
        return;
    }
    flow_ctrl_->set_direction(*selected_flow);
}

void ui::rgn_tool_panel::handle_selection_change() {
    handle_selection_change_brush();
    handle_selection_change_flow();
}

ui::rgn_map_tools::rgn_map_tools() : rgn_map_stack(nullptr), rgn_props(nullptr) {

}

bool ui::rgn_map_tools::has_rgn_maps() const
{
    if (!rgn_map_stack) {
        return false;
    }
    auto children = rgn_map_stack->findChildren<QWidget*>(QString(), Qt::FindDirectChildrenOnly);
    return !children.empty();
}

void ui::rgn_map_tools::clear() {
    if (!rgn_map_stack) {
        return;
    }
    auto children = rgn_map_stack->findChildren<QWidget*>("",Qt::FindDirectChildrenOnly);
    for (auto child : children) {
        rgn_map_stack->removeWidget(child);
        delete child;
    }
}
