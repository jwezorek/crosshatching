#include "rgn_tool_panel.hpp"
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
/*
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
*/

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
