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

ui::rgn_tool_panel::rgn_tool_panel() {

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
    
}

void ui::rgn_tool_panel::set(int num_layers, const std::vector<std::string>& brushes) {
    set_layers(num_layers);
    brush_lbl_->setText(k_no_selection);
    select_brush_btn_->set_items(brushes);
}

void ui::rgn_tool_panel::set_layers(int n) {
    layer_cbo_->clear();
    if (n == 0) {
        return;
    }
    for (int i = 0; i < n; ++i) {
        std::string lbl = "layer " + std::to_string(i);
        layer_cbo_->insertItem(i, lbl.c_str());
    }
}

void ui::rgn_tool_panel::set_brush_name(const std::string& brush_name) {
    if (brush_name.empty()) {
        brush_lbl_->setText(k_no_selection);
    } else {
        brush_lbl_->setText(std::string("brush: " + brush_name).c_str());
    }
}

void ui::rgn_tool_panel::set_hetero_brush() {
    brush_lbl_->setText(k_heterogeneous);
}

void ui::rgn_tool_panel::set_flow(std::optional<float> flow) {
    if (flow) {
        flow_ctrl_->set_direction(*flow);
    } else {
        flow_ctrl_->clear();
    }
}

void ui::rgn_tool_panel::connect_to_tools(rgn_map_tools* tools) {

    connect(layer_cbo_, &QComboBox::currentIndexChanged,
        tools->rgn_map_stack(), &QStackedWidget::setCurrentIndex);

    connect(brush_cbx_, &QCheckBox::stateChanged,
        [tools](bool checked) {
            auto rm = tools->current_rgn_map();
            rm->show_brushes(checked);
            rm->update();
        }
    );
    connect(flow_cbx_, &QCheckBox::stateChanged,
        [tools](bool checked) {
            auto rm = tools->current_rgn_map();
            rm->show_flow(checked);
            rm->update();
        }
    );
    connect(select_brush_btn_, &select_button::item_clicked,
        [tools](QString str) {tools->handle_brush_assigned(str); });
}

/*------------------------------------------------------------------------------------------------*/

void ui::rgn_tool_panel::handle_brush_name_change(std::string old_name, std::string new_name) {
    qDebug() << "brush name change";
}

ui::rgn_map_tools::rgn_map_tools(main_window* parent) : 
    parent_(parent),
    rgn_map_stack_(nullptr), 
    rgn_props_(nullptr) {

}

bool ui::rgn_map_tools::has_rgn_maps() const
{
    if (!rgn_map_stack_) {
        return false;
    }
    auto children = rgn_map_stack_->findChildren<QWidget*>(QString(), Qt::FindDirectChildrenOnly);
    return !children.empty();
}

void ui::rgn_map_tools::clear() {
    if (!rgn_map_stack_) {
        return;
    }
    auto children = rgn_map_stack_->findChildren<QWidget*>("",Qt::FindDirectChildrenOnly);
    for (auto child : children) {
        rgn_map_stack_->removeWidget(child);
        delete child;
    }
}

void ui::rgn_map_tools::populate(ui::main_window* parent) {
    rgn_map_stack_ = new QStackedWidget();
    rgn_props_ = new rgn_tool_panel();
    rgn_props_->connect_to_tools(this);
    parent->connect(&(parent->brush_panel()), &brush_panel::brush_name_changed,
        rgn_props_, &rgn_tool_panel::handle_brush_name_change);
    rgn_props_->set(0, parent->brush_names());
}

void ui::rgn_map_tools::set_layers(double scale, ch::ink_layers* layers) {
    rgn_props_->set_layers(layers->content.size());

    int old_num_layers = rgn_map_stack_->count();
    if (old_num_layers > 0) {
        for (int i = old_num_layers - 1; i >= 0; --i) {
            auto rgn_map = rgn_map_stack_->widget(i);
            rgn_map_stack_->removeWidget(rgn_map);
            delete rgn_map;
        }
    }
    auto def_brushes = get_brush_defaults();
    auto brushes = parent_->brush_panel().brushes();
    for (int i = 0; i < layers->content.size(); ++i) {
        auto rgn_map = new rgn_map_ctrl(brushes);
        rgn_map->set(def_brushes[i], scale, layers->sz, &(layers->content[i]));
        auto scroller = new QScrollArea();
        scroller->setWidget(rgn_map);
        rgn_map_stack_->addWidget(scroller);
        rgn_map->connect(
            rgn_map, &rgn_map_ctrl::selection_changed,
            [this]() {handle_selection_change(); }
        );
        rgn_map->connect(
            rgn_map, &rgn_map_ctrl::flow_assigned,
            [this](double theta) { handle_flow_assigned(theta); }
        );
        rgn_map->setVisible(true);
    }
    int test = rgn_map_stack_->count();
    rgn_map_stack_->setCurrentIndex(0);
    rgn_map_stack_->setVisible(true);
}

QStackedWidget* ui::rgn_map_tools::rgn_map_stack() const {
    return rgn_map_stack_;
}

ui::rgn_tool_panel* ui::rgn_map_tools::rgn_props() const {
    return rgn_props_;
}

std::vector<ch::brush_expr_ptr> ui::rgn_map_tools::get_brush_defaults() const {
    return std::get<0>(parent_->brush_per_intervals());
}

void  ui::rgn_map_tools::set_scale(double scale) {
    auto rgn_maps = rgn_map_stack_->findChildren<rgn_map_ctrl*>();
    for (auto rm : rgn_maps) {
        rm->set_scale(scale);
    }
}

ui::rgn_map_ctrl* ui::rgn_map_tools::current_rgn_map() const {
    return static_cast<rgn_map_ctrl*>(
        static_cast<QScrollArea*>(rgn_map_stack_->currentWidget())->widget()
    );
}

std::unordered_map<ch::brush_expr*, std::string> ui::rgn_map_tools::brush_to_name_dict() const {
    auto name_to_brush_ = parent_->brush_panel().brush_dictionary();
    return name_to_brush_ |
        rv::transform(
            [](auto&& v)->std::unordered_map<ch::brush_expr*, std::string>::value_type {
                auto [name, br] = v;
                return { br.get(), name };
            }
    ) | r::to<std::unordered_map<ch::brush_expr*, std::string>>();
}

void ui::rgn_map_tools::handle_brush_assigned(QString qbrush_name) {
    auto name_to_brush_ = parent_->brush_panel().brush_dictionary();
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
    rgn_props_->set_brush_name(brush_name);
}

void ui::rgn_map_tools::handle_flow_assigned(double theta) {
    current_rgn_map()->set_flow_of_selection(theta);
}

void ui::rgn_map_tools::handle_selection_change() {
    handle_selection_change_brush();
    handle_selection_change_flow();
}

void ui::rgn_map_tools::handle_selection_change_brush() {
    std::optional<ch::brush_expr_ptr> selected_brush = {};

    auto rgn_map = current_rgn_map();
    const auto& sel = rgn_map->current_selection();

    if (!sel.has_selection()) {
        rgn_props_->set_brush_name();
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
        rgn_props_->set_hetero_brush();
        return;
    }
    auto dict = brush_to_name_dict();
    auto brush_name = dict[selected_brush->get()];
    rgn_props_->set_brush_name(brush_name);
}

void ui::rgn_map_tools::handle_selection_change_flow() {
    std::optional<float> selected_flow = {};

    auto rgn_map = current_rgn_map();
    const auto& sel = rgn_map->current_selection();

    if (!sel.has_selection()) {
        rgn_props_->set_flow({});
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
    rgn_props_->set_flow(selected_flow);
}